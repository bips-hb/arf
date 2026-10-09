# Shared helpers for the bench/ scripts (sourced, not part of the package build).
#
# Each backend measurement runs in its OWN fresh subprocess (callr spawns a new
# R process rather than forking). This is deliberate: mirai/nanonext start
# background threads and forking afterwards is unsafe (SIGILL), and forked
# workers can remove the shared session tempdir on exit. Isolating each cell in
# a spawned process avoids both: a child does exactly one backend (mirai OR
# fork, never both), and the parent samples its memory via /proc without forking.

stopifnot(requireNamespace("callr", quietly = TRUE))

bench_make_data <- function(n, p) {
  X <- as.data.frame(matrix(stats::rnorm(n * p), n, p))
  X$grp <- factor(sample(letters[1:6], n, replace = TRUE))
  X
}

## ---- memory metric ----------------------------------------------------------
# Preference order:
# 1. "cgroup-anon+shmem": sample the `anon` and `shmem` counters of our cgroup
#    v2 memory.stat. The kernel charges each page ONCE per cgroup, so COW pages
#    shared by fork workers count once -- exactly the unique-memory total we
#    want, and fairer than summed PSS (which only approximates sharing). mori's
#    regions live in /dev/shm (tmpfs), which the kernel files under `shmem`,
#    not `anon`; without that counter the shared input copy is invisible and
#    mirai is undercounted.
#    Every process a cell spawns (callr child, mirai dispatcher + daemons,
#    foreach fork workers) inherits the cgroup, so short-lived workers cannot
#    escape the measurement and no pgid matching is needed. `anon` excludes
#    page cache (readRDS pulling the data file in would otherwise inflate small
#    cells). One small file read per sample, so we can poll at 5ms instead of
#    the >=50ms-plus-VMA-walk sweeps of the PSS path, shrinking the missed-
#    spike window. Measured as a delta against the pre-spawn baseline because
#    the orchestrator shares the cgroup (slurm gives each job one cgroup).
#    Upgrade path (not needed yet): a dedicated per-cell child cgroup would
#    give kernel-exact memory.peak with no sampling at all, but requires
#    cgroupfs write delegation, which cluster nodes rarely grant.
# 2. "PSS": sum Pss over the child's process group, sampled. smaps_rollup
#    reads force VMA walks, so the effective sampling interval grows with
#    memory size and forked workers can spawn and die between sweeps.
#    Fallback for hosts without cgroup v2 memory accounting.
# 3. "RSS": last resort without smaps_rollup; overcounts shared pages.

# Read a /proc file, quietly tolerating the race where a pid vanishes between
# being listed and being read (returns character(0) then).
.bench_read <- function(path) {
  tryCatch(suppressWarnings(readLines(path, warn = FALSE)), error = function(e) character(0))
}

# cgroup v2 dir of this process, or NULL if memory accounting is unavailable.
.bench_cgroup_dir <- function() {
  cg <- .bench_read("/proc/self/cgroup")
  line <- cg[startsWith(cg, "0::")]
  if (!length(line)) {
    return(NULL)
  }
  dir <- file.path("/sys/fs/cgroup", sub("^0::/?", "", line[1]))
  if (file.exists(file.path(dir, "memory.stat"))) dir else NULL
}
.bench_cgroup_anon_kb <- function(dir) {
  st <- .bench_read(file.path(dir, "memory.stat"))
  a <- st[startsWith(st, "anon ")]
  h <- st[startsWith(st, "shmem ")]
  if (!length(a)) {
    return(NA_real_)
  }
  if (!length(h)) {
    h <- "shmem 0"
  }
  (as.numeric(sub("^anon ", "", a[1])) + as.numeric(sub("^shmem ", "", h[1]))) / 1024 # bytes -> kB
}

# The memory ceiling this process runs under, in MB, or NA when unlimited.
# Walk UP the hierarchy: slurm sets the limit on the job's cgroup while the
# process often sits in a child whose own memory.max reads "max". A peak that
# reaches this is a measurement of the cap, not of the workload, and must never
# be reported as the latter.
.bench_cgroup_limit_mb <- function(dir = .bench_cgroup_dir()) {
  if (is.null(dir)) {
    return(NA_real_)
  }
  root <- "/sys/fs/cgroup"
  lim <- Inf
  d <- dir
  repeat {
    for (nm in c("memory.max", "memory.high")) {
      v <- .bench_read(file.path(d, nm))
      if (length(v) && !identical(v[1], "max")) {
        n <- suppressWarnings(as.numeric(v[1]))
        if (!is.na(n) && n > 0) {
          lim <- min(lim, n / 1024 / 1024)
        }
      }
    }
    if (nchar(d) <= nchar(root)) {
      break
    }
    d <- dirname(d)
  }
  if (is.finite(lim)) lim else NA_real_
}

BENCH_CGROUP <- .bench_cgroup_dir()
BENCH_MEM_LIMIT_MB <- .bench_cgroup_limit_mb()
BENCH_USE_PSS <- file.exists("/proc/self/smaps_rollup")
BENCH_METRIC <- if (!is.null(BENCH_CGROUP)) {
  "cgroup-anon+shmem"
} else if (BENCH_USE_PSS) {
  "PSS"
} else {
  "RSS"
}
.bench_pid_kb <- function(pp) {
  if (BENCH_USE_PSS) {
    l <- .bench_read(sprintf("/proc/%d/smaps_rollup", pp))
    x <- l[grepl("^Pss:", l)]
    if (length(x)) {
      return(as.numeric(sub("[^0-9]*([0-9]+).*", "\\1", x[1])))
    }
    0
  } else {
    st <- .bench_read(sprintf("/proc/%d/statm", pp))
    if (!length(st)) {
      return(0)
    }
    as.numeric(strsplit(st, " ")[[1]][2]) * 4 # RSS pages -> kB (4k pages)
  }
}
# Process-group id of a pid (field after "pid (comm) state ppid" in /proc/stat).
# mirai daemons reparent to init (ppid=1) but KEEP the launcher's process group,
# so ppid-walking misses them while foreach fork workers are caught -- an unfair
# undercount. Grouping by pgid catches both. callr r_bg gives each cell child its
# own pgid, so this never sweeps in the parent or other cells.
.bench_pgid <- function(pid) {
  st <- .bench_read(sprintf("/proc/%d/stat", pid))
  if (!length(st)) {
    return(NA_integer_)
  }
  as.integer(strsplit(sub("^.*\\) \\S+ ", "", st[1]), " ", fixed = TRUE)[[1]][3])
}
# pids sharing `root`'s process group: the orchestrator, mirai dispatcher +
# daemons, and foreach fork workers.
.bench_group <- function(root) {
  g <- .bench_pgid(root)
  if (is.na(g)) {
    return(root)
  }
  pids <- suppressWarnings(as.integer(list.files("/proc")))
  pids <- pids[!is.na(pids)]
  pids[vapply(pids, function(p) isTRUE(.bench_pgid(p) == g), logical(1))]
}
.bench_tree_kb <- function(root) sum(vapply(.bench_group(root), .bench_pid_kb, numeric(1)))

# .libPaths() silently DROPS a directory that does not exist, so a ref whose
# install failed would quietly load whatever arf is installed in the session and
# report its number under the ref's label. Every process that loads a ref's arf
# asserts provenance: the master, the PSOCK workers and the mirai daemons.
.bench_assert_lib <- function(lib) {
  from <- normalizePath(find.package("arf"), mustWork = FALSE)
  want <- normalizePath(lib, mustWork = FALSE)
  if (!startsWith(from, want)) {
    stop("arf loaded from ", from, ", not from the ref library ", want, call. = FALSE)
  }
  invisible(TRUE)
}

# Function executed in the CHILD process (fresh R): load the package, run one
# operation `iters` times for one backend, return the elapsed seconds per iter.
# dt_threads pins data.table (per worker); ranger_threads sets ranger's default
# threads, which forde()'s terminalNodes prediction picks up (it does not pass
# num.threads itself).
#
# `op` selects the pipeline stage; `op_args` carries its knobs. The cached data
# object holds arf/X/psi/evidence; each op reads what it needs:
#   forde          arf, X                       parallelizes over trees
#   forge          psi, evidence                over evidence steps
#   expct          psi, evidence                over evidence steps
#   lik            psi, X, arf, batch           over query folds
#   adversarial_rf X, trees                     ranger train + prune (over trees)
# For mirai, the arf-loaded daemons are set up once here so per-call reloading
# doesn't dominate; see arf_load_on_daemons().
.bench_cell_fn <- function(
  lib,
  data_path,
  backend,
  n_workers,
  dt_threads,
  ranger_threads,
  iters,
  op = "forde",
  op_args = list()
) {
  # Prepend, never replace: the ref's library holds only arf, and the session's
  # own libpaths are the shared dependency layer, so a comparison cannot be
  # confounded by a different data.table version.
  .libPaths(c(lib, .libPaths()))
  suppressWarnings(suppressMessages(library(arf)))
  .bench_assert_lib(lib)
  data.table::setDTthreads(dt_threads)
  options(ranger.num.threads = ranger_threads)
  # speed-for-memory knob (see ?arf-options), gridable via env
  br <- Sys.getenv("ARF_BENCH_BLOCK_ROWS", "")
  if (nzchar(br)) {
    options(arf.block_rows = as.numeric(br))
  }
  d <- readRDS(data_path)
  arf <- d$arf
  X <- d$X
  psi <- d$psi
  evidence <- d$evidence
  if (backend == "sequential") {
    options(arf.backend = NULL)
    par <- FALSE
  } else if (backend == "foreach") {
    options(arf.backend = "foreach")
    par <- TRUE
    doParallel::registerDoParallel(cores = n_workers)
  } else if (backend == "psock") {
    # clean (non-fork) foreach workers: separates fork's heap duplication from
    # the per-worker input copies that mori removes
    options(arf.backend = "foreach")
    par <- TRUE
    cl <- parallel::makeCluster(n_workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    # Workers are fresh Rscript processes whose .libPaths() comes from env vars
    # only, so send the master's whole set: otherwise the shared dependency
    # layer can differ between master and worker even when `lib` is fine.
    parallel::clusterCall(
      cl,
      function(paths, l, t, assert) {
        .libPaths(paths)
        suppressMessages(library(arf))
        assert(l)
        data.table::setDTthreads(t)
      },
      .libPaths(),
      lib,
      dt_threads,
      .bench_assert_lib
    )
    doParallel::registerDoParallel(cl)
  } else if (backend == "mirai") {
    options(arf.backend = "mirai")
    par <- TRUE
    mirai::daemons(n_workers)
    on.exit(mirai::daemons(0), add = TRUE) # always stop, even if the op errors
    mirai::everywhere(data.table::setDTthreads(dt_threads))
    # load the ref's arf on daemons once (forge/expct/lik/cforde workers call
    # arf internals); prune passes its worker as an object so needs no load.
    if (op %in% c("forge", "expct", "lik")) {
      mirai::everywhere(
        {
          .libPaths(paths)
          suppressMessages(library(arf))
          assert(lib)
        },
        paths = .libPaths(),
        lib = lib,
        assert = .bench_assert_lib
      )
    }
  } else {
    stop("unknown backend: ", backend)
  }
  run <- switch(
    op,
    forde = function() forde(arf, X, parallel = par),
    forge = function() {
      forge(
        psi,
        n_synth = op_args$n_synth %||% 1L,
        evidence = evidence,
        parallel = par,
        evidence_row_mode = op_args$rowmode %||% "separate",
        stepsize = op_args$stepsize %||% 0L,
        verbose = FALSE
      )
    },
    expct = function() {
      expct(
        psi,
        evidence = evidence,
        parallel = par,
        evidence_row_mode = op_args$rowmode %||% "separate",
        stepsize = op_args$stepsize %||% 0L,
        verbose = FALSE
      )
    },
    lik = function() lik(psi, X, arf = arf, batch = op_args$batch, parallel = par),
    adversarial_rf = function() adversarial_rf(X, num_trees = op_args$trees, parallel = par, verbose = FALSE),
    stop("unknown op: ", op)
  )
  kind <- unname(BENCH_DIGEST_KIND[[op]])
  # Digest the first iteration only: one hash per cell is enough, and a fixed
  # seed makes the stochastic ops comparable across refs.
  set.seed(1)
  first <- NULL
  secs <- vapply(
    seq_len(iters),
    function(i) {
      t <- system.time(res <- run())[["elapsed"]]
      if (i == 1L) {
        first <<- res
      }
      t
    },
    numeric(1)
  )
  list(
    seconds = secs,
    digest = tryCatch(bench_digest(first, kind), error = function(e) NA_character_),
    digest_kind = kind,
    arf_version = as.character(utils::packageVersion("arf"))
  )
}
`%||%` <- function(a, b) if (is.null(a)) b else a

# Parent-side: launch a cell in a fresh subprocess and track its peak memory
# while it runs (cgroup-anon delta when available, else /proc PSS/RSS sweeps;
# see the metric note above). Returns list(seconds = median, peak_mb).
bench_measure_cell <- function(
  backend,
  data_path,
  n_workers,
  dt_threads,
  lib,
  ranger_threads = 1L,
  iters = 1L,
  interval = 0.05,
  op = "forde",
  op_args = list(),
  mem_reps = 1L
) {
  args <- .bench_cell_args(
    lib,
    data_path,
    backend,
    n_workers,
    dt_threads,
    ranger_threads,
    iters,
    op,
    op_args
  )
  runs <- lapply(seq_len(max(1L, mem_reps)), function(i) {
    .bench_run_sampled(.bench_cell_payload(), args, interval)
  })
  r <- runs[[1]]
  if (isTRUE(r$ok)) {
    # Peak memory depends on when R's GC decides to grow the heap, so repeated
    # runs of identical code differ by tens of percent on the marginal. The
    # minimum is the least GC-inflated observation of the real requirement.
    peaks <- vapply(runs, function(x) if (isTRUE(x$ok)) x$peak_mb else NA_real_, numeric(1))
    r$peak_mb <- min(peaks, na.rm = TRUE)
    secs <- unlist(lapply(runs, function(x) if (isTRUE(x$ok)) x$value$seconds else NULL))
    r$value$seconds <- secs
  }
  if (!isTRUE(r$ok)) {
    # Whatever the sampler accumulated before the child died is not a
    # measurement of anything: reported as a number it renders as a large
    # memory "improvement" with verdict "real".
    return(list(
      seconds = NA_real_,
      peak_mb = NA_real_,
      digest = NA_character_,
      digest_kind = NA_character_,
      arf_version = NA_character_
    ))
  }
  list(
    seconds = r$value$seconds,
    peak_mb = r$peak_mb,
    digest = r$value$digest,
    digest_kind = r$value$digest_kind,
    arf_version = r$value$arf_version
  )
}

# Every child in a cell shares one cgroup, and the counter is read as a delta
# against a baseline taken before the child spawns. A previous child's pages
# are not reclaimed the instant it exits, so a baseline read too early is
# inflated and the next child's delta comes out too small -- which produced
# floor measurements LARGER than the peak they were subtracted from, i.e. a
# negative marginal. A fixed sleep is not enough; wait for the counter itself
# to stop moving.
.bench_settle <- function(max_wait = 3, tol_kb = 1024) {
  if (is.null(BENCH_CGROUP)) {
    return(invisible(NULL))
  }
  invisible(gc(FALSE))
  prev <- .bench_cgroup_anon_kb(BENCH_CGROUP)
  deadline <- Sys.time() + max_wait
  repeat {
    Sys.sleep(0.05)
    cur <- .bench_cgroup_anon_kb(BENCH_CGROUP)
    if (is.na(cur) || abs(cur - prev) < tol_kb || Sys.time() > deadline) {
      break
    }
    prev <- cur
  }
  invisible(NULL)
}

# Spawn `fn` in a fresh child (callr spawns, it does not fork) and track the
# cell's peak memory while it runs. Shared by the cell measurement and the
# floor measurement so both use one sampling path.
.bench_run_sampled <- function(fn, args, interval = 0.05) {
  use_cgroup <- !is.null(BENCH_CGROUP)
  baseline_kb <- NA_real_
  if (use_cgroup) {
    # Stabilize the orchestrator's share before taking the baseline: a GC
    # during the cell would deflate the delta's floor (harmless for a max),
    # but unreclaimed garbage at baseline time would inflate every sample.
    .bench_settle()
    baseline_kb <- .bench_cgroup_anon_kb(BENCH_CGROUP)
    interval <- 0.005 # one counter read per sample; poll fast
  }
  # callr ships the function and its arguments, not the parent's globals, so the
  proc <- callr::r_bg(fn, args = args)
  pid <- proc$get_pid()
  peak_kb <- 0
  repeat {
    alive <- proc$is_alive()
    cur_kb <- if (use_cgroup) {
      .bench_cgroup_anon_kb(BENCH_CGROUP) - baseline_kb
    } else {
      .bench_tree_kb(pid)
    }
    if (!is.na(cur_kb)) {
      peak_kb <- max(peak_kb, cur_kb)
    }
    if (!alive) {
      break
    }
    Sys.sleep(interval)
  }
  res <- tryCatch(list(ok = TRUE, value = proc$get_result()), error = function(e) {
    list(ok = FALSE, value = conditionMessage(e))
  })
  list(ok = res$ok, value = res$value, peak_mb = peak_kb / 1024)
}

# callr ships the function and its arguments, not the parent's globals, so the
# cell function and the digest helpers travel as an explicit payload.
.bench_cell_payload <- function() {
  function(lib, data_path, backend, n_workers, dt_threads, ranger_threads, iters, op, op_args, helpers) {
    for (nm in names(helpers)) {
      assign(nm, helpers[[nm]], envir = globalenv())
    }
    .bench_cell_fn(lib, data_path, backend, n_workers, dt_threads, ranger_threads, iters, op, op_args)
  }
}

.bench_cell_args <- function(lib, data_path, backend, n_workers, dt_threads, ranger_threads, iters, op, op_args) {
  list(
    lib,
    data_path,
    backend,
    n_workers,
    dt_threads,
    ranger_threads,
    iters,
    op,
    op_args,
    helpers = list(
      .bench_cell_fn = .bench_cell_fn,
      .bench_assert_lib = .bench_assert_lib,
      bench_digest = bench_digest,
      .bench_digest_exact = .bench_digest_exact,
      .bench_digest_summary = .bench_digest_summary,
      .bench_round_num = .bench_round_num,
      BENCH_DIGEST_KIND = BENCH_DIGEST_KIND,
      `%||%` = `%||%`
    )
  )
}

# The floor of a cell: an R interpreter, arf, and the fixture, with the op never
# run. Peak memory is dominated by it (hundreds of MB), so a percentage taken
# on the raw peak is diluted several-fold and the 3% "real" threshold cannot
# see a genuine regression. Subtracting the floor gives the op's own marginal.
bench_measure_floor <- function(lib, data_path, interval = 0.05, reps = 1L) {
  if (reps > 1L) {
    peaks <- vapply(
      seq_len(reps),
      function(i) {
        bench_measure_floor(lib, data_path, interval, reps = 1L)
      },
      numeric(1)
    )
    return(min(peaks, na.rm = TRUE))
  }
  r <- .bench_run_sampled(
    function(lib, data_path, assert) {
      .libPaths(c(lib, .libPaths()))
      suppressWarnings(suppressMessages(library(arf)))
      assert(lib)
      invisible(readRDS(data_path))
      NULL
    },
    args = list(lib, data_path, .bench_assert_lib),
    interval = interval
  )
  if (!isTRUE(r$ok)) NA_real_ else r$peak_mb
}

## ---- correctness digest -----------------------------------------------------
# A performance number from a ref that computes something different is worse
# than no number, so every cell fingerprints its own result. "exact" hashes the
# rounded values (forde parameters, lik values). "summary" hashes rounded column
# means and standard deviations, for ops whose output is a random sample
# (forge, expct): an exact hash there also moves when a refactor merely
# consumes RNG draws in a different order.
#
# adversarial_rf is deliberately "none": its result is a fitted ranger forest,
# where an exact hash would flag benign RNG-consumption differences and the
# summary path has no as.data.frame() method to stand on. The fit's correctness
# belongs to the test suite, not to a benchmark fingerprint.
BENCH_DIGEST_KIND <- c(
  forde = "exact",
  lik = "exact",
  adversarial_rf = "none",
  forge = "summary",
  expct = "summary"
)

.bench_round_num <- function(x, digits) {
  if (is.numeric(x)) round(x, digits) else x
}

.bench_digest_exact <- function(x) {
  if (is.null(x)) {
    return("NULL")
  }
  if (is.list(x) && !is.data.frame(x)) {
    ordered <- if (is.null(names(x))) x else x[order(names(x))]
    return(lapply(ordered, .bench_digest_exact))
  }
  if (is.data.frame(x)) {
    return(lapply(as.list(x)[order(names(x))], .bench_round_num, digits = 8))
  }
  .bench_round_num(x, 8)
}

# Order-free on purpose: row order is not part of the distribution, and sample
# order is not stable across refactors.
.bench_digest_summary <- function(x) {
  if (is.null(x)) {
    return("NULL")
  }
  x <- as.data.frame(x)
  if (!nrow(x) || !ncol(x)) {
    return("empty")
  }
  cols <- as.list(x)[order(names(x))]
  lapply(cols, function(col) {
    if (is.numeric(col)) {
      round(c(mean(col), stats::sd(col)), 4)
    } else {
      tab <- table(as.character(col))
      round(as.numeric(tab[order(names(tab))]) / length(col), 4)
    }
  })
}

bench_digest <- function(x, kind = c("exact", "summary", "none")) {
  kind <- match.arg(kind)
  if (kind == "none") {
    return(NA_character_)
  }
  payload <- if (kind == "exact") .bench_digest_exact(x) else .bench_digest_summary(x)
  digest::digest(payload, algo = "xxhash64")
}
