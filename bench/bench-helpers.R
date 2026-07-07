# Shared helpers for the bench/ scripts (sourced, not part of the package build).
#
# Each backend measurement runs in its OWN fresh subprocess (callr spawns a new
# R process rather than forking). This is deliberate: mirai/nanonext start
# background threads and forking afterwards is unsafe (SIGILL), and forked
# workers can remove the shared session tempdir on exit. Isolating each cell in
# a spawned process avoids both: a child does exactly one backend (mirai OR
# fork, never both), and the parent samples its memory via /proc without forking.

suppressWarnings(suppressMessages(pkgload::load_all(quiet = TRUE)))
stopifnot(requireNamespace("callr", quietly = TRUE))

# Abort unless every backend package is available -- a backend comparison is
# meaningless otherwise, and this fails BEFORE any expensive model fit.
bench_require_backends <- function() {
  missing <- c(
    if (!requireNamespace("doParallel", quietly = TRUE)) "doParallel",
    if (!requireNamespace("mirai", quietly = TRUE)) "mirai",
    if (!requireNamespace("mori", quietly = TRUE)) "mori"
  )
  if (length(missing)) {
    stop("benchmark needs all backends installed; missing: ",
         paste(missing, collapse = ", "),
         ". Install them (comparing backends is the whole point).", call. = FALSE)
  }
  invisible(TRUE)
}

# Comma-separated integer env override, e.g. ARF_BENCH_WORKERS=1,2,4,8.
bench_ints <- function(env, default) {
  v <- Sys.getenv(env, "")
  if (nzchar(v)) as.integer(strsplit(v, ",")[[1]]) else as.integer(default)
}

# Short commit hash of the benchmarked arf tree, "+dirty" if it has uncommitted
# changes -- recorded per CSV row so results from different code states are
# never silently compared. "unknown" outside a git checkout.
bench_git_commit <- function() {
  hash <- tryCatch(
    system2("git", c("rev-parse", "--short", "HEAD"),
            stdout = TRUE, stderr = FALSE)[1],
    error = function(e) NA_character_, warning = function(w) NA_character_)
  if (is.na(hash) || !nzchar(hash)) return("unknown")
  dirty <- tryCatch(
    length(system2("git", c("status", "--porcelain", "--untracked-files=no"),
                   stdout = TRUE, stderr = FALSE)) > 0L,
    error = function(e) FALSE, warning = function(w) FALSE)
  paste0(hash, if (dirty) "+dirty" else "")
}
BENCH_COMMIT <- bench_git_commit()

bench_make_data <- function(n, p) {
  X <- as.data.frame(matrix(stats::rnorm(n * p), n, p))
  X$grp <- factor(sample(letters[1:6], n, replace = TRUE))
  X
}

## ---- memory metric ----------------------------------------------------------
# Preference order:
# 1. "cgroup-anon": sample the `anon` counter of our cgroup v2 memory.stat.
#    The kernel charges each page ONCE per cgroup, so COW pages shared by fork
#    workers and mori-shared pages count once -- exactly the unique-memory total
#    we want, and fairer than summed PSS (which only approximates sharing).
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
  tryCatch(suppressWarnings(readLines(path, warn = FALSE)),
           error = function(e) character(0))
}

# cgroup v2 dir of this process, or NULL if memory accounting is unavailable.
.bench_cgroup_dir <- function() {
  cg <- .bench_read("/proc/self/cgroup")
  line <- cg[startsWith(cg, "0::")]
  if (!length(line)) return(NULL)
  dir <- file.path("/sys/fs/cgroup", sub("^0::/?", "", line[1]))
  if (file.exists(file.path(dir, "memory.stat"))) dir else NULL
}
.bench_cgroup_anon_kb <- function(dir) {
  st <- .bench_read(file.path(dir, "memory.stat"))
  x <- st[startsWith(st, "anon ")]
  if (!length(x)) return(NA_real_)
  as.numeric(sub("^anon ", "", x[1])) / 1024  # bytes -> kB
}

BENCH_CGROUP  <- .bench_cgroup_dir()
BENCH_USE_PSS <- file.exists("/proc/self/smaps_rollup")
BENCH_METRIC  <- if (!is.null(BENCH_CGROUP)) "cgroup-anon" else if (BENCH_USE_PSS) "PSS" else "RSS"
.bench_pid_kb <- function(pp) {
  if (BENCH_USE_PSS) {
    l <- .bench_read(sprintf("/proc/%d/smaps_rollup", pp))
    x <- l[grepl("^Pss:", l)]
    if (length(x)) return(as.numeric(sub("[^0-9]*([0-9]+).*", "\\1", x[1])))
    0
  } else {
    st <- .bench_read(sprintf("/proc/%d/statm", pp))
    if (!length(st)) return(0)
    as.numeric(strsplit(st, " ")[[1]][2]) * 4  # RSS pages -> kB (4k pages)
  }
}
# Process-group id of a pid (field after "pid (comm) state ppid" in /proc/stat).
# mirai daemons reparent to init (ppid=1) but KEEP the launcher's process group,
# so ppid-walking misses them while foreach fork workers are caught -- an unfair
# undercount. Grouping by pgid catches both. callr r_bg gives each cell child its
# own pgid, so this never sweeps in the parent or other cells.
.bench_pgid <- function(pid) {
  st <- .bench_read(sprintf("/proc/%d/stat", pid))
  if (!length(st)) return(NA_integer_)
  as.integer(strsplit(sub("^.*\\) \\S+ ", "", st[1]), " ", fixed = TRUE)[[1]][3])
}
# pids sharing `root`'s process group: the orchestrator, mirai dispatcher +
# daemons, and foreach fork workers.
.bench_group <- function(root) {
  g <- .bench_pgid(root)
  if (is.na(g)) return(root)
  pids <- suppressWarnings(as.integer(list.files("/proc")))
  pids <- pids[!is.na(pids)]
  pids[vapply(pids, function(p) isTRUE(.bench_pgid(p) == g), logical(1))]
}
.bench_tree_kb <- function(root) sum(vapply(.bench_group(root), .bench_pid_kb, numeric(1)))

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
.bench_cell_fn <- function(pkgdir, data_path, backend, n_workers, dt_threads,
                           ranger_threads, iters, op = "forde",
                           op_args = list()) {
  suppressWarnings(suppressMessages(pkgload::load_all(pkgdir, quiet = TRUE)))
  data.table::setDTthreads(dt_threads)
  options(ranger.num.threads = ranger_threads)
  d <- readRDS(data_path)
  arf <- d$arf; X <- d$X; psi <- d$psi; evidence <- d$evidence
  if (backend == "sequential") {
    options(arf.backend = NULL); par <- FALSE
  } else if (backend == "foreach") {
    options(arf.backend = "foreach"); par <- TRUE
    doParallel::registerDoParallel(cores = n_workers)
  } else if (backend == "mirai") {
    options(arf.backend = "mirai"); par <- TRUE
    mirai::daemons(n_workers)
    on.exit(mirai::daemons(0), add = TRUE)  # always stop, even if the op errors
    mirai::everywhere(data.table::setDTthreads(dt_threads))
    # load the dev build on daemons once (forge/expct/lik/cforde workers call arf
    # internals); prune passes its worker as an object so needs no load.
    if (op %in% c("forge", "expct", "lik") &&
        isTRUE(tryCatch(pkgload::is_dev_package("arf"), error = function(e) FALSE))) {
      mirai::everywhere(suppressMessages(pkgload::load_all(pkgdir, quiet = TRUE)),
                        pkgdir = pkgdir)
    }
  } else {
    stop("unknown backend: ", backend)
  }
  run <- switch(op,
    forde = function() forde(arf, X, parallel = par),
    forge = function() forge(psi, n_synth = op_args$n_synth %||% 1L,
                             evidence = evidence, parallel = par,
                             evidence_row_mode = op_args$rowmode %||% "separate",
                             stepsize = op_args$stepsize %||% 0L, verbose = FALSE),
    expct = function() expct(psi, evidence = evidence, parallel = par,
                             evidence_row_mode = op_args$rowmode %||% "separate",
                             stepsize = op_args$stepsize %||% 0L, verbose = FALSE),
    lik   = function() lik(psi, X, arf = arf, batch = op_args$batch, parallel = par),
    adversarial_rf = function() adversarial_rf(X, num_trees = op_args$trees,
                                               parallel = par, verbose = FALSE),
    stop("unknown op: ", op))
  secs <- vapply(seq_len(iters),
                 function(i) system.time(run())[["elapsed"]], numeric(1))
  secs
}
`%||%` <- function(a, b) if (is.null(a)) b else a

# Parent-side: launch a cell in a fresh subprocess and track its peak memory
# while it runs (cgroup-anon delta when available, else /proc PSS/RSS sweeps;
# see the metric note above). Returns list(seconds = median, peak_mb).
bench_measure_cell <- function(backend, data_path, n_workers, dt_threads,
                               pkgdir, ranger_threads = 1L, iters = 1L,
                               interval = 0.05, op = "forde", op_args = list()) {
  use_cgroup <- !is.null(BENCH_CGROUP)
  if (use_cgroup) {
    # Stabilize the orchestrator's share before taking the baseline: a GC
    # during the cell would deflate the delta's floor (harmless for a max),
    # but unreclaimed garbage at baseline time would inflate every sample.
    invisible(gc(FALSE))
    baseline_kb <- .bench_cgroup_anon_kb(BENCH_CGROUP)
    interval <- 0.005  # one counter read per sample; poll fast
  }
  proc <- callr::r_bg(.bench_cell_fn,
                      args = list(pkgdir, data_path, backend, n_workers,
                                  dt_threads, ranger_threads, iters, op, op_args))
  pid <- proc$get_pid()
  peak_kb <- 0
  repeat {
    alive <- proc$is_alive()
    cur_kb <- if (use_cgroup) {
      .bench_cgroup_anon_kb(BENCH_CGROUP) - baseline_kb
    } else {
      .bench_tree_kb(pid)
    }
    if (!is.na(cur_kb)) peak_kb <- max(peak_kb, cur_kb)
    if (!alive) break
    Sys.sleep(interval)
  }
  secs <- tryCatch(proc$get_result(), error = function(e) NA_real_)
  list(seconds = stats::median(secs), peak_mb = peak_kb / 1024)
}
