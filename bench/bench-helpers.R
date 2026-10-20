# Shared helpers for the bench/ scripts (sourced, never shipped).
# Every measurement runs in a fresh R process via callr (spawned, not forked):
# mirai starts background threads, and forking after that is unsafe.

stopifnot(requireNamespace("callr", quietly = TRUE))

# Absolute, so a measurement child can source this file whatever its cwd.
BENCH_HELPERS <- normalizePath("bench/bench-helpers.R", mustWork = TRUE)

bench_make_data <- function(n, p) {
  X <- as.data.frame(matrix(stats::rnorm(n * p), n, p))
  X$grp <- factor(sample(letters[1:6], n, replace = TRUE))
  X
}

## ---- memory metric ----------------------------------------------------------
# cgroup v2 `anon` + `shmem` of the job's cgroup: each page is charged once, so
# pages shared by fork workers count once and mori's /dev/shm segments count.
# Every process a cell spawns inherits the cgroup.

# A file can vanish between listing and reading; treat that as an empty read.
.bench_read <- function(path) {
  tryCatch(suppressWarnings(readLines(path, warn = FALSE)), error = function(e) character(0))
}

# cgroup v2 dir of this process, or NULL without memory accounting.
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
  (as.numeric(sub("^anon ", "", a[1])) + as.numeric(sub("^shmem ", "", h[1]))) / 1024
}

# Walk up: slurm sets the limit on the job cgroup while the process sits in a
# child whose own memory.max reads "max". NA when unlimited.
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
BENCH_METRIC <- "cgroup-anon+shmem"

# No fallback metric: PSS and RSS count shared pages differently, so a number
# from either would not be comparable with the rest of the results.
bench_assert_cgroup <- function() {
  if (is.null(BENCH_CGROUP)) {
    stop(
      "no cgroup v2 memory accounting for this process; run under slurm or a systemd slice",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

## ---- child-side setup ------------------------------------------------------

# .libPaths() silently drops a missing directory, so a ref whose install failed
# would load whatever arf the session has. Every process that loads arf checks.
.bench_assert_lib <- function(lib) {
  from <- normalizePath(find.package("arf"), mustWork = FALSE)
  want <- normalizePath(lib, mustWork = FALSE)
  if (!startsWith(from, want)) {
    stop("arf loaded from ", from, ", not from the ref library ", want, call. = FALSE)
  }
  invisible(TRUE)
}

# everywhere() evaluates on the daemon, so the count travels as a value.
.bench_mirai_dt_threads <- function(dt_threads) {
  mirai::everywhere(data.table::setDTthreads(.n), .n = dt_threads)
}

# Daemons take their library paths from the environment, so the ref library
# has to be set on each; worker bodies resolve arf internals against it.
.bench_mirai_load_arf <- function(lib) {
  mirai::everywhere(
    {
      .libPaths(.paths)
      suppressMessages(library(arf))
      .assert(.lib)
    },
    .paths = .libPaths(),
    .lib = lib,
    .assert = .bench_assert_lib
  )
}

# PSOCK workers get the master's whole libpath set, so the shared dependency
# layer matches the master even when the ref library alone would.
.bench_psock_setup <- function(cl, lib, dt_threads) {
  setup <- function(paths, l, t, assert) {
    .libPaths(paths)
    suppressMessages(library(arf))
    assert(l)
    data.table::setDTthreads(t)
  }
  # A closure ships its frame, and this frame holds the fixture.
  environment(setup) <- globalenv()
  parallel::clusterCall(
    cl,
    setup,
    .libPaths(),
    lib,
    dt_threads,
    .bench_assert_lib
  )
}

# everywhere() is fire-and-forget, so a failed daemon setup is checked here.
.bench_mirai_verify <- function(lib, dt_threads) {
  n <- max(1L, as.integer(mirai::status()$connections))
  probe <- function(i) c(dirname(find.package("arf")), data.table::getDTthreads())
  # A closure ships its frame, and this frame holds the fixture.
  environment(probe) <- globalenv()
  got <- mirai::mirai_map(seq_len(n), probe)[]
  want <- normalizePath(lib)
  for (g in got) {
    if (!startsWith(normalizePath(g[1]), want)) {
      stop("daemon loaded arf from ", g[1], ", not from ", want, call. = FALSE)
    }
    if (as.integer(g[2]) != as.integer(dt_threads)) {
      stop(
        "daemon has ",
        g[2],
        " data.table thread(s), not ",
        dt_threads,
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}

# Runs in the measurement child: load the ref's arf, set up one backend, time
# `iters` calls of `op`. ranger_threads sets ranger's default, which forde's
# terminalNodes prediction picks up.
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
  # Prepend: the ref library holds only arf, the session's paths hold the deps.
  .libPaths(c(lib, .libPaths()))
  suppressWarnings(suppressMessages(library(arf)))
  .bench_assert_lib(lib)
  data.table::setDTthreads(dt_threads)
  options(ranger.num.threads = ranger_threads)
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
    # Non-fork foreach workers: no copy-on-write sharing with the parent.
    options(arf.backend = "foreach")
    par <- TRUE
    cl <- parallel::makeCluster(n_workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    .bench_psock_setup(cl, lib, dt_threads)
    doParallel::registerDoParallel(cl)
  } else if (backend == "mirai") {
    options(arf.backend = "mirai")
    par <- TRUE
    mirai::daemons(n_workers)
    on.exit(mirai::daemons(0), add = TRUE)
    .bench_mirai_dt_threads(dt_threads)
    .bench_mirai_load_arf(lib)
    .bench_mirai_verify(lib, dt_threads)
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
  # One digest per cell, from the first call under a fixed seed.
  set.seed(1)
  digest <- NA_character_
  out_rows <- NA_integer_
  secs <- vapply(
    seq_len(iters),
    function(i) {
      t <- system.time(res <- run())[["elapsed"]]
      if (i == 1L) {
        digest <<- tryCatch(bench_digest(res, kind), error = function(e) NA_character_)
        # A worker that dies drops its chunk and arf returns the short result.
        out_rows <<- tryCatch(
          if (is.null(res)) NA_integer_ else as.integer(NROW(res)),
          error = function(e) NA_integer_
        )
      }
      # Dropped before the next call, or it sits in every later peak.
      rm(res)
      invisible(gc(FALSE))
      t
    },
    numeric(1)
  )
  list(
    seconds = secs,
    digest = digest,
    digest_kind = kind,
    out_rows = out_rows,
    arf_version = as.character(utils::packageVersion("arf"))
  )
}
`%||%` <- function(a, b) if (is.null(a)) b else a

## ---- parent-side measurement -------------------------------------------------

# Launch one cell measurement in a fresh child and track its peak memory.
bench_measure_cell <- function(
  backend,
  data_path,
  n_workers,
  dt_threads,
  lib,
  ranger_threads = 1L,
  iters = 1L,
  op = "forde",
  op_args = list()
) {
  r <- .bench_run_sampled(
    function(helpers, args) {
      source(helpers)
      do.call(.bench_cell_fn, args)
    },
    args = list(
      BENCH_HELPERS,
      list(
        lib = lib,
        data_path = data_path,
        backend = backend,
        n_workers = n_workers,
        dt_threads = dt_threads,
        ranger_threads = ranger_threads,
        iters = iters,
        op = op,
        op_args = op_args
      )
    )
  )
  if (!isTRUE(r$ok)) {
    # Without the reason a dead child looks like a slow one in the job log.
    message("  child failed: ", paste(r$value, collapse = " "))
    return(list(
      seconds = NA_real_,
      peak_mb = NA_real_,
      digest = NA_character_,
      digest_kind = NA_character_,
      out_rows = NA_integer_,
      arf_version = NA_character_
    ))
  }
  list(
    seconds = r$value$seconds,
    peak_mb = r$peak_mb,
    digest = r$value$digest,
    digest_kind = r$value$digest_kind,
    out_rows = r$value$out_rows,
    arf_version = r$value$arf_version
  )
}

# A previous child's pages are reclaimed after it exits, not at exit; a
# baseline read before the counter stops moving inflates it.
.bench_settle <- function(max_wait = 3, tol_kb = 1024) {
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

# Spawn `fn` in a fresh child and sample the cgroup delta while it runs.
.bench_run_sampled <- function(fn, args) {
  bench_assert_cgroup()
  .bench_settle()
  baseline_kb <- .bench_cgroup_anon_kb(BENCH_CGROUP)
  proc <- callr::r_bg(fn, args = args)
  peak_kb <- 0
  repeat {
    alive <- proc$is_alive()
    cur_kb <- .bench_cgroup_anon_kb(BENCH_CGROUP) - baseline_kb
    if (!is.na(cur_kb)) {
      peak_kb <- max(peak_kb, cur_kb)
    }
    if (!alive) {
      break
    }
    Sys.sleep(0.005)
  }
  res <- tryCatch(list(ok = TRUE, value = proc$get_result()), error = function(e) {
    list(ok = FALSE, value = conditionMessage(e))
  })
  list(ok = res$ok, value = res$value, peak_mb = peak_kb / 1024)
}

# The floor: R, arf and the fixture with the op never run. The op's marginal is
# the peak minus this, so a 3% threshold is not diluted by the fixed cost.
bench_measure_floor <- function(lib, data_path) {
  r <- .bench_run_sampled(
    function(lib, data_path, assert) {
      .libPaths(c(lib, .libPaths()))
      suppressWarnings(suppressMessages(library(arf)))
      assert(lib)
      invisible(readRDS(data_path))
      NULL
    },
    args = list(lib, data_path, .bench_assert_lib)
  )
  if (!isTRUE(r$ok)) NA_real_ else r$peak_mb
}

## ---- correctness digest -----------------------------------------------------
# "exact" hashes rounded values. "summary" hashes rounded column means and sds,
# for ops whose output is a random sample: an exact hash there also moves when
# a refactor merely consumes RNG draws in a different order. adversarial_rf
# returns a ranger forest, which neither path can fingerprint usefully.
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

# Order-free: row order is not part of the distribution.
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
