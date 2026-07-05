# Shared helpers for the bench/ scripts (sourced, not part of the package build).
# Provides: backend requirement check, env grid parsing, data generation,
# peak-memory sampling across the R process tree, and a single-backend runner.

suppressWarnings(suppressMessages(pkgload::load_all(quiet = TRUE)))

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

bench_make_data <- function(n, p) {
  X <- as.data.frame(matrix(stats::rnorm(n * p), n, p))
  X$grp <- factor(sample(letters[1:6], n, replace = TRUE))
  X
}

## ---- peak memory across this user's R process tree -------------------------
# PSS (Proportional Set Size) on Linux -- fair for shared vs copied memory --
# else RSS. See bench/README.md for the rationale.
BENCH_USE_PSS <- file.exists("/proc/self/smaps_rollup")
BENCH_METRIC  <- if (BENCH_USE_PSS) "PSS" else "RSS"

.bench_r_pids <- function() {
  pids <- suppressWarnings(as.integer(list.files("/proc")))
  pids <- pids[!is.na(pids)]
  keep <- vapply(pids, function(pp) {
    comm <- tryCatch(readLines(sprintf("/proc/%d/comm", pp), warn = FALSE),
                     error = function(e) "")
    length(comm) && grepl("^R", comm)
  }, logical(1))
  pids[keep]
}
.bench_pid_kb <- function(pp) {
  if (BENCH_USE_PSS) {
    l <- tryCatch(readLines(sprintf("/proc/%d/smaps_rollup", pp), warn = FALSE),
                  error = function(e) character(0))
    x <- l[grepl("^Pss:", l)]
    if (length(x)) return(as.numeric(sub("[^0-9]*([0-9]+).*", "\\1", x[1])))
    0
  } else {
    st <- tryCatch(readLines(sprintf("/proc/%d/statm", pp), warn = FALSE),
                   error = function(e) character(0))
    if (!length(st)) return(0)
    as.numeric(strsplit(st, " ")[[1]][2]) * 4  # RSS pages -> kB (4k pages)
  }
}
bench_total_mem_kb <- function() {
  sum(vapply(.bench_r_pids(), .bench_pid_kb, numeric(1)))
}

# Run `expr` while a forked sampler tracks peak memory of the R process tree;
# returns list(seconds, peak_mb).
bench_with_peak <- function(expr, interval = 0.02) {
  stopifnot(requireNamespace("parallel", quietly = TRUE))
  stopf <- tempfile()
  sampler <- parallel::mcparallel({
    peak <- 0
    repeat {
      peak <- max(peak, bench_total_mem_kb())
      if (file.exists(stopf)) break
      Sys.sleep(interval)
    }
    peak
  })
  t <- system.time(force(expr))[["elapsed"]]
  file.create(stopf)
  peak_kb <- tryCatch(parallel::mccollect(sampler)[[1]], error = function(e) NA_real_)
  unlink(stopf)
  list(seconds = t, peak_mb = peak_kb / 1024)
}

# Run forde() once for `backend` at `n_workers`, measuring time + peak memory.
# Registers/tears down the backend around the call; pins data.table threads on
# mirai daemons to `dt_threads` (main + forks inherit the caller's setting).
bench_run_backend <- function(backend, arf, X, n_workers, dt_threads = 1L,
                              interval = 0.02) {
  options(arf.backend = NULL)
  on.exit(options(arf.backend = NULL), add = TRUE)
  if (backend == "sequential") {
    return(bench_with_peak(forde(arf, X, parallel = FALSE), interval))
  }
  if (backend == "foreach") {
    options(arf.backend = "foreach")
    doParallel::registerDoParallel(cores = n_workers)
    on.exit(try(doParallel::stopImplicitCluster(), silent = TRUE), add = TRUE)
    return(bench_with_peak(forde(arf, X, parallel = TRUE), interval))
  }
  if (backend == "mirai") {
    options(arf.backend = "mirai")
    mirai::daemons(n_workers)
    mirai::everywhere(data.table::setDTthreads(dt_threads))
    on.exit(mirai::daemons(0), add = TRUE)
    return(bench_with_peak(forde(arf, X, parallel = TRUE), interval))
  }
  stop("unknown backend: ", backend)
}
