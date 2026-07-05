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

bench_make_data <- function(n, p) {
  X <- as.data.frame(matrix(stats::rnorm(n * p), n, p))
  X$grp <- factor(sample(letters[1:6], n, replace = TRUE))
  X
}

## ---- memory: peak PSS (Linux) / RSS across a process subtree ---------------
BENCH_USE_PSS <- file.exists("/proc/self/smaps_rollup")
BENCH_METRIC  <- if (BENCH_USE_PSS) "PSS" else "RSS"

# Read a /proc file, quietly tolerating the race where a pid vanishes between
# being listed and being read (returns character(0) then).
.bench_read <- function(path) {
  tryCatch(suppressWarnings(readLines(path, warn = FALSE)),
           error = function(e) character(0))
}
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
# pids of `root` plus all descendants (walks /proc ppid links); catches foreach
# fork workers and mirai daemons spawned by the child.
.bench_descendants <- function(root) {
  pids <- suppressWarnings(as.integer(list.files("/proc")))
  pids <- pids[!is.na(pids)]
  ppid <- setNames(rep(NA_integer_, length(pids)), pids)
  for (pp in pids) {
    st <- .bench_read(sprintf("/proc/%d/stat", pp))
    if (!length(st)) next
    after <- sub("^\\d+ \\(.*\\) \\S+ ", "", st)  # strip "pid (comm) state "
    ppid[as.character(pp)] <- as.integer(strsplit(after, " ", fixed = TRUE)[[1]][1])
  }
  out <- root; frontier <- root
  repeat {
    kids <- as.integer(names(ppid)[ppid %in% frontier])
    kids <- setdiff(kids, out)
    if (!length(kids)) break
    out <- c(out, kids); frontier <- kids
  }
  out
}
.bench_tree_kb <- function(root) sum(vapply(.bench_descendants(root), .bench_pid_kb, numeric(1)))

# Function executed in the CHILD process (fresh R): load the package, run forde()
# `iters` times for one backend, return the elapsed seconds per iteration.
# dt_threads pins data.table (per worker); ranger_threads sets ranger's default
# threads, which forde()'s terminalNodes prediction picks up (it does not pass
# num.threads itself).
.bench_cell_fn <- function(pkgdir, data_path, backend, n_workers, dt_threads,
                           ranger_threads, iters) {
  suppressWarnings(suppressMessages(pkgload::load_all(pkgdir, quiet = TRUE)))
  data.table::setDTthreads(dt_threads)
  options(ranger.num.threads = ranger_threads)
  d <- readRDS(data_path); arf <- d$arf; X <- d$X
  if (backend == "sequential") {
    options(arf.backend = NULL); par <- FALSE
  } else if (backend == "foreach") {
    options(arf.backend = "foreach"); par <- TRUE
    doParallel::registerDoParallel(cores = n_workers)
  } else if (backend == "mirai") {
    options(arf.backend = "mirai"); par <- TRUE
    mirai::daemons(n_workers)
    mirai::everywhere(data.table::setDTthreads(dt_threads))
  } else {
    stop("unknown backend: ", backend)
  }
  secs <- vapply(seq_len(iters),
                 function(i) system.time(forde(arf, X, parallel = par))[["elapsed"]],
                 numeric(1))
  if (backend == "mirai") mirai::daemons(0)
  secs
}

# Parent-side: launch a cell in a fresh subprocess and sample its peak memory
# via /proc while it runs. Returns list(seconds = median, peak_mb).
bench_measure_cell <- function(backend, data_path, n_workers, dt_threads,
                               pkgdir, ranger_threads = 1L, iters = 1L,
                               interval = 0.05) {
  proc <- callr::r_bg(.bench_cell_fn,
                      args = list(pkgdir, data_path, backend, n_workers,
                                  dt_threads, ranger_threads, iters))
  pid <- proc$get_pid()
  peak_kb <- 0
  repeat {
    alive <- proc$is_alive()
    peak_kb <- max(peak_kb, .bench_tree_kb(pid))
    if (!alive) break
    Sys.sleep(interval)
  }
  secs <- tryCatch(proc$get_result(), error = function(e) NA_real_)
  list(seconds = stats::median(secs), peak_mb = peak_kb / 1024)
}
