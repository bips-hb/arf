#!/usr/bin/env Rscript
# Peak-MEMORY benchmark for forde() backends. Companion to bench-backends.R,
# which measures TIME (the `bench` package cannot profile memory for parallel
# code). See bench/README.md. Not part of the package build.
#
# Metric: peak total PSS (Proportional Set Size) across this user's R process
# tree while forde() runs. PSS splits shared pages across the processes that map
# them, so mori's shared data is counted once (split), not once per daemon --
# the fair way to compare a copy-per-worker backend (foreach) against a
# shared-memory one (mirai+mori). PSS is Linux-only (/proc/<pid>/smaps_rollup);
# on other platforms the script falls back to RSS and says so.
#
# Usage:
#   Rscript bench/mem-backends.R
#   ARF_BENCH_N=20000 ARF_BENCH_TREES=200 ARF_BENCH_WORKERS=8 Rscript bench/mem-backends.R
#
# Run in an otherwise-idle session: the sampler sums over ALL of this user's R
# processes, so concurrent unrelated R work would be attributed here.

suppressWarnings(suppressMessages(pkgload::load_all(quiet = TRUE)))
stopifnot(requireNamespace("parallel", quietly = TRUE))

n         <- as.integer(Sys.getenv("ARF_BENCH_N", "5000"))
p         <- as.integer(Sys.getenv("ARF_BENCH_P", "30"))
trees     <- as.integer(Sys.getenv("ARF_BENCH_TREES", "100"))
n_workers <- as.integer(Sys.getenv("ARF_BENCH_WORKERS", "4"))

## ---- memory sampling: peak total PSS (or RSS) over our R process tree -------
use_pss <- file.exists("/proc/self/smaps_rollup")
metric  <- if (use_pss) "PSS" else "RSS"

r_proc_pids <- function() {
  pids <- suppressWarnings(as.integer(list.files("/proc")))
  pids <- pids[!is.na(pids)]
  keep <- vapply(pids, function(pp) {
    comm <- tryCatch(readLines(sprintf("/proc/%d/comm", pp), warn = FALSE),
                     error = function(e) "")
    # readable + an R process (R, Rscript, /usr/lib/R/bin/exec/R, ...)
    length(comm) && grepl("^R", comm)
  }, logical(1))
  pids[keep]
}
pid_mem_kb <- function(pp) {
  if (use_pss) {
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
total_mem_kb <- function() sum(vapply(r_proc_pids(), pid_mem_kb, numeric(1)))

with_peak_mb <- function(expr, interval = 0.02) {
  stopf <- tempfile()
  sampler <- parallel::mcparallel({
    peak <- 0
    repeat {
      peak <- max(peak, total_mem_kb())
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

## ---- backends ---------------------------------------------------------------
have_foreach <- requireNamespace("doParallel", quietly = TRUE)
have_mirai   <- requireNamespace("mirai", quietly = TRUE) &&
  requireNamespace("mori", quietly = TRUE)

set.seed(1)
X <- as.data.frame(matrix(stats::rnorm(n * p), n, p))
X$grp <- factor(sample(letters[1:6], n, replace = TRUE))
message(sprintf("Data %d x %d | trees %d | workers %d | metric %s",
                nrow(X), ncol(X), trees, n_workers, metric))
arf <- adversarial_rf(X, num_trees = trees, verbose = FALSE, parallel = FALSE)

bench1 <- function(label, setup = function() NULL, teardown = function() NULL) {
  ok <- tryCatch({ setup(); TRUE },
                 error = function(e) { message("  [skip] ", label, ": ",
                                                conditionMessage(e)); FALSE })
  if (!ok) return(NULL)
  on.exit(teardown(), add = TRUE)
  m <- with_peak_mb(forde(arf, X, parallel = !identical(label, "sequential")))
  data.frame(backend = label, seconds = round(m$seconds, 2),
             peak_mb = round(m$peak_mb, 1))
}

rows <- list()
rows$seq <- bench1("sequential", function() options(arf.backend = NULL))
if (have_foreach) {
  rows$fe <- bench1("foreach",
    function() { options(arf.backend = "foreach")
                 doParallel::registerDoParallel(cores = n_workers) },
    function() try(doParallel::stopImplicitCluster(), silent = TRUE))
} else message("  [skip] foreach: doParallel not installed")
if (have_mirai) {
  rows$mi <- bench1("mirai",
    function() { options(arf.backend = "mirai"); mirai::daemons(n_workers) },
    function() mirai::daemons(0))
} else message("  [skip] mirai: mirai and/or mori not installed")

out <- do.call(rbind, rows)
out$metric <- metric
cat(sprintf("\n=== forde() backend peak memory (%s) ===\n", metric))
print(out, row.names = FALSE)

dir.create("bench/results", showWarnings = FALSE, recursive = TRUE)
stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
write.csv(out, sprintf("bench/results/mem-%s.csv", stamp), row.names = FALSE)
cat("\nWritten to bench/results/mem-", stamp, ".csv\n", sep = "")
