#!/usr/bin/env Rscript
# Single-config peak-memory + time snapshot for forde() backends. For scaling
# across worker counts and problem sizes, use bench/sweep.R instead.
# Each backend runs in its own subprocess (see bench/bench-helpers.R).
# See bench/README.md (esp. the PSS metric and topology notes).
#
# Usage:
#   Rscript bench/mem-backends.R
#   ARF_BENCH_N=20000 ARF_BENCH_TREES=200 ARF_BENCH_WORKERS=8 Rscript bench/mem-backends.R

source("bench/bench-helpers.R")
bench_require_backends()

pkgdir     <- normalizePath(".")
n          <- as.integer(Sys.getenv("ARF_BENCH_N", "5000"))
p          <- as.integer(Sys.getenv("ARF_BENCH_P", "30"))
trees      <- as.integer(Sys.getenv("ARF_BENCH_TREES", "100"))
n_workers  <- as.integer(Sys.getenv("ARF_BENCH_WORKERS", "4"))
dt_threads <- as.integer(Sys.getenv("ARF_BENCH_DT_THREADS", "1"))
rgr_threads <- as.integer(Sys.getenv("ARF_BENCH_RANGER_THREADS", "1"))
iters      <- as.integer(Sys.getenv("ARF_BENCH_ITERS", "1"))

set.seed(1)
X <- bench_make_data(n, p)
message(sprintf("Data %d x %d | trees %d | workers %d | dt.threads/worker %d | metric %s",
                nrow(X), ncol(X), trees, n_workers, dt_threads, BENCH_METRIC))
arf <- adversarial_rf(X, num_trees = trees, verbose = FALSE, parallel = FALSE)
data_path <- tempfile(fileext = ".rds")
saveRDS(list(arf = arf, X = X), data_path)
rm(arf, X); invisible(gc())

out <- do.call(rbind, lapply(c("sequential", "foreach", "mirai"), function(be) {
  m <- bench_measure_cell(be, data_path, n_workers, dt_threads, pkgdir,
                          ranger_threads = rgr_threads, iters = iters)
  data.frame(backend = be, seconds = round(m$seconds, 2),
             peak_mb = round(m$peak_mb, 1), metric = BENCH_METRIC,
             commit = BENCH_COMMIT)
}))
unlink(data_path)
cat(sprintf("\n=== forde() backend peak memory (%s) ===\n", BENCH_METRIC))
print(out, row.names = FALSE)

dir.create("bench/results", showWarnings = FALSE, recursive = TRUE)
stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
write.csv(cbind(n = n, trees = trees, workers = n_workers, dt_threads = dt_threads,
                ranger_threads = rgr_threads, out),
          sprintf("bench/results/mem-%s.csv", stamp), row.names = FALSE)
cat("\nWritten to bench/results/mem-", stamp, ".csv\n", sep = "")
