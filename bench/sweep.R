#!/usr/bin/env Rscript
# Worker/size sweep for forde() backends -> tidy CSV in bench/results/.
# Measures wall-clock time and peak memory (PSS on Linux) for sequential /
# foreach / mirai across a grid of worker counts and problem sizes. The ARF is
# fit ONCE per (n, trees) and reused across all worker counts.
#
# Run on RESERVED cores; see bench/README.md for the topology rationale
# (data.table threads are pinned so "N workers" means N cores).
#
# Usage (env-overridable, comma-separated grids):
#   Rscript bench/sweep.R
#   ARF_BENCH_WORKERS=1,2,4,8,16 ARF_BENCH_N=5000,20000 ARF_BENCH_TREES=100,200 \
#     ARF_BENCH_DT_THREADS=1 Rscript bench/sweep.R

source("bench/bench-helpers.R")
bench_require_backends()

dt_threads  <- as.integer(Sys.getenv("ARF_BENCH_DT_THREADS", "1"))
data.table::setDTthreads(dt_threads)
worker_grid <- bench_ints("ARF_BENCH_WORKERS", c(1, 2, 4, 8))
n_grid      <- bench_ints("ARF_BENCH_N", c(5000, 20000))
trees_grid  <- bench_ints("ARF_BENCH_TREES", c(100, 200))
p           <- as.integer(Sys.getenv("ARF_BENCH_P", "30"))
backends    <- c("sequential", "foreach", "mirai")

message(sprintf(
  "Sweep | workers {%s} | n {%s} | trees {%s} | p %d | dt.threads/worker %d | metric %s",
  paste(worker_grid, collapse = ","), paste(n_grid, collapse = ","),
  paste(trees_grid, collapse = ","), p, dt_threads, BENCH_METRIC))

rows <- list()
for (n in n_grid) {
  for (trees in trees_grid) {
    set.seed(1)
    X <- bench_make_data(n, p)
    arf <- adversarial_rf(X, num_trees = trees, verbose = FALSE, parallel = FALSE)
    for (w in worker_grid) {
      for (be in backends) {
        # sequential is worker-independent: run it once (at the first worker count)
        if (be == "sequential" && w != worker_grid[1]) next
        m <- bench_run_backend(be, arf, X, w, dt_threads)
        rows[[length(rows) + 1L]] <- data.frame(
          n = n, trees = trees,
          workers = if (be == "sequential") NA_integer_ else w,
          dt_threads = dt_threads, backend = be,
          seconds = round(m$seconds, 2), peak_mb = round(m$peak_mb, 1),
          metric = BENCH_METRIC)
        message(sprintf("  n=%-6d trees=%-4d w=%-3s %-10s %8.1fs %9.1f MB",
                        n, trees, if (be == "sequential") "-" else as.character(w),
                        be, m$seconds, m$peak_mb))
      }
    }
  }
}

out <- do.call(rbind, rows)
dir.create("bench/results", showWarnings = FALSE, recursive = TRUE)
stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
f <- sprintf("bench/results/sweep-%s.csv", stamp)
write.csv(out, f, row.names = FALSE)
cat(sprintf("\n=== forde() backend sweep (%s) ===\n", BENCH_METRIC))
print(out, row.names = FALSE)
cat("\nWritten to", f, "\n")
