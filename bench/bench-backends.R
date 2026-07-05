#!/usr/bin/env Rscript
# Benchmark forde() parallel backends across task sizes with the `bench`
# package. Results are written to bench/results/ (gitignored).
# See bench/README.md. Not part of the package build (.Rbuildignore ^bench$).
#
# Usage:
#   Rscript bench/bench-backends.R
#   ARF_BENCH_WORKERS=8 Rscript bench/bench-backends.R
#
# MEMORY CAVEAT: bench::mark's `mem_alloc` measures allocations in the MAIN R
# process only. The parallel backends do their per-tree work in separate
# worker/daemon processes, so mem_alloc does NOT capture the per-worker data
# copies that motivate mirai+mori (it can even make the parallel backends look
# lighter than sequential). Treat mem_alloc/gc as main-process signals; for the
# cross-process memory story use OS-level RSS/PSS of the process tree.

source("bench/bench-helpers.R")
stopifnot(requireNamespace("bench", quietly = TRUE))
bench_require_backends()

n_workers  <- as.integer(Sys.getenv("ARF_BENCH_WORKERS", "4"))
# data.table is multi-threaded by default, which oversubscribes cores when N
# workers each spin up threads (and even makes "sequential" multi-threaded).
# Pin threads per worker so "N workers" means N cores -- the fair comparison.
dt_threads <- as.integer(Sys.getenv("ARF_BENCH_DT_THREADS", "1"))
data.table::setDTthreads(dt_threads)

backends <- c("sequential", "foreach", "mirai")
message("Backends: ", paste(backends, collapse = ", "),
        "  |  workers = ", n_workers, "  |  data.table threads/worker = ", dt_threads)

# Register the parallel backends once; they persist for the whole run.
# doParallel forks inherit the main setDTthreads; mirai daemons need it set too.
doParallel::registerDoParallel(cores = n_workers)
mirai::daemons(n_workers)
mirai::everywhere(data.table::setDTthreads(dt_threads))

# One forde() call for a given backend; the option toggles the code path.
forde_be <- function(be, arf, X) {
  if (be == "sequential") return(forde(arf, X, parallel = FALSE))
  options(arf.backend = be)
  on.exit(options(arf.backend = NULL), add = TRUE)
  forde(arf, X, parallel = TRUE)
}

# Grid is env-overridable (comma-separated), e.g.
#   ARF_BENCH_N=500 ARF_BENCH_TREES=10 Rscript bench/bench-backends.R
grid <- expand.grid(
  n = bench_ints("ARF_BENCH_N", c(1000, 5000, 20000)),
  trees = bench_ints("ARF_BENCH_TREES", c(50, 200)),
  KEEP.OUT.ATTRS = FALSE
)

results <- bench::press(
  .grid = grid,
  {
    set.seed(1)
    X <- bench_make_data(n, p = 30L)
    arf <- adversarial_rf(X, num_trees = trees, verbose = FALSE, parallel = FALSE)
    exprs <- setNames(
      lapply(backends, function(be) bquote(forde_be(.(be), arf, X))),
      backends)
    # memory = FALSE is required: bench cannot profile memory for parallel
    # code (it errors otherwise), and per-process alloc would not capture the
    # cross-process copies anyway. This measures TIME; see README for memory.
    # check = tolerant all.equal doubles as a cross-backend correctness gate.
    bench::mark(
      exprs = exprs,
      check = function(a, b) isTRUE(all.equal(a, b, check.attributes = FALSE)),
      memory = FALSE, filter_gc = FALSE, iterations = 3
    )
  }
)

mirai::daemons(0)

# Persist to the gitignored results directory.
dir.create("bench/results", showWarnings = FALSE, recursive = TRUE)
stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
saveRDS(results, sprintf("bench/results/backends-%s.rds", stamp))
flat <- data.frame(
  n        = results$n,
  trees    = results$trees,
  backend  = as.character(results$expression),
  min_s    = as.numeric(results$min),
  median_s = as.numeric(results$median),
  itr_per_s = as.numeric(results$`itr/sec`)
)
write.csv(flat, sprintf("bench/results/backends-%s.csv", stamp), row.names = FALSE)
cat("\n=== forde() backend benchmark (bench::press) ===\n")
print(flat, row.names = FALSE)
cat("\nWritten to bench/results/backends-", stamp, ".{rds,csv}\n", sep = "")
