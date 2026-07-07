#!/usr/bin/env Rscript
# Pipeline sweep: sequential / foreach / mirai across the arf operations that
# now share the backend (forde, forge, expct, lik, adversarial_rf), reporting
# wall-clock and peak memory (PSS on Linux). Companion to sweep.R (forde-only);
# gives the fuller picture once forge/expct/lik/cforde/prune all use mirai/mori.
#
# Each cell runs in its OWN fresh subprocess (see bench/bench-helpers.R) -- mixing
# mirai with forked parallelism in one process crashes. Per (n, trees) the arf +
# psi + evidence are built ONCE, cached to disk; children read them. The
# adversarial_rf op is the exception -- it refits, since its cost IS train+prune.
#
# Run on RESERVED cores; see bench/README.md.
#
# Usage (env-overridable, comma-separated grids):
#   ARF_BENCH_OPS=forde,forge,expct,lik,adversarial_rf \
#   ARF_BENCH_WORKERS=1,2,4,8,16 ARF_BENCH_N=5000,20000 ARF_BENCH_TREES=100 \
#   ARF_BENCH_NEVIDENCE=100 ARF_BENCH_NSYNTH=1 ARF_BENCH_NFOLDS=8 \
#   ARF_BENCH_DT_THREADS=1 ARF_BENCH_ITERS=1 Rscript bench/sweep-ops.R

source("bench/bench-helpers.R")
bench_require_backends()

pkgdir      <- normalizePath(".")
dt_threads  <- as.integer(Sys.getenv("ARF_BENCH_DT_THREADS", "1"))
rgr_threads <- as.integer(Sys.getenv("ARF_BENCH_RANGER_THREADS", "1"))
iters       <- as.integer(Sys.getenv("ARF_BENCH_ITERS", "1"))
worker_grid <- bench_ints("ARF_BENCH_WORKERS", c(1, 2, 4, 8))
n_grid      <- bench_ints("ARF_BENCH_N", c(5000, 20000))
trees_grid  <- bench_ints("ARF_BENCH_TREES", 100)
p           <- as.integer(Sys.getenv("ARF_BENCH_P", "30"))
n_evidence  <- as.integer(Sys.getenv("ARF_BENCH_NEVIDENCE", "100"))
n_synth     <- as.integer(Sys.getenv("ARF_BENCH_NSYNTH", "1"))
n_folds     <- as.integer(Sys.getenv("ARF_BENCH_NFOLDS", "8"))
# evidence_row_mode for forge/expct: "separate" parallelizes forge/expct over
# steps; "or" delegates parallelism to cforde (benchmarks the cforde backend).
rowmode     <- match.arg(Sys.getenv("ARF_BENCH_ROWMODE", "separate"),
                         c("separate", "or"))
ops_grid    <- strsplit(Sys.getenv("ARF_BENCH_OPS",
                 "forde,forge,expct,lik,adversarial_rf"), ",")[[1]]
backends    <- c("sequential", "foreach", "mirai")

message(sprintf(
  "Ops sweep | ops {%s} | workers {%s} | n {%s} | trees {%s} | p %d | n_evidence %d | n_synth %d | n_folds %d | rowmode %s | metric %s",
  paste(ops_grid, collapse = ","), paste(worker_grid, collapse = ","),
  paste(n_grid, collapse = ","), paste(trees_grid, collapse = ","),
  p, n_evidence, n_synth, n_folds, rowmode, BENCH_METRIC))

cores <- as.integer(Sys.getenv("ARF_BENCH_CORES", parallel::detectCores()))
if (max(worker_grid) * dt_threads > cores) {
  message(sprintf(
    "  WARNING: peak workers x dt.threads = %d exceeds %d cores -- OVERSUBSCRIBED (mirai worst).",
    max(worker_grid) * dt_threads, cores))
}

# op-specific knobs (stepsize/batch left auto -> sized to worker count).
op_args_for <- function(op, trees, n) switch(op,
  forge          = list(n_synth = n_synth, stepsize = 0L, rowmode = rowmode),
  expct          = list(stepsize = 0L, rowmode = rowmode),
  lik            = list(batch = ceiling(n / n_folds)),
  adversarial_rf = list(trees = trees),
  list())

rows <- list()
for (n in n_grid) {
  for (trees in trees_grid) {
    set.seed(1)
    X <- bench_make_data(n, p)
    arf <- adversarial_rf(X, num_trees = trees, verbose = FALSE, parallel = FALSE)
    psi <- forde(arf, X, parallel = FALSE)
    evidence <- data.frame(grp = sample(levels(X$grp), n_evidence, replace = TRUE))
    data_path <- tempfile(fileext = ".rds")
    saveRDS(list(arf = arf, X = X, psi = psi, evidence = evidence), data_path)
    rm(arf, psi); invisible(gc())
    for (op in ops_grid) {
      oa <- op_args_for(op, trees, n)
      for (w in worker_grid) {
        for (be in backends) {
          if (be == "sequential" && w != worker_grid[1]) next
          m <- bench_measure_cell(be, data_path, w, dt_threads, pkgdir,
                                  ranger_threads = rgr_threads, iters = iters,
                                  op = op, op_args = oa)
          rows[[length(rows) + 1L]] <- data.frame(
            op = op, n = n, trees = trees,
            workers = if (be == "sequential") NA_integer_ else w,
            backend = be, seconds = round(m$seconds, 2),
            peak_mb = round(m$peak_mb, 1), metric = BENCH_METRIC,
            commit = BENCH_COMMIT,
            # op-scale knobs recorded per row (NA where an op ignores them)
            n_evidence = if (op %in% c("forge", "expct")) n_evidence else NA_integer_,
            n_synth = if (op == "forge") n_synth else NA_integer_,
            n_folds = if (op == "lik") n_folds else NA_integer_,
            rowmode = if (op %in% c("forge", "expct")) rowmode else NA_character_)
          message(sprintf("  %-14s n=%-6d trees=%-4d w=%-3s %-10s %8.2fs %9.1f MB",
                          op, n, trees,
                          if (be == "sequential") "-" else as.character(w),
                          be, m$seconds, m$peak_mb))
        }
      }
    }
    unlink(data_path)
  }
}

out <- do.call(rbind, rows)
dir.create("bench/results", showWarnings = FALSE, recursive = TRUE)
stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
# ARF_BENCH_LABEL (set per job by bench/submit-ops.sh) keys the filename to the
# config, so parallel jobs stay distinguishable and viz.R can pick the latest
# run per config.
label <- Sys.getenv("ARF_BENCH_LABEL", "")
f <- sprintf("bench/results/sweep-ops-%s%s.csv",
             if (nzchar(label)) paste0(label, "-") else "", stamp)
write.csv(out, f, row.names = FALSE)
cat(sprintf("\n=== arf pipeline backend sweep (%s) ===\n", BENCH_METRIC))
print(out, row.names = FALSE)
cat("\nWritten to", f, "\n")
