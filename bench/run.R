#!/usr/bin/env Rscript
# Entry point for the version-comparing benchmark. See bench/DESIGN.md.
#
#   Rscript bench/run.R                       # quick tier, local
#   ARF_BENCH_TIER=full ARF_BENCH_CLUSTER=slurm Rscript bench/run.R
#
# Refs default to HEAD, main and the curated anchors.

source("bench/bench-helpers.R")
source("bench/refs.R")
source("bench/cells.R")
source("bench/run-cell.R")
source("bench/registry.R")
source("bench/collate.R")

bench_assert_deps()
bench_assert_clean_tree()

tier <- Sys.getenv("ARF_BENCH_TIER", "quick")
cluster <- Sys.getenv("ARF_BENCH_CLUSTER", "local")
refs <- unique(c("HEAD", "main", bench_anchors()))
refs <- bench_resolve_refs(refs[nzchar(refs)])

cells <- bench_cells(tier)
message(sprintf(
  "arf bench | tier %s | %d cells | refs %s | metric %s",
  tier,
  nrow(cells),
  paste(refs, collapse = ", "),
  BENCH_METRIC
))

dir <- file.path("bench", "registry", format(Sys.time(), "%Y%m%d-%H%M%S"))
reg <- bench_make_registry(dir, cluster)
# DESIGN "Settled" names these as the starting point. Submitting a 16-worker
# cell with the template's defaults would either cpu-throttle it, destroying
# the timing, or get it OOM-killed, destroying the cell.
resources <- if (identical(cluster, "slurm")) {
  list(ncpus = 17L, memory = 256000L, walltime = 8L * 3600L)
} else {
  list()
}
reg <- bench_submit_cells(reg, cells, refs, resources = resources)
batchtools::waitForJobs(reg = reg)

rows <- bench_collect(reg)
dir.create("bench/results", showWarnings = FALSE, recursive = TRUE)
out <- file.path(
  "bench/results",
  sprintf("bench-%s-%s.csv", tier, format(Sys.time(), "%Y%m%d-%H%M%S"))
)
write.csv(rows, out, row.names = FALSE)
message("Written to ", out)

# Anchor rows are appended automatically; only the commit is manual.
bench_history_append(rows)

report <- bench_report(rows)
cat(report, sep = "\n")
writeLines(report, "bench/results/report.md")
message("Report at bench/results/report.md")
