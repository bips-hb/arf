#!/usr/bin/env Rscript
# Collate a run that was submitted without waiting:
#
#   Rscript bench/collect.R                      # newest registry under bench/registry
#   Rscript bench/collect.R bench/registry/2026... # a specific one
#
# Safe to run while jobs are still going: it collates whatever has finished and
# says how many are outstanding.

source("bench/bench-helpers.R")
source("bench/refs.R")
source("bench/cells.R")
source("bench/run-cell.R")
source("bench/registry.R")
source("bench/collate.R")

args <- commandArgs(trailingOnly = TRUE)
dir <- if (length(args)) {
  args[1]
} else {
  dirs <- sort(list.dirs("bench/registry", recursive = FALSE), decreasing = TRUE)
  if (!length(dirs)) {
    stop("no registry under bench/registry", call. = FALSE)
  }
  dirs[1]
}

# writeable = FALSE: another session may still own this registry.
reg <- batchtools::loadRegistry(dir, writeable = FALSE)
st <- batchtools::getStatus(reg = reg)
print(st)

rows <- bench_collect(reg)
dir.create("bench/results", showWarnings = FALSE, recursive = TRUE)
stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
out <- file.path("bench/results", sprintf("bench-%s.csv", stamp))
write.csv(rows, out, row.names = FALSE)
message("Written to ", out)

bench_history_append(rows)
report <- bench_report(rows)
cat(report, sep = "\n")
writeLines(report, "bench/results/report.md")
message("Report at bench/results/report.md")
