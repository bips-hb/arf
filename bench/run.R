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
# Staging filter: validate the slurm path on a slice before committing the
# whole grid (the full tier is 352 cells).
ops <- Sys.getenv("ARF_BENCH_OPS", "")
if (nzchar(ops)) {
  cells <- cells[cells$op %in% strsplit(ops, ",")[[1]], ]
}
max_cells <- Sys.getenv("ARF_BENCH_MAX_CELLS", "")
if (nzchar(max_cells)) {
  cells <- utils::head(cells, as.integer(max_cells))
}
message(sprintf(
  "arf bench | tier %s | %d cells | refs %s | metric %s",
  tier,
  nrow(cells),
  paste(refs, collapse = ", "),
  BENCH_METRIC
))

# Installed ONCE, here, into a library on the shared filesystem rather than in
# each job: see the note in .bench_job(). tempdir() would be node-local and
# invisible to the compute nodes, so this lives in the repo (bench/lib is
# gitignored).
lib_root <- Sys.getenv("ARF_BENCH_LIB_ROOT", file.path("bench", "lib"))
dir.create(lib_root, recursive = TRUE, showWarnings = FALSE)
message("installing refs into ", normalizePath(lib_root), " ...")
installed <- lapply(refs, bench_install_ref, root = normalizePath(lib_root))
for (r in installed) {
  message("  ", r$label, " -> arf ", r$arf_version, " (", r$commit, ")")
}

dir <- file.path("bench", "registry", format(Sys.time(), "%Y%m%d-%H%M%S"))
reg <- bench_make_registry(dir, cluster)
# batchtools' slurm templates take `memory` as megabytes PER CPU
# (#SBATCH --mem-per-cpu), not per job as submit-ops.sh's --mem=256G did.
# Passing 256000 with ncpus=17 would request 4.3 TB and never schedule.
# Check against the template in use: these names and that semantics are
# template-specific, which is why every one is env-overridable.
resources <- if (identical(cluster, "slurm")) {
  # ncpus counts HYPERTHREADS on this cluster (its config notes 1 physical core
  # = 2 threads), so a 16-worker cell asking for 17 would get 8.5 cores and be
  # oversubscribed, which degrades mirai catastrophically. Two threads per
  # worker plus one for the orchestrator.
  peak_workers <- suppressWarnings(max(cells$workers, na.rm = TRUE))
  if (!is.finite(peak_workers)) {
    peak_workers <- 1L
  }
  list(
    ncpus = as.integer(Sys.getenv(
      "ARF_BENCH_SLURM_CPUS",
      as.character(2L * (as.integer(peak_workers) + 1L))
    )),
    # TOTAL megabytes, via --mem. Deliberately not mem_per_cpu: the site
    # default.resources already sets `memory`, and the template rejects both
    # together.
    memory = as.integer(Sys.getenv("ARF_BENCH_SLURM_MEM", "256000")),
    # 8h = 480 min, inside the default "medium" QoS ceiling of 1440 min.
    walltime = as.integer(Sys.getenv("ARF_BENCH_SLURM_WALLTIME", as.character(8L * 3600L)))
  )
} else {
  list()
}
if (nzchar(Sys.getenv("ARF_BENCH_SLURM_PARTITION"))) {
  resources$partition <- Sys.getenv("ARF_BENCH_SLURM_PARTITION")
}
message(sprintf(
  "submitting %d cells%s",
  nrow(cells),
  if (identical(cluster, "slurm")) {
    sprintf(
      " | slurm: %d cpus (~%d cores), %.0f GB/job total, %.1f h",
      resources$ncpus,
      resources$ncpus %/% 2L,
      resources$memory / 1024,
      resources$walltime / 3600
    )
  } else {
    ""
  }
))
reg <- bench_submit_cells(reg, cells, installed, resources = resources)
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
