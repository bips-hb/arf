# batchtools orchestrates; the cgroup sampler in bench-helpers.R measures.
# Slurm's own accounting (MaxRSS) is unusable for this comparison, which is why
# the metric never comes from the scheduler (bench/DESIGN.md).

bench_make_registry <- function(
  dir,
  cluster = c("local", "slurm"),
  slurm_template = Sys.getenv("ARF_BENCH_SLURM_TMPL", "")
) {
  cluster <- match.arg(cluster)
  reg <- batchtools::makeRegistry(
    file.dir = dir,
    source = c(
      "bench/bench-helpers.R",
      "bench/refs.R",
      "bench/cells.R",
      "bench/run-cell.R"
    ),
    packages = character(),
    make.default = FALSE,
    conf.file = NA_character_
  )
  reg$cluster.functions <- if (cluster == "slurm") {
    if (!nzchar(slurm_template)) {
      stop("set ARF_BENCH_SLURM_TMPL to the BIPS cluster slurm template", call. = FALSE)
    }
    batchtools::makeClusterFunctionsSlurm(template = slurm_template)
  } else {
    # Sequential and in-process on purpose: concurrent cells would contend for
    # CPU and memory and spoil both metrics. Cell isolation comes from the
    # callr child, not from the scheduler.
    batchtools::makeClusterFunctionsInteractive()
  }
  reg
}

# One job per (op, cell) with every ref inside it. Refs are installed INSIDE
# the job so each comparison uses libraries built on the node that measures it.
.bench_job <- function(i, cells, ref_specs) {
  # Computed HERE, not passed in: a lazy default evaluated in the orchestrator
  # bakes ITS tempdir into every job, so two slurm jobs on one node would
  # install into the same per-ref library concurrently and collide on
  # 00LOCK-arf. Under the interactive cluster functions the jobs run in-process,
  # so this is still the orchestrator's tempdir and the install cache is kept.
  lib_root <- file.path(tempdir(), "arf-bench-lib")
  cell <- cells[i, ]
  refs <- lapply(ref_specs, bench_install_ref, root = lib_root)
  # Must list every argument the run closures in .bench_cell_fn() actually
  # pass. adversarial_rf() has `...`, so an argument missing from an old ref is
  # absorbed and forwarded to ranger rather than erroring -- the one op where
  # this guard is load-bearing.
  calls <- list(
    forde = "parallel",
    lik = c("arf", "batch", "parallel"),
    forge = c("n_synth", "evidence", "stepsize", "evidence_row_mode", "parallel", "verbose"),
    expct = c("evidence", "stepsize", "evidence_row_mode", "parallel", "verbose"),
    adversarial_rf = c("num_trees", "parallel", "verbose")
  )
  for (ref in refs) {
    bench_assert_args(ref$lib, calls[cell$op])
  }
  bench_run_cell(cell, refs, job_id = Sys.getenv("SLURM_JOB_ID", NA_character_))
}

bench_submit_cells <- function(reg, cells, ref_specs, resources = list()) {
  batchtools::batchMap(
    .bench_job,
    i = seq_len(nrow(cells)),
    more.args = list(cells = cells, ref_specs = ref_specs),
    reg = reg
  )
  batchtools::submitJobs(resources = resources, reg = reg)
  reg
}

bench_collect <- function(reg) {
  done <- batchtools::findDone(reg = reg)
  if (!nrow(done)) {
    stop("no jobs finished; see batchtools::getErrorMessages()", call. = FALSE)
  }
  # A cell lost to a transient fetch failure, an OOM kill or an argument check
  # would otherwise vanish from the CSV, the report and history with no trace,
  # which is the half-populated report DESIGN forbids.
  err <- batchtools::findErrors(reg = reg)
  if (nrow(err)) {
    warning(
      nrow(err),
      " job(s) errored and are missing from these results: ids ",
      paste(err$job.id, collapse = ", "),
      "; see batchtools::getErrorMessages(reg = reg)",
      call. = FALSE
    )
  }
  rows <- batchtools::reduceResultsList(ids = done, reg = reg)
  out <- do.call(rbind, rows)
  out[, bench_schema()]
}
