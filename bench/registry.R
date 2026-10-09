# batchtools orchestrates; the cgroup sampler in bench-helpers.R measures.
# Slurm's own accounting (MaxRSS) is unusable for this comparison, which is why
# the metric never comes from the scheduler (bench/DESIGN.md).

bench_make_registry <- function(
  dir,
  cluster = c("local", "slurm"),
  slurm_template = Sys.getenv("ARF_BENCH_SLURM_TMPL", "")
) {
  cluster <- match.arg(cluster)
  # batchtools asserts the parent exists rather than creating it, and the
  # registry path is nested under bench/registry/<stamp>.
  dir.create(dirname(dir), recursive = TRUE, showWarnings = FALSE)
  srcs <- c(
    "bench/bench-helpers.R",
    "bench/refs.R",
    "bench/cells.R",
    "bench/run-cell.R"
  )
  # On the cluster, READ the system config: /etc/xdg/batchtools/config.R already
  # supplies cluster.functions with the site template plus default.resources
  # (qos, clusters, partition) and max.concurrent.jobs. Suppressing it with
  # conf.file = NA, which is right for the local tier's determinism, would throw
  # all of that away and demand a template we do not need to name.
  reg <- if (cluster == "slurm") {
    batchtools::makeRegistry(
      file.dir = dir,
      source = srcs,
      packages = character(),
      make.default = FALSE
    )
  } else {
    batchtools::makeRegistry(
      file.dir = dir,
      source = srcs,
      packages = character(),
      make.default = FALSE,
      conf.file = NA_character_
    )
  }
  if (cluster == "slurm") {
    # Only override what the site config set if a template is named explicitly.
    if (nzchar(slurm_template)) {
      reg$cluster.functions <- batchtools::makeClusterFunctionsSlurm(
        template = slurm_template,
        array.jobs = TRUE
      )
    }
    if (
      is.null(reg$cluster.functions) ||
        identical(reg$cluster.functions$name, "Interactive")
    ) {
      stop(
        "no slurm cluster functions found: batchtools read no site config and ",
        "no ARF_BENCH_SLURM_TMPL was given",
        call. = FALSE
      )
    }
    return(reg)
  }
  # Sequential and in-process on purpose: concurrent cells would contend for
  # CPU and memory and spoil both metrics. Cell isolation comes from the callr
  # child, not from the scheduler.
  reg$cluster.functions <- batchtools::makeClusterFunctionsInteractive()
  reg
}

# One job per (op, cell), with every ref measured inside it so each comparison
# stays on one node.
.bench_job <- function(i, cells, refs) {
  # Refs arrive already installed, from the orchestrator, into a library on the
  # shared filesystem. Installing inside the job instead would mean every job
  # running `git worktree add` against the one shared .git -- up to 352
  # concurrent mutations of .git/worktrees -- and 1056 redundant installs.
  # Safe here because arf is pure R: no src/, no NeedsCompilation, so the
  # installed tree is architecture-independent. A package with compiled code,
  # or heterogeneous nodes, would have to go back to installing per job.
  cell <- cells[i, ]
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

bench_submit_cells <- function(reg, cells, refs, resources = list()) {
  batchtools::batchMap(
    .bench_job,
    i = seq_len(nrow(cells)),
    more.args = list(cells = cells, refs = refs),
    reg = reg
  )
  # Submit in memory classes rather than once with a single figure: see
  # bench_cell_memory_mb(). batchtools takes resources per submitJobs call, so
  # one call per distinct request gets each cell what it needs without making
  # every cell queue for the largest.
  if (is.null(resources$memory) && is.null(resources$ncpus)) {
    req <- bench_cell_resources(cells)
    grp <- paste(req$ncpus, req$memory)
    for (g in unique(grp[order(req$fraction)])) {
      i <- which(grp == g)
      ids <- data.table::data.table(job.id = i)
      message(sprintf(
        "  submitting %d cell(s) at 1/%g node: %d cpus, %.0f GB",
        length(i),
        1 / req$fraction[i[1]],
        req$ncpus[i[1]],
        req$memory[i[1]] / 1024
      ))
      batchtools::submitJobs(
        ids = ids,
        resources = c(resources, list(ncpus = req$ncpus[i[1]], memory = req$memory[i[1]])),
        reg = reg
      )
    }
  } else {
    batchtools::submitJobs(resources = resources, reg = reg)
  }
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
  # Expired jobs are neither done nor errored, so they would otherwise leave no
  # trace: a walltime-killed cell simply vanishes and a 45-of-64 report reads
  # as a complete grid. A full run lost 19 cells to an 8 h walltime this way.
  missing <- setdiff(
    batchtools::findJobs(reg = reg)$job.id,
    c(done$job.id, err$job.id)
  )
  if (length(missing)) {
    warning(
      length(missing),
      " job(s) neither finished nor errored (expired, killed, or still ",
      "running) and are missing from these results: ids ",
      paste(missing, collapse = ", "),
      "; batchtools::findExpired(reg = reg) separates the dead from the live",
      call. = FALSE
    )
  }
  rows <- batchtools::reduceResultsList(ids = done, reg = reg)
  # Tolerate results that predate a schema change. A job stores whatever the
  # code produced when it RAN, so a grid running for days across an edit to
  # bench_schema() holds both shapes, and a strict column select would make the
  # whole run uncollectable at precisely the moment its results matter. Fill
  # and say so instead.
  out <- as.data.frame(data.table::rbindlist(rows, fill = TRUE, use.names = TRUE))
  missing <- setdiff(bench_schema(), names(out))
  if (length(missing)) {
    warning(
      "these results predate the current schema; filling with NA: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
    for (m in missing) {
      out[[m]] <- NA
    }
  }
  out[, bench_schema()]
}
