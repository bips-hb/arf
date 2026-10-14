# batchtools orchestrates; the cgroup sampler in bench-helpers.R measures.
# Slurm's own accounting (MaxRSS) is unusable for this comparison, which is why
# the metric never comes from the scheduler.

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
  req <- .bench_job_resources(cells, resources, n_refs = length(refs))
  # Submit in resource classes rather than once with a single figure: see
  # bench_cell_memory_mb(). batchtools takes resources per submitJobs call, so
  # one call per distinct request gets each cell what it needs without making
  # every cell queue for the largest.
  extra <- resources[setdiff(names(resources), c("ncpus", "memory", "walltime"))]
  grp <- paste(req$ncpus, req$memory, req$walltime)
  for (g in unique(grp[order(req$ncpus, req$memory, req$walltime)])) {
    i <- which(grp == g)
    message(sprintf(
      "  submitting %d cell(s): %d cpus, %.0f GB, %.1f h",
      length(i),
      req$ncpus[i[1]],
      req$memory[i[1]] / 1024,
      req$walltime[i[1]] / 3600
    ))
    batchtools::submitJobs(
      ids = data.table::data.table(job.id = i),
      resources = c(
        extra,
        list(
          ncpus = req$ncpus[i[1]],
          memory = req$memory[i[1]],
          walltime = req$walltime[i[1]]
        )
      ),
      reg = reg
    )
  }
  reg
}

# Per cell, not one figure for the grid: timed calls span two orders of
# magnitude, so a single walltime either reserves hours for a two-minute cell
# or kills the long ones. A requested dimension replaces only itself: cores and
# memory are independent asks, and dropping the derivation of the other one
# leaves it at the scheduler's default.
.bench_job_resources <- function(cells, resources = list(), n_refs) {
  req <- bench_cell_resources(cells)
  req$walltime <- if (is.null(resources$walltime)) {
    bench_cell_walltime(cells, n_refs = n_refs, fallback = 28800)
  } else {
    rep(as.integer(resources$walltime), nrow(cells))
  }
  for (nm in c("ncpus", "memory")) {
    if (!is.null(resources[[nm]])) {
      req[[nm]] <- rep(as.integer(resources[[nm]]), nrow(cells))
    }
  }
  req
}

# Newest registry unless one is named: the usual question is about the run just
# submitted, and both collect and resubmit ask it.
bench_registry_dir <- function(registry = NULL) {
  if (!is.null(registry)) {
    return(registry)
  }
  dirs <- sort(list.dirs("bench/registry", recursive = FALSE), decreasing = TRUE)
  if (!length(dirs)) {
    stop("no registry under bench/registry", call. = FALSE)
  }
  dirs[1]
}

# A job that expired got too little walltime, so resubmitting the same request
# repeats the same death. Escalation starts from the request the job actually
# carried, read back from the registry: the grid may have been edited since the
# run started, so recomputing it from bench_cells() can describe a cell that is
# no longer at that job id.
.bench_escalate_walltime <- function(walltime, factor) {
  w <- BENCH_WALLTIME_MIN_S * ceiling(factor * walltime / BENCH_WALLTIME_MIN_S)
  as.integer(pmin(BENCH_WALLTIME_MAX_S, pmax(BENCH_WALLTIME_MIN_S, w)))
}

# Resubmit the cells a run lost. Expiry has two causes and the registry cannot
# tell them apart: the scheduler killed the job at its walltime, or the kernel
# killed it for memory. Walltime grows by default and memory only when asked,
# since doubling both would double the node-hours spent to learn the same thing.
bench_resubmit <- function(
  reg,
  ids = NULL,
  factor = 2,
  mem_factor = 1,
  errors = FALSE,
  dry_run = FALSE
) {
  if (is.null(ids)) {
    ids <- batchtools::findExpired(reg = reg)$job.id
    if (isTRUE(errors)) {
      ids <- sort(union(ids, batchtools::findErrors(reg = reg)$job.id))
    }
  }
  req <- data.table::data.table(job.id = as.integer(ids))
  if (!nrow(req)) {
    message(
      "Nothing to resubmit: no expired",
      if (isTRUE(errors)) " or errored" else "",
      " jobs in this registry."
    )
    return(invisible(req))
  }
  res <- batchtools::getJobResources(ids = req$job.id, reg = reg)$resources
  req$walltime <- .bench_escalate_walltime(
    vapply(res, function(r) as.numeric(r$walltime %||% NA_real_), numeric(1)),
    factor
  )
  req$memory <- as.integer(
    mem_factor * vapply(res, function(r) as.numeric(r$memory %||% NA_real_), numeric(1))
  )
  req$ncpus <- vapply(res, function(r) as.integer(r$ncpus %||% NA_integer_), integer(1))
  # Name the cells, not just the ids: the job carries the grid it was built
  # from, so this says which measurement is being retried without trusting
  # that bench_cells() still returns the same rows in the same order.
  for (k in seq_len(nrow(req))) {
    cell <- tryCatch(
      {
        pars <- batchtools::makeJob(req$job.id[k], reg = reg)$pars
        pars$cells[pars$i, ]
      },
      error = function(e) NULL
    )
    message(sprintf(
      "  %4d  %-14s %-10s w=%-4s n=%-7s e=%-5s %-9s %.1f h -> %.1f h, %.0f GB",
      req$job.id[k],
      cell$op %||% "?",
      cell$backend %||% "?",
      cell$workers %||% "-",
      cell$n %||% "?",
      cell$n_evidence %||% "-",
      cell$rowmode %||% "-",
      as.numeric(res[[k]]$walltime %||% NA) / 3600,
      req$walltime[k] / 3600,
      req$memory[k] / 1024
    ))
  }
  if (isTRUE(dry_run)) {
    message(nrow(req), " job(s) would be resubmitted; --dry-run, so nothing was.")
    return(invisible(req))
  }
  # Per-job resource columns in `ids`, which submitJobs documents for walltime,
  # memory and ncpus: a lost cell keeps its own shape without one call per
  # distinct request.
  batchtools::submitJobs(ids = req, reg = reg)
  message(nrow(req), " job(s) resubmitted; collate with arf-bench collect when they finish.")
  invisible(req)
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
  # Expired jobs are neither done nor errored, so without this a walltime-killed
  # cell vanishes and a partial grid reads as a complete one.
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
