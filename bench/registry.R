# batchtools orchestrates; the cgroup sampler in bench-helpers.R measures.

bench_make_registry <- function(dir, cluster = c("local", "slurm")) {
  cluster <- match.arg(cluster)
  # batchtools asserts the parent exists rather than creating it.
  dir.create(dirname(dir), recursive = TRUE, showWarnings = FALSE)
  srcs <- c(
    "bench/bench-helpers.R",
    "bench/refs.R",
    "bench/cells.R",
    "bench/run-cell.R"
  )
  if (cluster == "slurm") {
    # The site config (/etc/xdg/batchtools/config.R) supplies the template,
    # partition, qos and max.concurrent.jobs.
    reg <- batchtools::makeRegistry(
      file.dir = dir,
      source = srcs,
      packages = character(),
      make.default = FALSE
    )
    if (
      is.null(reg$cluster.functions) ||
        identical(reg$cluster.functions$name, "Interactive")
    ) {
      stop("no slurm cluster functions: batchtools read no site config", call. = FALSE)
    }
    return(reg)
  }
  reg <- batchtools::makeRegistry(
    file.dir = dir,
    source = srcs,
    packages = character(),
    make.default = FALSE,
    conf.file = NA_character_
  )
  # Sequential and in-process: concurrent cells would contend for CPU and
  # memory. Isolation comes from the callr child, not the scheduler.
  reg$cluster.functions <- batchtools::makeClusterFunctionsInteractive()
  reg
}

# One job per cell, every ref measured inside it, so a comparison stays on one
# node. Refs arrive installed on the shared filesystem (arf is pure R).
.bench_job <- function(i, cells, refs) {
  cell <- cells[i, ]
  # Every argument the run closures in .bench_cell_fn() pass. adversarial_rf()
  # has `...`, so a missing argument there is absorbed rather than an error.
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
  # batchtools takes resources per submitJobs call: one call per distinct shape.
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

# Per cell; a forced dimension replaces only itself.
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

.bench_escalate_walltime <- function(walltime, factor) {
  w <- BENCH_WALLTIME_MIN_S * ceiling(factor * walltime / BENCH_WALLTIME_MIN_S)
  as.integer(pmin(BENCH_WALLTIME_MAX_S, pmax(BENCH_WALLTIME_MIN_S, w)))
}

# Expiry is a walltime kill or an OOM kill and the registry cannot tell them
# apart; walltime grows by default, memory only when asked. Requests are read
# back from the registry, since the grid may have been edited since.
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
  # Per-job resource columns in `ids`: each lost cell keeps its own shape.
  batchtools::submitJobs(ids = req, reg = reg)
  message(nrow(req), " job(s) resubmitted; collate with arf-bench collect when they finish.")
  invisible(req)
}

bench_collect <- function(reg) {
  done <- batchtools::findDone(reg = reg)
  if (!nrow(done)) {
    stop("no jobs finished; see batchtools::getErrorMessages()", call. = FALSE)
  }
  # Lost cells must not vanish silently: a partial grid would read as complete.
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
  # A job stores whatever the code produced when it ran, so a grid spanning a
  # schema edit holds both shapes; fill rather than refuse.
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
