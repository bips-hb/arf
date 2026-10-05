# One cell measured against every ref, in one process, on one node.
# The ref loop lives HERE rather than in the batchtools grid on purpose: if
# `ref` were a job dimension, refs would scatter across nodes and reintroduce
# the drift that same-node A/B exists to remove (bench/DESIGN.md).

bench_schema <- function() {
  c(
    "ref",
    "arf_version",
    "commit",
    "op",
    "backend",
    "workers",
    "n",
    "p",
    "trees",
    "iters",
    "dt_threads",
    "ranger_threads",
    "metric",
    "peak_mb",
    "floor_mb",
    "peak_delta_mb",
    "mem_reps",
    "time_median",
    "time_min",
    "time_max",
    "digest",
    "digest_kind",
    "tier",
    "n_evidence",
    "n_synth",
    "n_folds",
    "rowmode",
    "host",
    "kernel",
    "r_version",
    "job_id",
    "timestamp"
  )
}

# Built once per cell so every ref reads the same file: a comparison must never
# differ by its input data. Built in a CHILD against an explicit library, since
# the orchestrator's own `library(arf)` may resolve to any installed version.
.bench_cell_data <- function(cell, lib) {
  path <- tempfile(fileext = ".rds")
  callr::r(
    function(lib, make_data, n, p, trees, n_evidence, path) {
      .libPaths(c(lib, .libPaths()))
      library(arf)
      set.seed(1)
      X <- make_data(n, p)
      a <- adversarial_rf(X, num_trees = trees, verbose = FALSE, parallel = FALSE)
      psi <- forde(a, X, parallel = FALSE)
      evidence <- if (!is.na(n_evidence)) {
        data.frame(grp = sample(levels(X$grp), n_evidence, replace = TRUE))
      } else {
        NULL
      }
      saveRDS(list(arf = a, X = X, psi = psi, evidence = evidence), path)
    },
    args = list(
      lib = lib,
      make_data = bench_make_data,
      n = cell$n,
      p = cell$p,
      trees = cell$trees,
      n_evidence = cell$n_evidence,
      path = path
    )
  )
  path
}

.bench_op_args <- function(cell) {
  switch(
    cell$op,
    forge = list(n_synth = cell$n_synth, stepsize = 0L, rowmode = cell$rowmode),
    expct = list(stepsize = 0L, rowmode = cell$rowmode),
    lik = list(batch = ceiling(cell$n / cell$n_folds)),
    adversarial_rf = list(trees = cell$trees),
    list()
  )
}

bench_run_cell <- function(
  cell,
  refs,
  data_path = NULL,
  baseline = "main",
  dt_threads = 1L,
  ranger_threads = 1L,
  job_id = NA_character_,
  mem_reps = as.integer(Sys.getenv("ARF_BENCH_MEM_REPS", "1"))
) {
  stopifnot(nrow(cell) == 1L)
  labels <- vapply(refs, function(r) r$label, character(1))
  fixture <- refs[[match(baseline, labels, nomatch = 1L)]]
  if (is.null(data_path)) {
    data_path <- .bench_cell_data(cell, fixture$lib)
    on.exit(unlink(data_path), add = TRUE)
  }
  uname <- tryCatch(system2("uname", "-r", stdout = TRUE), error = function(e) NA_character_)
  rows <- lapply(refs, function(ref) {
    # A ref that dies, from an OOM kill or a missing library, must not cost the
    # other refs their measurements.
    m <- tryCatch(
      bench_measure_cell(
        cell$backend,
        data_path,
        if (is.na(cell$workers)) NA_integer_ else cell$workers,
        dt_threads,
        ref$lib,
        ranger_threads = ranger_threads,
        iters = cell$iters,
        op = cell$op,
        op_args = .bench_op_args(cell),
        mem_reps = mem_reps
      ),
      error = function(e) {
        message("  ref ", ref$label, " failed: ", conditionMessage(e))
        list(
          seconds = NA_real_,
          peak_mb = NA_real_,
          digest = NA_character_,
          digest_kind = NA_character_,
          arf_version = ref$arf_version
        )
      }
    )
    secs <- m$seconds[!is.na(m$seconds)]
    # The floor is this cell's cost with the op never run. Deltas are taken on
    # the marginal, since the floor dwarfs small ops and would dilute them.
    floor_mb <- if (is.na(m$peak_mb)) NA_real_ else bench_measure_floor(ref$lib, data_path)
    data.frame(
      ref = ref$label,
      arf_version = ref$arf_version,
      commit = ref$commit,
      op = cell$op,
      backend = cell$backend,
      workers = cell$workers,
      n = cell$n,
      p = cell$p,
      trees = cell$trees,
      iters = cell$iters,
      dt_threads = dt_threads,
      ranger_threads = ranger_threads,
      metric = BENCH_METRIC,
      peak_mb = round(m$peak_mb, 1),
      floor_mb = round(floor_mb, 1),
      mem_reps = mem_reps,
      peak_delta_mb = round(m$peak_mb - floor_mb, 1),
      time_median = if (length(secs)) round(stats::median(secs), 3) else NA_real_,
      time_min = if (length(secs)) round(min(secs), 3) else NA_real_,
      time_max = if (length(secs)) round(max(secs), 3) else NA_real_,
      digest = m$digest,
      digest_kind = m$digest_kind,
      tier = cell$tier,
      n_evidence = cell$n_evidence,
      n_synth = cell$n_synth,
      n_folds = cell$n_folds,
      rowmode = cell$rowmode,
      host = Sys.info()[["nodename"]],
      kernel = uname[1],
      r_version = paste0(R.version$major, ".", R.version$minor),
      job_id = job_id,
      timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[, bench_schema()]
}
