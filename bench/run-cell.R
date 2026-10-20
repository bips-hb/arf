# One cell measured against every ref, in one process, on one node. The ref
# loop lives here, not in the batchtools grid: refs on different nodes would
# reintroduce the drift that same-node A/B exists to remove.

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
    "mem_limit_mb",
    "time_median",
    "time_min",
    "time_max",
    "digest",
    "digest_kind",
    "digest_stable",
    "out_rows",
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

# One fixture per cell, so every ref reads the same file. Built in a child
# against an explicit library, since the session's arf may be any version.
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
  job_id = NA_character_
) {
  stopifnot(nrow(cell) == 1L)
  mem_reps <- if (is.null(cell$mem_reps) || is.na(cell$mem_reps)) 1L else cell$mem_reps
  # Before the fixture is built, so a cell with nothing to measure costs nothing.
  refs <- .bench_usable_refs(refs, cell$backend)
  if (!length(refs)) {
    return(NULL)
  }
  labels <- vapply(refs, function(r) r$label, character(1))
  fixture <- refs[[match(baseline, labels, nomatch = 1L)]]
  if (is.null(data_path)) {
    data_path <- .bench_cell_data(cell, fixture$lib)
    on.exit(unlink(data_path), add = TRUE)
  }
  uname <- tryCatch(system2("uname", "-r", stdout = TRUE), error = function(e) NA_character_)

  # Rounds interleave refs and alternate their order, so drift over the cell's
  # lifetime cancels instead of landing on one ref.
  rounds <- lapply(seq_len(mem_reps), function(round) {
    order <- if (round %% 2L == 0L) rev(seq_along(refs)) else seq_along(refs)
    measured <- lapply(order, function(j) {
      ref <- refs[[j]]
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
          op_args = .bench_op_args(cell)
        ),
        error = function(e) {
          message("  ref ", ref$label, " failed: ", conditionMessage(e))
          list(
            seconds = NA_real_,
            peak_mb = NA_real_,
            digest = NA_character_,
            digest_kind = NA_character_,
            out_rows = NA_integer_,
            arf_version = ref$arf_version
          )
        }
      )
      # The floor is the peak's difference partner, so it is measured in the same round.
      m$floor_mb <- if (is.na(m$peak_mb)) NA_real_ else bench_measure_floor(ref$lib, data_path)
      # One line per measurement, or a slurm log stays silent for an hour.
      message(sprintf(
        "  %s n=%g trees=%g %s w=%s | round %d/%d %-12s %6.2fs peak %6.1f floor %6.1f MB",
        cell$op,
        cell$n,
        cell$trees,
        cell$backend,
        if (is.na(cell$workers)) "-" else cell$workers,
        round,
        mem_reps,
        ref$label,
        if (length(m$seconds)) stats::median(m$seconds, na.rm = TRUE) else NA_real_,
        m$peak_mb,
        m$floor_mb
      ))
      m
    })
    measured[order(order)]
  })

  rows <- lapply(seq_along(refs), function(j) {
    ref <- refs[[j]]
    per_round <- lapply(rounds, function(r) r[[j]])
    mid <- function(v) if (all(is.na(v))) NA_real_ else stats::median(v, na.rm = TRUE)
    field <- function(nm) vapply(per_round, function(x) x[[nm]], numeric(1))
    first <- per_round[[1]]
    # A ref whose output moves between identical rounds (parallel forge
    # ignores the seed) cannot support a cross-ref output comparison.
    digests <- stats::na.omit(vapply(
      per_round,
      function(x) x$digest %||% NA_character_,
      character(1)
    ))
    seconds <- unlist(lapply(per_round, function(x) x$seconds))
    floor_mb <- mid(field("floor_mb"))
    # The marginal is a paired difference per round, then the median: a round
    # whose baseline was still contaminated yields floor > peak and is dropped.
    paired <- field("peak_mb") - field("floor_mb")
    paired <- paired[!is.na(paired) & paired > 0]
    marginal_mb <- if (length(paired)) stats::median(paired) else NA_real_
    secs <- seconds[!is.na(seconds)]
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
      peak_mb = round(mid(field("peak_mb")), 1),
      floor_mb = round(floor_mb, 1),
      mem_reps = mem_reps,
      mem_limit_mb = round(BENCH_MEM_LIMIT_MB, 1),
      peak_delta_mb = round(marginal_mb, 1),
      time_median = if (length(secs)) round(stats::median(secs), 3) else NA_real_,
      time_min = if (length(secs)) round(min(secs), 3) else NA_real_,
      time_max = if (length(secs)) round(max(secs), 3) else NA_real_,
      digest = first$digest,
      digest_kind = first$digest_kind,
      digest_stable = length(digests) > 0L && length(unique(digests)) == 1L,
      out_rows = first$out_rows %||% NA_integer_,
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
