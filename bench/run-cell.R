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
  mem_reps = NULL
) {
  stopifnot(nrow(cell) == 1L)
  # Tier default, overridable per call, with the env var winning for one-offs.
  env_reps <- Sys.getenv("ARF_BENCH_MEM_REPS", "")
  mem_reps <- if (nzchar(env_reps)) {
    as.integer(env_reps)
  } else if (!is.null(mem_reps)) {
    mem_reps
  } else if (!is.null(cell$mem_reps) && !is.na(cell$mem_reps)) {
    cell$mem_reps
  } else {
    1L
  }
  labels <- vapply(refs, function(r) r$label, character(1))
  fixture <- refs[[match(baseline, labels, nomatch = 1L)]]
  if (is.null(data_path)) {
    data_path <- .bench_cell_data(cell, fixture$lib)
    on.exit(unlink(data_path), add = TRUE)
  }
  uname <- tryCatch(system2("uname", "-r", stdout = TRUE), error = function(e) NA_character_)

  # Interleaved, not blocked: measuring all of ref A's replicates before all of
  # ref B's lets any drift over the cell's lifetime land entirely on B. Measured
  # same-commit, blocked replication left a recurring one-sided outlier of
  # -24% to -42% at n=1e3 that more replicates did not shrink. One round per
  # replicate, every ref inside it, cancels drift that is linear in time.
  rounds <- lapply(seq_len(mem_reps), function(round) {
    # Counterbalanced: interleaving rounds alone leaves the within-round order
    # intact, so ref A is always measured before ref B and any drift inside a
    # round still lands one-sided. Alternate the order and the bias cancels
    # across rounds instead of accumulating.
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
          op_args = .bench_op_args(cell),
          mem_reps = 1L
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
      # The floor is a difference partner of the peak, so it is replicated in
      # the same rounds: at small n the two are nearly equal and un-replicated
      # floor noise lands undiluted in peak_mb - floor_mb.
      m$floor_mb <- if (is.na(m$peak_mb)) NA_real_ else bench_measure_floor(ref$lib, data_path)
      m
    })
    measured[order(order)]
  })

  rows <- lapply(seq_along(refs), function(j) {
    ref <- refs[[j]]
    per_round <- lapply(rounds, function(r) r[[j]])
    # The marginal is a PAIRED difference: peak and floor are measured in the
    # same round, under the same conditions, and subtracted there. Taking a
    # minimum on each end separately instead (min(peak) - min(floor)) inflates
    # the spread rather than reducing it, which is what measurement showed.
    # The median over rounds is then robust to the heavy tail that GC timing
    # gives peak memory.
    mid <- function(v) if (all(is.na(v))) NA_real_ else stats::median(v, na.rm = TRUE)
    field <- function(nm) vapply(per_round, function(x) x[[nm]], numeric(1))
    first <- per_round[[1]]
    m <- list(
      seconds = unlist(lapply(per_round, function(x) x$seconds)),
      peak_mb = mid(field("peak_mb")),
      digest = first$digest,
      digest_kind = first$digest_kind,
      arf_version = first$arf_version
    )
    floor_mb <- mid(field("floor_mb"))
    # A round whose baseline was still contaminated can yield floor > peak.
    # Such a round measured nothing, so drop it rather than let it become the
    # median at a small replicate count.
    paired <- field("peak_mb") - field("floor_mb")
    paired <- paired[!is.na(paired) & paired > 0]
    marginal_mb <- if (length(paired)) stats::median(paired) else NA_real_
    secs <- m$seconds[!is.na(m$seconds)]
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
      peak_delta_mb = round(marginal_mb, 1),
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
