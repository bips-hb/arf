# The benchmark grid as data: which (op, backend, workers, size, variant)
# combinations get measured. See bench/README.md to run it.

.bench_op_knobs <- function(cells, n_evidence = 100L, n_synth = 1L, rowmode = "separate") {
  ev <- cells$op %in% c("forge", "expct")
  cells$n_evidence <- ifelse(ev, n_evidence, NA_integer_)
  cells$n_synth <- ifelse(cells$op == "forge", n_synth, NA_integer_)
  cells$n_folds <- ifelse(cells$op == "lik", 8L, NA_integer_)
  cells$rowmode <- ifelse(ev, rowmode, NA_character_)
  cells
}

.bench_ops <- c("adversarial_rf", "forde", "lik", "forge", "expct")

# quick: sequential, small, minutes. What `make bench` runs.
.bench_cells_quick <- function() {
  cells <- expand.grid(
    op = .bench_ops,
    n = c(1e3, 1e4),
    trees = c(10L, 50L),
    stringsAsFactors = FALSE
  )
  cells$p <- 10L
  cells$backend <- "sequential"
  cells$workers <- NA_integer_
  cells$iters <- 5L
  # Peak memory varies with GC timing, so one round cannot carry a verdict.
  cells$mem_reps <- 3L
  .bench_op_knobs(cells)
}

# full: backends, large data, the worker grid. Cluster only: a slurm cgroup
# caps memory, a workstation cannot.
.bench_cells_full <- function() {
  parallel_cells <- expand.grid(
    op = .bench_ops,
    n = c(1e4, 5e4),
    trees = 200L,
    backend = c("foreach", "psock", "mirai"),
    workers = c(1L, 2L, 4L, 8L, 16L),
    stringsAsFactors = FALSE
  )
  seq_cells <- expand.grid(
    op = .bench_ops,
    n = c(1e4, 5e4),
    trees = 200L,
    backend = "sequential",
    workers = NA_integer_,
    stringsAsFactors = FALSE
  )
  cells <- rbind(parallel_cells, seq_cells)
  cells$p <- 10L
  cells$iters <- 5L
  cells$mem_reps <- 5L
  # "or" delegates parallelism to cforde, the mirai/mori path this suite guards.
  plain <- .bench_op_knobs(cells[!cells$op %in% c("forge", "expct"), ])
  ev <- cells[cells$op %in% c("forge", "expct"), ]
  # Evidence stays at 200: at 1000 the 0.2.5 anchor needs more than a node
  # for expct and the release-trend comparison caps instead of measuring.
  variants <- list(
    list(n_evidence = 100L, n_synth = 1L, rowmode = "separate"),
    list(n_evidence = 200L, n_synth = 100L, rowmode = "separate"),
    list(n_evidence = 100L, n_synth = 1L, rowmode = "or"),
    list(n_evidence = 200L, n_synth = 100L, rowmode = "or")
  )
  ev <- do.call(
    rbind,
    lapply(variants, function(v) {
      .bench_op_knobs(ev, v$n_evidence, v$n_synth, v$rowmode)
    })
  )
  rbind(plain, ev)
}

bench_cells <- function(tier = c("quick", "full")) {
  tier <- match.arg(tier)
  cells <- switch(tier, quick = .bench_cells_quick(), full = .bench_cells_full())
  cells$tier <- tier
  cells[order(cells$op, cells$n, cells$trees, cells$backend), ]
}

# Curated: an anchor is a cost someone chose to pay on every full run.
bench_anchors <- function(path = "bench/anchors.csv") {
  read.csv(path, stringsAsFactors = FALSE)$spec
}

# Shared by run and plan, so the cell set that plans is the cell set that runs.
# Filters are independent: a backend named here keeps that backend at every size.
bench_select_cells <- function(
  cells,
  ops = NA_character_,
  sizes = NA_character_,
  rowmode = NA_character_,
  evidence = NA_character_,
  backends = NA_character_,
  max_cells = NA_integer_
) {
  parts <- function(x) trimws(strsplit(x, ",")[[1]])
  if (!is.na(ops)) {
    cells <- cells[cells$op %in% parts(ops), ]
  }
  if (!is.na(sizes)) {
    cells <- cells[cells$n %in% as.numeric(parts(sizes)), ]
  }
  if (!is.na(backends)) {
    cells <- cells[cells$backend %in% parts(backends), ]
  }
  # The variants sit after the worker grid, out of --max-cells' reach.
  if (!is.na(rowmode)) {
    cells <- cells[!is.na(cells$rowmode) & cells$rowmode %in% parts(rowmode), ]
  }
  if (!is.na(evidence)) {
    cells <- cells[
      !is.na(cells$n_evidence) &
        cells$n_evidence %in% as.integer(parts(evidence)),
    ]
  }
  if (!is.na(max_cells)) {
    cells <- utils::head(cells, max_cells)
  }
  cells
}

# BIPS compute node: 96 cores x 2 SMT, 1.14 TB RealMemory (slurm rejects a
# request above it, so this sits below).
BENCH_NODE_THREADS <- 192L
BENCH_NODE_MEM_MB <- 1100000L

# Cores and memory as matched fractions of a node: requesting them
# independently strands whichever is left over.
BENCH_NODE_FRACTIONS <- c(1 / 16, 1 / 8, 1 / 4, 1 / 2, 1)

# Memory a cell needs in MB, before rounding to a node fraction. Classes from
# measured peaks of cran-0.2.5, the hungriest ref, with headroom: a capped
# measurement is worthless, an over-request only costs queue position.
bench_cell_memory_mb <- function(cells) {
  heavy <- cells$op %in% c("expct", "forge")
  ifelse(
    cells$n <= 1000,
    16000L,
    ifelse(
      cells$n <= 10000,
      ifelse(heavy, 192000L, 48000L),
      ifelse(heavy, .bench_heavy_mb(cells), 192000L)
    )
  )
}

# Per measurement, on top of the timed calls: child start, arf load, fixture
# read, sampler settling. Calibrated against full-tier jobs, within 5%.
BENCH_MEASURE_OVERHEAD_S <- 120
BENCH_WALLTIME_SAFETY <- 1.3
# The medium QoS ceiling; a cell estimated above it cannot complete.
BENCH_WALLTIME_MAX_S <- 86400L
BENCH_WALLTIME_MIN_S <- 1800L

.bench_timing_key <- function(x) {
  paste(x$op, x$backend, x$workers, x$n, x$p, x$trees, x$n_evidence, x$rowmode, sep = "|")
}

# Seconds per timed call from bench/timings.csv. Exact cell, else the slowest
# cell sharing op, size, evidence and row mode: the worker axis transfers
# worst, mirai and foreach scale in opposite directions on expct.
bench_call_seconds <- function(cells, path = "bench/timings.csv") {
  if (!file.exists(path)) {
    return(rep(NA_real_, nrow(cells)))
  }
  t <- utils::read.csv(path, stringsAsFactors = FALSE)
  # Committed, so it can arrive hand-edited or mangled by a merge; the timings
  # are a hint, so a damaged file costs the fallback walltime, not the run.
  need <- c("op", "n", "n_evidence", "rowmode", "call_seconds")
  if (!all(need %in% names(t)) || !is.numeric(t$call_seconds)) {
    warning(
      path,
      " is not usable (needs columns ",
      paste(need, collapse = ", "),
      " with a numeric call_seconds); falling back to the default walltime. ",
      "Regenerate it with arf-bench collect.",
      call. = FALSE
    )
    return(rep(NA_real_, nrow(cells)))
  }
  exact <- t$call_seconds[match(.bench_timing_key(cells), .bench_timing_key(t))]
  cls <- function(x) paste(x$op, x$n, x$n_evidence, x$rowmode, sep = "|")
  worst <- data.table::data.table(k = cls(t), s = t$call_seconds)[,
    list(s = max(s, na.rm = TRUE)),
    by = "k"
  ]
  ifelse(is.na(exact), worst$s[match(cls(cells), worst$k)], exact)
}

# mem_reps rounds, each measuring every ref once: iters timed calls plus a
# floor measurement.
bench_cell_walltime <- function(
  cells,
  n_refs,
  mem_reps = cells$mem_reps,
  iters = cells$iters,
  fallback = NA_real_,
  path = "bench/timings.csv"
) {
  t <- bench_call_seconds(cells, path)
  est <- BENCH_WALLTIME_SAFETY *
    mem_reps *
    n_refs *
    (iters * t + BENCH_MEASURE_OVERHEAD_S)
  est[is.na(est)] <- fallback
  # Half-hour grid: batchtools submits once per distinct resource figure.
  est <- BENCH_WALLTIME_MIN_S * ceiling(est / BENCH_WALLTIME_MIN_S)
  as.integer(pmax(BENCH_WALLTIME_MIN_S, pmin(BENCH_WALLTIME_MAX_S, est)))
}

# expct on cran-0.2.5 at n = 5e4 peaks near 250 GB per 100 evidence rows and
# scales about linearly in the evidence count; row mode moves it far less.
.bench_heavy_mb <- function(cells) {
  much <- !is.na(cells$n_evidence) & cells$n_evidence >= 200
  ifelse(much, 500000L, 275000L)
}

# Two threads per worker plus two spare cores for mirai's dispatcher and the
# sampler: a starved dispatcher makes mirai lose comparisons it should win.
bench_cell_threads <- function(cells) {
  w <- ifelse(is.na(cells$workers), 1L, as.integer(cells$workers))
  2L * (w + 2L)
}

# Per-cell slurm request: the smallest node fraction satisfying both dimensions.
bench_cell_resources <- function(cells) {
  need <- pmax(
    bench_cell_threads(cells) / BENCH_NODE_THREADS,
    bench_cell_memory_mb(cells) / BENCH_NODE_MEM_MB
  )
  frac <- vapply(
    need,
    function(x) {
      ok <- BENCH_NODE_FRACTIONS[BENCH_NODE_FRACTIONS >= x]
      if (length(ok)) min(ok) else 1
    },
    numeric(1)
  )
  data.frame(
    fraction = frac,
    ncpus = as.integer(ceiling(frac * BENCH_NODE_THREADS)),
    memory = as.integer(ceiling(frac * BENCH_NODE_MEM_MB))
  )
}
