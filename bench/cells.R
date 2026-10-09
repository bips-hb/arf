# The benchmark grid as data. Tiers and sizes come from bench/DESIGN.md.

.bench_op_knobs <- function(cells, n_evidence = 100L, n_synth = 1L, rowmode = "separate") {
  ev <- cells$op %in% c("forge", "expct")
  cells$n_evidence <- ifelse(ev, n_evidence, NA_integer_)
  cells$n_synth <- ifelse(cells$op == "forge", n_synth, NA_integer_)
  cells$n_folds <- ifelse(cells$op == "lik", 8L, NA_integer_)
  cells$rowmode <- ifelse(ev, rowmode, NA_character_)
  cells
}

.bench_ops <- c("adversarial_rf", "forde", "lik", "forge", "expct")

# quick: sequential, small, minutes. What `make bench` runs and what a routine
# PR claim cites.
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
  # Peak memory is a stochastic function of GC scheduling, so one sample cannot
  # support a percentage claim (see bench/DESIGN.md). Three cell runs is the
  # floor for a quick tier that stays cheap enough to run on bertha mid-work.
  cells$mem_reps <- 3L
  .bench_op_knobs(cells)
}

# full: adds the backend comparison, large data and the worker grid. Cluster
# only, because a slurm job's cgroup caps memory and bertha cannot.
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
  # The cluster has the budget the quick tier does not, and these are the runs
  # whose numbers get published.
  cells$mem_reps <- 5L
  # The evidence ops get both row modes and a large variant, as
  # submit-ops.sh already defines: "or" is the branch that delegates
  # parallelism to cforde, which is exactly the mirai/mori code this suite
  # exists to guard. Ops without evidence appear once.
  plain <- .bench_op_knobs(cells[!cells$op %in% c("forge", "expct"), ])
  ev <- cells[cells$op %in% c("forge", "expct"), ]
  variants <- list(
    list(n_evidence = 100L, n_synth = 1L, rowmode = "separate"),
    list(n_evidence = 1000L, n_synth = 100L, rowmode = "separate"),
    list(n_evidence = 100L, n_synth = 1L, rowmode = "or"),
    list(n_evidence = 1000L, n_synth = 100L, rowmode = "or")
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

# Curated: adding an anchor is a decision with a cost someone chose to pay.
bench_anchors <- function(path = "bench/anchors.csv") {
  read.csv(path, stringsAsFactors = FALSE)$spec
}

# Compute-node profile for the BIPS cluster: 192 threads (1 socket x 96 cores
# x 2 SMT) and about 1.1 TB of RAM, so roughly 5.6 GB per thread.
#
# Three numbers describe a node's memory and only one of them is the one that
# matters for a request:
#   MemTotal      1,160,597 MB  what the kernel sees (arf-bench diag on node03)
#   RealMemory    1,120,665 MB  what slurm enforces (scontrol show node node03)
#   sinfo MEMORY    112,066 MB  an order of magnitude low, not to be trusted
# A request above RealMemory is rejected rather than queued, so the figure
# below sits under it with room to spare: the cluster's RealMemory is set
# conservatively on purpose, and the harness should not need it raised.
BENCH_NODE_THREADS <- 192L
BENCH_NODE_MEM_MB <- 1100000L

# Allocate in matched fractions of a node. Requesting cores and memory
# independently strands whichever one is left over: a 36-thread job asking for
# 750 GB blocks 65% of a node's memory behind 19% of its cores, so nothing else
# of the same shape fits and the node is effectively half idle. Taking the same
# fraction of both means a job's footprint is "an eighth of a node" in each
# dimension and the remainder stays usable.
BENCH_NODE_FRACTIONS <- c(1 / 16, 1 / 8, 1 / 4, 1 / 2, 1)

# Memory a cell needs, in MB, before rounding to a node fraction.
#
# Measured peaks, cluster run 2026-10-08 (cran-0.2.5, the hungriest ref):
#   n=1e3, any op                      < 1 GB
#   n=1e4, adversarial_rf/forde/lik   <= 18 GB
#   n=1e4, expct/forge                 39-91 GB
#   n=5e4, expct/forge                 >= 250 GB, true value unknown: these
#                                      cells were clamped at a 256 GB cap
# The heavy classes take headroom rather than a fitted estimate, because a
# clamped measurement is worthless while an over-request only costs queue
# position.
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

# Measured peaks of the hungriest ref in an expct cell at n = 5e4:
#   evidence =  100, separate:  183 GB (cran-0.2.5; HEAD and main need 27)
#   evidence = 1000, separate:  267 GB (no anchor can run it at all)
#   evidence = 1000, "or":      545 GB at 8 workers, and over 550 at 16
# Evidence count and row mode move the peak twentyfold, and the allocator saw
# neither: one figure for every heavy cell both wasted half a node on the light
# corner and lost the heavy one to its own cap. Each figure here is the
# measurement rounded up to the next node fraction.
.bench_heavy_mb <- function(cells) {
  wide <- !is.na(cells$rowmode) & cells$rowmode == "or"
  much <- !is.na(cells$n_evidence) & cells$n_evidence >= 1000
  ifelse(much & wide, 1100000L, ifelse(much | wide, 500000L, 275000L))
}

# Threads a cell needs: two per worker, plus two spare cores for mirai's
# dispatcher and the orchestrator's 5 ms cgroup sampler, since a starved
# dispatcher makes mirai lose comparisons it should win.
bench_cell_threads <- function(cells) {
  w <- ifelse(is.na(cells$workers), 1L, as.integer(cells$workers))
  2L * (w + 2L)
}

# Per-cell slurm request, rounded up to the smallest node fraction that
# satisfies both dimensions.
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
