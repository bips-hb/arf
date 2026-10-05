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
