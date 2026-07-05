#!/usr/bin/env Rscript
# Benchmark forde() parallel backends on wall-clock time and peak memory.
# See bench/README.md. Not part of the package build (.Rbuildignore ^bench$).
#
# Usage: Rscript bench/bench-backends.R [n] [p] [num_trees] [n_workers]

suppressWarnings(suppressMessages(pkgload::load_all(quiet = TRUE)))

args <- commandArgs(trailingOnly = TRUE)
n         <- if (length(args) >= 1) as.integer(args[1]) else 3000L
p         <- if (length(args) >= 2) as.integer(args[2]) else 30L
num_trees <- if (length(args) >= 3) as.integer(args[3]) else 100L
n_workers <- if (length(args) >= 4) as.integer(args[4]) else 4L

## --- peak memory sampling (Linux PSS across the process subtree) ------------
# Returns peak total PSS in MB while `expr` runs, plus elapsed seconds. On
# non-Linux (no /proc/<pid>/smaps_rollup) the memory value is NA.
proc_descendants <- function(root) {
  pids <- suppressWarnings(as.integer(list.files("/proc")))
  pids <- pids[!is.na(pids)]
  ppid <- integer(0)
  for (pp in pids) {
    st <- tryCatch(readLines(sprintf("/proc/%d/stat", pp), warn = FALSE),
                   error = function(e) NA_character_)
    if (is.na(st[1])) next
    # fields after "pid (comm) state ": ppid is the first
    tail <- sub("^\\d+ \\(.*\\) \\S+ ", "", st)
    ppid[as.character(pp)] <- as.integer(strsplit(tail, " ", fixed = TRUE)[[1]][1])
  }
  out <- root; frontier <- root
  repeat {
    kids <- as.integer(names(ppid)[ppid %in% frontier])
    kids <- setdiff(kids, out)
    if (!length(kids)) break
    out <- c(out, kids); frontier <- kids
  }
  out
}

total_pss_kb <- function(pids) {
  tot <- 0
  for (pp in pids) {
    sr <- tryCatch(readLines(sprintf("/proc/%d/smaps_rollup", pp), warn = FALSE),
                   error = function(e) character(0))
    line <- sr[grepl("^Pss:", sr)]
    if (length(line)) {
      tot <- tot + as.numeric(sub("[^0-9]*([0-9]+).*", "\\1", line[1]))
    }
  }
  tot
}

with_peak_mem <- function(expr, interval = 0.02) {
  linux <- file.exists("/proc/self/smaps_rollup")
  root <- Sys.getpid()
  if (!linux || !requireNamespace("parallel", quietly = TRUE)) {
    t <- system.time(force(expr))["elapsed"]
    return(list(seconds = as.numeric(t), peak_mb = NA_real_))
  }
  stopf <- tempfile()
  sampler <- parallel::mcparallel({
    peak <- 0
    repeat {
      peak <- max(peak, total_pss_kb(proc_descendants(root)))
      if (file.exists(stopf)) break
      Sys.sleep(interval)
    }
    peak
  })
  t <- system.time(force(expr))["elapsed"]
  file.create(stopf)
  peak_kb <- tryCatch(parallel::mccollect(sampler)[[1]], error = function(e) NA_real_)
  unlink(stopf)
  list(seconds = as.numeric(t), peak_mb = peak_kb / 1024)
}

## --- data + fitted forest (shared across backends) --------------------------
set.seed(1)
X <- as.data.frame(matrix(rnorm(n * p), n, p))
X$grp <- factor(sample(letters[1:6], n, replace = TRUE))
message(sprintf("Data: %d x %d  |  num_trees = %d  |  n_workers = %d",
                nrow(X), ncol(X), num_trees, n_workers))
arf <- adversarial_rf(X, num_trees = num_trees, verbose = FALSE, parallel = FALSE)

run <- function(label, setup, teardown = function() NULL) {
  ok <- tryCatch({ setup(); TRUE }, error = function(e) {
    message(sprintf("  [skip] %s: %s", label, conditionMessage(e))); FALSE })
  if (!ok) return(NULL)
  on.exit(teardown(), add = TRUE)
  res <- with_peak_mem(forde(arf, X, parallel = !identical(label, "sequential")))
  data.frame(backend = label, seconds = round(res$seconds, 3),
             peak_mb = round(res$peak_mb, 1))
}

results <- list()

results$seq <- run("sequential", function() options(arf.backend = NULL))

if (requireNamespace("doParallel", quietly = TRUE)) {
  results$fe <- run("foreach",
    setup = function() {
      options(arf.backend = "foreach")
      doParallel::registerDoParallel(cores = n_workers)
    },
    teardown = function() try(doParallel::stopImplicitCluster(), silent = TRUE))
} else message("  [skip] foreach: doParallel not installed")

if (requireNamespace("mirai", quietly = TRUE) &&
    requireNamespace("mori", quietly = TRUE)) {
  results$mi <- run("mirai",
    setup = function() { options(arf.backend = "mirai"); mirai::daemons(n_workers) },
    teardown = function() mirai::daemons(0))
} else message("  [skip] mirai: mirai and/or mori not installed")

cat("\n=== forde() backend benchmark ===\n")
print(do.call(rbind, results), row.names = FALSE)
