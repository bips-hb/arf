#!/usr/bin/env Rscript
# Standalone tests for the benchmark suite. Run from the repo root:
#
#   Rscript bench/tests.R
#
# Not part of the package's testthat suite: bench/ is in .Rbuildignore and
# never ships, so these run on demand rather than in R CMD check.

library(testthat)
source("bench/bench-helpers.R")

test_that("harness helpers are available", {
  expect_true(is.function(bench_measure_cell))
  expect_true(is.function(bench_make_data))
  expect_true(BENCH_METRIC %in% c("cgroup-anon+shmem", "PSS", "RSS"))
})

test_that("generated benchmark data has the documented shape", {
  set.seed(1)
  x <- bench_make_data(50, 10)
  expect_equal(nrow(x), 50)
  expect_true("grp" %in% names(x))
  expect_true(any(vapply(x, is.factor, logical(1))))
})

source("bench/refs.R")

test_that("ref specs parse into kind and value", {
  expect_equal(bench_parse_ref("git:main"), list(kind = "git", value = "main", label = "main"))
  expect_equal(bench_parse_ref("cran:0.2.5"), list(kind = "cran", value = "0.2.5", label = "cran-0.2.5"))
  expect_equal(bench_parse_ref("HEAD")$kind, "git")
  expect_equal(bench_parse_ref("HEAD")$value, "HEAD")
})

test_that("unparseable ref specs are refused", {
  expect_error(bench_parse_ref("svn:trunk"), "unknown ref kind")
  expect_error(bench_parse_ref("cran:"), "needs a version")
  expect_error(bench_parse_ref(""), "empty")
})

# Hermetic on purpose: the live tree holds staged, uncommitted benchmark files
# while this suite runs, so the fixture is a throwaway repo.
.bench_test_repo <- function() {
  repo <- file.path(tempdir(), paste0("benchrepo-", Sys.getpid(), "-", sample.int(1e6, 1)))
  dir.create(repo, recursive = TRUE)
  system2("git", c("-C", repo, "init", "-q"))
  writeLines("x", file.path(repo, "tracked.txt"))
  system2("git", c("-C", repo, "add", "tracked.txt"))
  system2(
    "git",
    c("-C", repo, "-c", "user.email=t@example.org", "-c", "user.name=t", "commit", "-qm", "init"),
    stdout = FALSE
  )
  repo
}

test_that("untracked files do not count as a dirty tree", {
  repo <- .bench_test_repo()
  on.exit(unlink(repo, recursive = TRUE), add = TRUE)
  writeLines("y", file.path(repo, "untracked.txt"))
  expect_equal(bench_tree_dirty(repo), character())
  expect_silent(bench_assert_clean_tree(repo))
})

test_that("a modified tracked file is refused by name", {
  repo <- .bench_test_repo()
  on.exit(unlink(repo, recursive = TRUE), add = TRUE)
  writeLines("changed", file.path(repo, "tracked.txt"))
  expect_equal(bench_tree_dirty(repo), "tracked.txt")
  expect_error(bench_assert_clean_tree(repo), "modified tracked files")
})

test_that("a git ref installs and reports its version", {
  skip_on_cran()
  root <- file.path(tempdir(), paste0("benchlib-", Sys.getpid()))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  ref <- bench_install_ref("git:HEAD", root)
  expect_true(dir.exists(file.path(ref$lib, "arf")))
  expect_match(ref$arf_version, "^[0-9]")
  expect_match(ref$commit, "^[0-9a-f]{7}")
})

test_that("an argument the ref lacks aborts with the ref named", {
  skip_on_cran()
  root <- file.path(tempdir(), paste0("benchlib-args-", Sys.getpid()))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  ref <- bench_install_ref("git:HEAD", root)
  expect_silent(bench_assert_args(ref$lib, list(forde = c("parallel"))))
  expect_error(
    bench_assert_args(ref$lib, list(forde = c("no_such_argument"))),
    "no_such_argument"
  )
})

test_that("a cran version that does not exist fails before any measurement", {
  skip_on_cran()
  skip_if_offline()
  root <- file.path(tempdir(), paste0("benchlib-404-", Sys.getpid()))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  expect_error(bench_install_ref("cran:0.0.0.1", root), "could not download")
})

test_that("every benchmark dependency is present", {
  expect_true(bench_assert_deps())
})

source("bench/cells.R")

test_that("the quick tier is sequential, small, and covers the five ops", {
  cells <- bench_cells("quick")
  expect_setequal(cells$op, c("adversarial_rf", "forde", "lik", "forge", "expct"))
  expect_true(all(cells$backend == "sequential"))
  expect_true(all(cells$n <= 1e4))
  expect_true(all(cells$iters == 5))
})

test_that("no cell uses p < 4, where the 0.2.5 mtry default would differ", {
  expect_true(all(bench_cells("quick")$p >= 4))
  expect_true(all(bench_cells("full")$p >= 4))
})

test_that("the full tier adds backends and a worker grid", {
  cells <- bench_cells("full")
  expect_true(all(c("foreach", "mirai", "psock") %in% cells$backend))
  expect_true(max(cells$workers, na.rm = TRUE) >= 16)
  expect_true(any(cells$n >= 5e4))
})

test_that("sequential cells carry no worker count", {
  cells <- bench_cells("full")
  expect_true(all(is.na(cells$workers[cells$backend == "sequential"])))
})

test_that("anchors are read from the curated file", {
  expect_true("cran:0.2.5" %in% bench_anchors("bench/anchors.csv"))
})

test_that("a fitted model is not fingerprinted", {
  expect_true(is.na(bench_digest(list(forest = 1), "none")))
  expect_equal(unname(BENCH_DIGEST_KIND[["adversarial_rf"]]), "none")
})

test_that("digests are defined for empty and NULL results", {
  expect_type(bench_digest(NULL, "exact"), "character")
  expect_type(bench_digest(data.frame(), "summary"), "character")
  expect_type(bench_digest(numeric(0), "exact"), "character")
})

test_that("digests ignore noise below the rounding tolerance but catch real change", {
  a <- c(1.0000000001, 2)
  b <- c(1, 2)
  expect_equal(bench_digest(a, "exact"), bench_digest(b, "exact"))
  expect_false(identical(bench_digest(c(1.1, 2), "exact"), bench_digest(b, "exact")))
})

test_that("summary digests survive row permutation but not a distribution shift", {
  set.seed(1)
  d <- data.frame(x = rnorm(200), y = rnorm(200))
  expect_equal(bench_digest(d, "summary"), bench_digest(d[sample(200), ], "summary"))
  expect_false(identical(bench_digest(d, "summary"), bench_digest(d * 2, "summary")))
})

test_that("a cell measures an installed ref and returns all three timings", {
  skip_on_cran()
  root <- file.path(tempdir(), paste0("benchcell-", Sys.getpid()))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  ref <- bench_install_ref("git:HEAD", root)
  data_path <- file.path(tempdir(), paste0("benchdata-", Sys.getpid(), ".rds"))
  on.exit(unlink(data_path), add = TRUE)
  # Built in a child against the ref's library: this session resolves arf 0.2.4,
  # which must never be what builds a fixture.
  callr::r(
    function(lib, make_data, path) {
      .libPaths(c(lib, .libPaths()))
      library(arf)
      set.seed(1)
      X <- make_data(300, 10)
      a <- adversarial_rf(X, num_trees = 5, verbose = FALSE, parallel = FALSE)
      saveRDS(list(arf = a, X = X, psi = forde(a, X, parallel = FALSE), evidence = NULL), path)
    },
    args = list(lib = ref$lib, make_data = bench_make_data, path = data_path)
  )

  m <- bench_measure_cell("sequential", data_path, NA_integer_, 1L, ref$lib, iters = 2L, op = "forde")
  expect_length(m$seconds, 2L)
  expect_true(m$peak_mb >= 0)
  expect_equal(m$arf_version, ref$arf_version)
  expect_equal(m$digest_kind, "exact")
  m2 <- bench_measure_cell("sequential", data_path, NA_integer_, 1L, ref$lib, iters = 1L, op = "forde")
  expect_equal(m$digest, m2$digest)
})

source("bench/run-cell.R")

test_that("a cell produces one schema-conforming row per ref", {
  skip_on_cran()
  root <- file.path(tempdir(), paste0("benchrun-", Sys.getpid()))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  refs <- list(bench_install_ref("git:HEAD", root))
  refs[[2]] <- refs[[1]]
  refs[[2]]$label <- "head-copy"

  cell <- bench_cells("quick")[1, ]
  cell$n <- 300
  cell$trees <- 5L
  cell$iters <- 1L
  rows <- bench_run_cell(cell, refs)

  expect_equal(nrow(rows), 2L)
  expect_equal(names(rows), bench_schema())
  expect_equal(rows$digest[1], rows$digest[2])
  expect_true(all(!is.na(rows$metric)))
})

test_that("one ref failing still records the others", {
  skip_on_cran()
  root <- file.path(tempdir(), paste0("benchfail-", Sys.getpid()))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  good <- bench_install_ref("git:HEAD", root)
  bad <- list(
    label = "broken",
    lib = file.path(root, "nonexistent"),
    arf_version = "0.0.0",
    commit = "0000000"
  )

  cell <- bench_cells("quick")[1, ]
  cell$n <- 300
  cell$trees <- 5L
  cell$iters <- 1L
  rows <- bench_run_cell(cell, list(good, bad))

  expect_equal(nrow(rows), 2L)
  expect_false(is.na(rows$time_median[rows$ref == good$label]))
  expect_true(is.na(rows$time_median[rows$ref == "broken"]))
  # Whatever the sampler accumulated before the child died is not a measurement:
  # reported as a number it renders as a large memory "improvement".
  expect_true(is.na(rows$peak_mb[rows$ref == "broken"]))
  expect_true(is.na(rows$peak_delta_mb[rows$ref == "broken"]))
})

source("bench/registry.R")

test_that("a local registry runs a one-cell grid end to end", {
  skip_on_cran()
  skip_if_not_installed("batchtools")
  dir <- file.path(tempdir(), paste0("benchreg-", Sys.getpid()))
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)

  cells <- bench_cells("quick")
  cells <- cells[cells$op == "forde", ][1, ]
  cells$n <- 300
  cells$trees <- 5L
  cells$iters <- 1L

  root <- file.path(tempdir(), paste0("submitlib-", Sys.getpid()))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  reg <- bench_make_registry(dir, "local")
  reg <- bench_submit_cells(reg, cells, list(bench_install_ref("git:HEAD", root)))
  batchtools::waitForJobs(reg = reg)
  rows <- bench_collect(reg)

  expect_equal(nrow(rows), 1L)
  expect_equal(names(rows), bench_schema())
  expect_false(is.na(rows$time_median))
})

source("bench/collate.R")

.bench_fake <- function(...) {
  base <- data.frame(
    ref = NA_character_,
    arf_version = "0.2.5.9000",
    commit = "abc1234",
    op = "lik",
    backend = "sequential",
    workers = NA_integer_,
    n = 1e4,
    p = 10L,
    trees = 50L,
    iters = 5L,
    dt_threads = 1L,
    ranger_threads = 1L,
    metric = "cgroup-anon+shmem",
    peak_mb = NA_real_,
    floor_mb = NA_real_,
    peak_delta_mb = NA_real_,
    mem_reps = 3L,
    mem_limit_mb = NA_real_,
    time_median = NA_real_,
    time_min = NA_real_,
    time_max = NA_real_,
    digest = NA_character_,
    digest_kind = "exact",
    tier = "quick",
    n_evidence = NA_integer_,
    n_synth = NA_integer_,
    n_folds = 8L,
    rowmode = NA_character_,
    host = "h",
    kernel = "k",
    r_version = "4.5",
    job_id = NA_character_,
    timestamp = "2026-10-05T00:00:00",
    stringsAsFactors = FALSE
  )
  rows <- lapply(list(...), function(o) {
    r <- base
    for (nm in names(o)) {
      r[[nm]] <- o[[nm]]
    }
    r
  })
  do.call(rbind, rows)
}

test_that("a memory win above the threshold reads as real and a small time move does not", {
  rows <- .bench_fake(
    list(ref = "main", peak_mb = 412, peak_delta_mb = 200, time_median = 1.84, digest = "a3f1"),
    list(ref = "HEAD", peak_mb = 371, peak_delta_mb = 180, time_median = 1.79, digest = "a3f1")
  )
  d <- bench_deltas(rows, baseline = "main")
  head_row <- d[d$ref == "HEAD", ]
  expect_equal(round(head_row$delta_mem_pct, 1), -10.0)
  expect_equal(head_row$mem_verdict, "real")
  expect_equal(head_row$time_verdict, "inconclusive")
  expect_false(head_row$digest_differs)
})

test_that("the report shows absolute marginals and the floor, not only percentages", {
  rows <- .bench_fake(
    list(
      ref = "main",
      peak_mb = 1450,
      floor_mb = 216,
      peak_delta_mb = 1234,
      mem_reps = 3L,
      time_median = 12.34,
      digest = "a",
      op = "expct"
    ),
    list(
      ref = "cran-0.2.5",
      peak_mb = 7840,
      floor_mb = 216,
      peak_delta_mb = 7623,
      mem_reps = 3L,
      time_median = 12.6,
      digest = "a",
      op = "expct"
    )
  )
  txt <- paste(bench_report(rows, baseline = "main"), collapse = "\n")
  # a +517% delta is unreadable without knowing whether it is MB or GB
  expect_match(txt, "7623.0 vs 1234.0")
  expect_match(txt, "216.0")
  expect_match(txt, "12.600 vs 12.340")
})

test_that("a digest mismatch is flagged in the report", {
  rows <- .bench_fake(
    list(ref = "main", peak_mb = 412, peak_delta_mb = 200, time_median = 1.84, digest = "a3f1"),
    list(ref = "HEAD", peak_mb = 371, peak_delta_mb = 180, time_median = 1.79, digest = "9c2e")
  )
  txt <- paste(bench_report(rows, baseline = "main"), collapse = "\n")
  expect_match(txt, "OUTPUT DIFFERS")
  expect_match(txt, "not apples-to-apples")
})

test_that("a ref identical to the baseline yields no finding and no division by zero", {
  rows <- .bench_fake(
    list(ref = "main", peak_mb = 400, peak_delta_mb = 200, time_median = 2, digest = "a"),
    list(ref = "HEAD", peak_mb = 400, peak_delta_mb = 200, time_median = 2, digest = "a")
  )
  d <- bench_deltas(rows, baseline = "main")
  expect_equal(d$delta_mem_pct[d$ref == "HEAD"], 0)
  expect_equal(d$mem_verdict[d$ref == "HEAD"], "unchanged")
  expect_false(any(is.nan(d$delta_mem_pct)))
})

test_that("a failed baseline measurement does not render as Inf or NaN", {
  rows <- .bench_fake(
    list(ref = "main", peak_mb = 0, peak_delta_mb = 0, time_median = NA_real_, digest = NA_character_),
    list(ref = "HEAD", peak_mb = 371, peak_delta_mb = 180, time_median = 1.79, digest = "a3f1")
  )
  d <- bench_deltas(rows, baseline = "main")
  expect_true(is.na(d$delta_mem_pct[d$ref == "HEAD"]))
  expect_equal(d$mem_verdict[d$ref == "HEAD"], "no baseline")
  txt <- paste(bench_report(rows, baseline = "main"), collapse = "\n")
  expect_false(grepl("Inf|NaN", txt))
})

test_that("history keeps anchor rows only", {
  path <- file.path(tempdir(), paste0("hist-", Sys.getpid(), ".csv"))
  on.exit(unlink(path), add = TRUE)
  rows <- .bench_fake(
    list(ref = "HEAD", peak_mb = 371, peak_delta_mb = 180, time_median = 1.79, digest = "a"),
    list(ref = "cran-0.2.5", peak_mb = 412, peak_delta_mb = 200, time_median = 1.84, digest = "b")
  )
  h1 <- bench_history_append(rows, path)
  expect_equal(nrow(h1), 1L)
  expect_equal(h1$ref, "cran-0.2.5")
  # Re-appending the same configuration replaces it rather than duplicating:
  # this file is committed and plotted.
  h2 <- bench_history_append(rows, path)
  expect_equal(nrow(h2), 1L)
  expect_equal(names(read.csv(path)), bench_schema())
})

test_that("installing the same ref twice reuses the library instead of rebuilding", {
  skip_on_cran()
  root <- file.path(tempdir(), paste0("benchcache-", Sys.getpid()))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  first <- bench_install_ref("git:HEAD", root)
  stamp <- file.path(first$lib, ".bench-stamp")
  expect_true(file.exists(stamp))
  mtime <- file.mtime(stamp)
  second <- bench_install_ref("git:HEAD", root)
  expect_equal(second$commit, first$commit)
  expect_equal(second$arf_version, first$arf_version)
  expect_equal(file.mtime(stamp), mtime)
})

test_that("a zero on either side of a delta is not a finding", {
  expect_true(is.na(.bench_pct(0, 400)))
  expect_true(is.na(.bench_pct(400, 0)))
  expect_true(is.na(.bench_pct(NA_real_, 400)))
  expect_equal(.bench_pct(360, 400), -10)
})

test_that("the memory floor is measured and the marginal excludes it", {
  skip_on_cran()
  root <- file.path(tempdir(), paste0("benchfloor-", Sys.getpid()))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  ref <- bench_install_ref("git:HEAD", root)
  cell <- bench_cells("quick")[1, ]
  cell$op <- "forde"
  cell$n <- 2000
  cell$trees <- 20L
  cell$iters <- 1L
  rows <- bench_run_cell(cell, list(ref))

  # The floor is an R interpreter plus arf plus the fixture: hundreds of MB here,
  # and it dwarfs a small op, so percentages computed on peak_mb are diluted.
  expect_true(rows$floor_mb > 50)
  expect_true(rows$peak_mb >= rows$floor_mb)
  # peak_delta_mb is the median of PAIRED per-round differences, so it is
  # deliberately not the difference of the two reported medians. The contract
  # is that it is a positive quantity strictly smaller than the raw peak.
  expect_gt(rows$peak_delta_mb, 0)
  expect_lt(rows$peak_delta_mb, rows$peak_mb)
  expect_true(rows$peak_delta_mb < rows$peak_mb)
})

test_that("deltas are computed on the marginal, not the floor-inflated total", {
  rows <- .bench_fake(
    list(ref = "main", peak_mb = 412, peak_delta_mb = 200, time_median = 1.84, digest = "a"),
    list(ref = "HEAD", peak_mb = 392, peak_delta_mb = 180, time_median = 1.84, digest = "a")
  )
  d <- bench_deltas(rows, baseline = "main")
  # -20 MB of 200 MB of actual work is -10%, not the -4.9% the totals suggest.
  expect_equal(round(d$delta_mem_pct[d$ref == "HEAD"], 1), -10.0)
  expect_equal(d$mem_verdict[d$ref == "HEAD"], "real")
})

test_that("an unnamed list digests instead of erroring", {
  expect_type(bench_digest(list(1, 2), "exact"), "character")
  expect_type(bench_digest(list(a = 1, 2), "exact"), "character")
})

test_that("provenance is asserted against the library actually loaded", {
  expect_error(.bench_assert_lib("/nonexistent/lib"), "not from the ref library")
  expect_true(.bench_assert_lib(dirname(find.package("arf"))))
})

test_that("ref specs resolving to the same commit collapse to one", {
  skip_on_cran()
  # Branch-independent: on main, "HEAD" and "main" are the same commit and
  # collapse; on a feature branch they differ and both are kept, which is the
  # comparison the suite exists for. Pin the duplicate explicitly instead.
  sha <- system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE)[1]
  resolved <- bench_resolve_refs(c("HEAD", paste0("git:", sha)))
  expect_length(resolved, 1L)
  expect_equal(resolved, "HEAD")
  expect_length(bench_resolve_refs(c("HEAD", "cran:0.2.5")), 2L)
})

test_that("the committed history file does not dirty the tree for the next run", {
  repo <- .bench_test_repo()
  on.exit(unlink(repo, recursive = TRUE), add = TRUE)
  dir.create(file.path(repo, "bench"))
  writeLines("a", file.path(repo, "bench/history.csv"))
  system2("git", c("-C", repo, "add", "bench/history.csv"))
  system2(
    "git",
    c("-C", repo, "-c", "user.email=t@example.org", "-c", "user.name=t", "commit", "-qm", "history"),
    stdout = FALSE
  )
  writeLines("a,b", file.path(repo, "bench/history.csv"))
  # Benchmark OUTPUT, not benchmark input: it cannot make a result irreproducible.
  expect_silent(bench_assert_clean_tree(repo))
  writeLines("changed", file.path(repo, "tracked.txt"))
  expect_error(bench_assert_clean_tree(repo), "tracked.txt")
})

test_that("the full tier covers both evidence row modes and a large variant", {
  cells <- bench_cells("full")
  rowmodes <- cells$rowmode[!is.na(cells$rowmode)]
  expect_true(all(c("separate", "or") %in% rowmodes))
  expect_true(any(cells$n_evidence == 1000, na.rm = TRUE))
  expect_true(any(cells$n_synth == 100, na.rm = TRUE))
})

test_that("errored jobs are surfaced rather than silently dropped", {
  skip_on_cran()
  skip_if_not_installed("batchtools")
  dir <- file.path(tempdir(), paste0("bencherr-", Sys.getpid()))
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  reg <- bench_make_registry(dir, "local")
  batchtools::batchMap(
    function(i) {
      if (i == 2L) {
        stop("boom")
      }
      row <- setNames(as.data.frame(as.list(rep(NA, length(bench_schema())))), bench_schema())
      row$ref <- "main"
      row
    },
    i = 1:2,
    reg = reg
  )
  suppressWarnings(batchtools::submitJobs(reg = reg))
  batchtools::waitForJobs(reg = reg, stop.on.error = FALSE)
  expect_warning(rows <- bench_collect(reg), "1 job")
  expect_equal(nrow(rows), 1L)
})

test_that("the report names the worker count and refuses to mix metrics", {
  rows <- .bench_fake(
    list(
      ref = "main",
      peak_mb = 412,
      peak_delta_mb = 200,
      time_median = 1.84,
      digest = "a",
      backend = "mirai",
      workers = 8L
    ),
    list(
      ref = "HEAD",
      peak_mb = 371,
      peak_delta_mb = 180,
      time_median = 1.79,
      digest = "a",
      backend = "mirai",
      workers = 8L
    )
  )
  txt <- paste(bench_report(rows, baseline = "main"), collapse = "\n")
  expect_match(txt, "\\| 8 \\|")

  mixed <- rows
  mixed$metric[2] <- "PSS"
  expect_warning(bench_report(mixed, baseline = "main"), "metric")
})

test_that("a single memory sample is never reported as a finding", {
  rows <- .bench_fake(
    list(ref = "main", peak_mb = 412, peak_delta_mb = 200, mem_reps = 1L, time_median = 1.84, digest = "a"),
    list(ref = "HEAD", peak_mb = 371, peak_delta_mb = 180, mem_reps = 1L, time_median = 1.84, digest = "a")
  )
  d <- bench_deltas(rows, baseline = "main")
  expect_equal(round(d$delta_mem_pct[d$ref == "HEAD"], 1), -10.0)
  expect_equal(d$mem_verdict[d$ref == "HEAD"], "single sample, not a finding")

  rows$mem_reps <- 3L
  d3 <- bench_deltas(rows, baseline = "main")
  expect_equal(d3$mem_verdict[d3$ref == "HEAD"], "real")
})

test_that("a registry path whose parent does not exist is created", {
  skip_on_cran()
  skip_if_not_installed("batchtools")
  # batchtools asserts the dirname exists rather than creating it, and the CLI
  # nests the registry under bench/registry/<stamp>.
  root <- file.path(tempdir(), paste0("regparent-", Sys.getpid()))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  nested <- file.path(root, "registry", "20260101-000000")
  expect_false(dir.exists(dirname(nested)))
  reg <- bench_make_registry(nested, "local")
  expect_true(dir.exists(nested))
})

test_that("each tier carries its own memory replication budget", {
  expect_equal(unique(bench_cells("quick")$mem_reps), 3L)
  expect_equal(unique(bench_cells("full")$mem_reps), 5L)
})

test_that("a marginal too small relative to its floor carries no memory verdict", {
  thin <- .bench_fake(
    list(
      ref = "main",
      peak_mb = 255,
      floor_mb = 200,
      peak_delta_mb = 55,
      mem_reps = 3L,
      time_median = 1.84,
      digest = "a"
    ),
    list(
      ref = "HEAD",
      peak_mb = 250,
      floor_mb = 200,
      peak_delta_mb = 48,
      mem_reps = 3L,
      time_median = 1.84,
      digest = "a"
    )
  )
  d <- bench_deltas(thin, baseline = "main")
  expect_equal(d$mem_verdict[d$ref == "HEAD"], "cell too small to resolve memory")

  # ratio 0.73 measured a 4.2% same-commit delta, so it carries no verdict either
  borderline <- .bench_fake(
    list(
      ref = "main",
      peak_mb = 370,
      floor_mb = 214,
      peak_delta_mb = 156,
      mem_reps = 3L,
      time_median = 1.84,
      digest = "a"
    ),
    list(
      ref = "HEAD",
      peak_mb = 364,
      floor_mb = 214,
      peak_delta_mb = 150,
      mem_reps = 3L,
      time_median = 1.84,
      digest = "a"
    )
  )
  db <- bench_deltas(borderline, baseline = "main")
  expect_equal(db$mem_verdict[db$ref == "HEAD"], "cell too small to resolve memory")

  fat <- .bench_fake(
    list(
      ref = "main",
      peak_mb = 825,
      floor_mb = 275,
      peak_delta_mb = 550,
      mem_reps = 3L,
      time_median = 1.84,
      digest = "a"
    ),
    list(
      ref = "HEAD",
      peak_mb = 770,
      floor_mb = 275,
      peak_delta_mb = 495,
      mem_reps = 3L,
      time_median = 1.84,
      digest = "a"
    )
  )
  d2 <- bench_deltas(fat, baseline = "main")
  expect_equal(d2$mem_verdict[d2$ref == "HEAD"], "real")
})

test_that("the resolution guard applies to the whole cell, not one row", {
  # One ref thin, the other not: the comparison is still unresolvable, and the
  # cluster's first run reported +10.3% "real" for exactly this shape.
  rows <- .bench_fake(
    list(
      ref = "main",
      peak_mb = 300,
      floor_mb = 200,
      peak_delta_mb = 100,
      mem_reps = 3L,
      time_median = 2,
      digest = NA_character_,
      digest_kind = "none"
    ),
    list(
      ref = "cran-0.2.5",
      peak_mb = 530,
      floor_mb = 200,
      peak_delta_mb = 330,
      mem_reps = 3L,
      time_median = 2,
      digest = NA_character_,
      digest_kind = "none"
    )
  )
  d <- bench_deltas(rows, baseline = "main")
  expect_equal(
    d$mem_verdict[d$ref == "cran-0.2.5"],
    "cell too small to resolve memory"
  )
})

test_that("history replaces a re-run configuration instead of duplicating it", {
  path <- file.path(tempdir(), paste0("hist2-", Sys.getpid(), ".csv"))
  on.exit(unlink(path), add = TRUE)
  mk <- function(ts, mb) {
    .bench_fake(list(
      ref = "cran-0.2.5",
      peak_mb = mb + 200,
      floor_mb = 200,
      peak_delta_mb = mb,
      mem_reps = 3L,
      time_median = 1.5,
      digest = "a",
      timestamp = ts
    ))
  }
  bench_history_append(mk("2026-10-01T00:00:00", 100), path)
  h <- bench_history_append(mk("2026-10-02T00:00:00", 120), path)
  expect_equal(nrow(h), 1L)
  expect_equal(h$peak_delta_mb, 120)

  # a different cell is a different row, not a replacement
  other <- mk("2026-10-02T00:00:00", 300)
  other$n <- 5e4
  h2 <- bench_history_append(other, path)
  expect_equal(nrow(h2), 2L)
})

test_that("history refuses rows that measured nothing", {
  path <- file.path(tempdir(), paste0("hist3-", Sys.getpid(), ".csv"))
  on.exit(unlink(path), add = TRUE)
  dead <- .bench_fake(list(
    ref = "cran-0.2.5",
    peak_mb = NA_real_,
    floor_mb = NA_real_,
    peak_delta_mb = NA_real_,
    time_median = NA_real_,
    digest = NA_character_
  ))
  expect_equal(nrow(bench_history_append(dead, path)), 0L)
})

test_that("a run too large to tabulate is summarised with only real findings", {
  many <- do.call(
    rbind,
    lapply(seq_len(50), function(i) {
      .bench_fake(
        list(
          ref = "main",
          peak_mb = 1000,
          floor_mb = 200,
          peak_delta_mb = 800,
          mem_reps = 3L,
          time_median = 10,
          digest = "a",
          n = i * 1000
        ),
        list(
          ref = "HEAD",
          peak_mb = 1000,
          floor_mb = 200,
          peak_delta_mb = if (i == 1L) 400 else 800,
          mem_reps = 3L,
          time_median = 10,
          digest = "a",
          n = i * 1000
        )
      )
    })
  )
  txt <- paste(bench_report(many, baseline = "main"), collapse = "\n")
  expect_match(txt, "50 comparisons across 50 cells")
  expect_match(txt, "1 comparison\\(s\\) with a verdict")
  # the one real finding is shown, the 49 unchanged cells are not
  expect_equal(lengths(regmatches(txt, gregexpr("\\| `HEAD` \\|", txt))), 1L)
})

test_that("a percentage on a millisecond measurement is not a time finding", {
  # From a real quick-tier run: adversarial_rf at n=1e3 reported +12.7% "real"
  # on a 62 ms call, while every cell above a second read inconclusive.
  fast <- .bench_fake(
    list(
      ref = "main",
      peak_mb = 215,
      floor_mb = 200,
      peak_delta_mb = 14.3,
      mem_reps = 3L,
      time_median = 0.055,
      digest = "a"
    ),
    list(
      ref = "HEAD",
      peak_mb = 215,
      floor_mb = 200,
      peak_delta_mb = 14.5,
      mem_reps = 3L,
      time_median = 0.062,
      digest = "a"
    )
  )
  d <- bench_deltas(fast, baseline = "main")
  expect_equal(d$time_verdict[d$ref == "HEAD"], "too fast to time reliably")

  slow <- fast
  slow$time_median <- c(24.9, 27.5)
  d2 <- bench_deltas(slow, baseline = "main")
  expect_equal(d2$time_verdict[d2$ref == "HEAD"], "real")
})

test_that("a peak at the cgroup cap is not reported as a measurement", {
  # From a real cluster run: six cells sat at 99.75% of a 256 GB cap, identical
  # to the megabyte across a fivefold difference in n.
  rows <- .bench_fake(
    list(
      ref = "main",
      peak_mb = 13500,
      floor_mb = 482,
      peak_delta_mb = 13018,
      mem_reps = 5L,
      mem_limit_mb = 256000,
      time_median = 50,
      digest = "a"
    ),
    list(
      ref = "cran-0.2.5",
      peak_mb = 255354,
      floor_mb = 482,
      peak_delta_mb = 254872,
      mem_reps = 5L,
      mem_limit_mb = 256000,
      time_median = 220,
      digest = "a"
    )
  )
  d <- bench_deltas(rows, baseline = "main")
  expect_equal(
    d$mem_verdict[d$ref == "cran-0.2.5"],
    "hit the memory cap, not a measurement"
  )
  # the capped verdict wins over the resolution and single-sample guards
  expect_equal(
    d$mem_verdict[d$ref == "main"],
    "hit the memory cap, not a measurement"
  )

  roomy <- rows
  roomy$mem_limit_mb <- 900000
  d2 <- bench_deltas(roomy, baseline = "main")
  expect_equal(d2$mem_verdict[d2$ref == "cran-0.2.5"], "real")
})

test_that("an unlimited cgroup yields no cap and no cap verdict", {
  expect_true(is.na(.bench_cgroup_limit_mb(NULL)))
  rows <- .bench_fake(
    list(
      ref = "main",
      peak_mb = 1000,
      floor_mb = 200,
      peak_delta_mb = 800,
      mem_reps = 3L,
      mem_limit_mb = NA_real_,
      time_median = 10,
      digest = "a"
    ),
    list(
      ref = "HEAD",
      peak_mb = 900,
      floor_mb = 200,
      peak_delta_mb = 700,
      mem_reps = 3L,
      mem_limit_mb = NA_real_,
      time_median = 10,
      digest = "a"
    )
  )
  expect_equal(bench_deltas(rows, baseline = "main")$mem_verdict[2], "real")
})

test_that("cores and memory are matched fractions of a node", {
  cells <- bench_cells("full")
  req <- bench_cell_resources(cells)
  # every request is a whole fraction of the node in BOTH dimensions, so
  # neither cores nor memory are stranded behind the other
  expect_true(all(req$ncpus / BENCH_NODE_THREADS - req$memory / BENCH_NODE_MEM_MB < 1e-9))
  expect_true(all(req$ncpus >= bench_cell_threads(cells)))
  expect_true(all(req$memory >= bench_cell_memory_mb(cells)))
  expect_true(all(req$fraction %in% BENCH_NODE_FRACTIONS))
  # a 16-worker cell needs 36 threads, which no 1/16 node (12) can serve
  w16 <- cells$backend != "sequential" & cells$workers == 16L
  expect_true(all(req$ncpus[w16] >= 36L))
})

test_that("memory is requested per cell, from measured peaks", {
  cells <- bench_cells("full")
  mem <- bench_cell_memory_mb(cells)
  # full tier: 48 GB (n=1e4 cheap), 192 GB (n=1e4 heavy and n=5e4 cheap),
  # 500 GB (n=5e4 heavy)
  expect_length(unique(mem), 3L)
  expect_length(unique(bench_cell_memory_mb(bench_cells("quick"))), 3L)
  expect_equal(unique(mem[cells$n == 5e4 & cells$op == "expct"]), 500000L)
  # n=1e4 cheap ops peaked at 17.9 GB; n=5e4 heavy ops were capped at 250 GB
  expect_equal(unique(mem[cells$n == 1e4 & cells$op == "lik"]), 48000L)

  expect_true(all(mem[cells$op %in% c("expct", "forge")] >= mem[match(TRUE, !cells$op %in% c("expct", "forge"))]))
})
