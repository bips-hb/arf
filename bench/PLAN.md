# arf benchmark suite implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Turn the existing backend-comparison harness into a version-comparing performance and memory suite that can back a claim like "this PR cuts peak memory 10% against main".

**Architecture:** batchtools orchestrates one job per `(op, cell)` with the ref loop inside the job, so every comparison is same-node. Each ref is installed into its own library over a shared dependency library, and each cell measurement runs in a fresh `callr` child whose memory is sampled from the cgroup by the existing helpers. Results land in one schema whose `ref` column drives both the PR report and the committed anchor history.

**Tech Stack:** R, data.table, batchtools, callr, digest, testthat (standalone), GNU make, slurm.

**Spec:** `bench/DESIGN.md`

## Global Constraints

- Everything lives under `bench/`, which is in `.Rbuildignore`. Nothing here may enter the package tarball, `NAMESPACE`, `DESCRIPTION` `Imports`, or the testthat suite in `tests/`.
- Ref set is `HEAD`, `main`, `cran:0.2.5`. Anchors are declared in `bench/anchors.csv`, never derived from all tags.
- Dependencies are installed once into a shared library; a per-ref library holds only `arf`. Child processes set `.libPaths(c(ref_lib, shared_lib))`.
- Install flags must be byte-for-byte identical across refs and must never include `--no-byte-compile`, which would distort timings.
- A tracked-file modification aborts the run. Untracked files are fine.
- The memory metric comes from the existing `bench-helpers.R` sampler unchanged. Every row carries `metric`; rows with different `metric` values are never compared.
- One cell process yields both metrics: 5 iterations inside one `callr` child give `time_median`/`time_min`/`time_max`, and the cgroup peak over that child gives `peak_mb`.
- Thresholds: memory delta above 3% is real, time delta below 10% is inconclusive.
- Digest is exact for `forde` and `lik`, and rounded summary statistics under a fixed seed for `forge` and `expct`. A mismatch is flagged in the report and never aborts.
- New code uses double quotes; run `air format bench/` before finishing a task. Source files end with a newline.
- Lukas commits. Each task ends by leaving the change in the working tree with the suggested commit message recorded, not by running `git commit`.
- Tests are standalone: `bench/tests.R`, run with `Rscript bench/tests.R`, using `testthat::test_that()` and `expect_*`. They must not need a cluster, a network, or more than a few seconds per block except where a task says otherwise.

## Review Focus

- A cell that fails for one ref only, from an out-of-memory kill or an unsupported argument, must still record the other refs' rows rather than losing the whole job. Task 6.
- A `cran:` install with no network or a 404 from both the current and archive URLs must abort before any measurement rather than produce a half-populated report. Task 3.
- Two refs resolving to the same commit, which is the normal case when `HEAD` is `main`, must not yield a divide-by-zero or a 0% delta presented as a finding. Task 8.
- A digest over an op result that is `NULL` or has zero rows, which `forge` can produce under impossible evidence, must return a defined hash rather than error. Task 5.
- A baseline row whose `peak_mb` is `0` or `NA` because measurement failed must not render as `Inf%` or `NaN%` in the report. Task 8.

---

### Task 1: Land the harness on main with correct ignore rules

**Files:**
- Create: `bench/` tracked files, from the `bench-harness` branch
- Modify: `bench/.gitignore`
- Test: `bench/tests.R`

**Interfaces:**
- Consumes: nothing
- Produces: `bench/bench-helpers.R` sourceable from the repo root, exposing `bench_require_backends()`, `bench_ints(env, default)`, `bench_git_commit()`, `bench_make_data(n, p)`, `bench_measure_cell(...)`, and the constants `BENCH_METRIC`, `BENCH_COMMIT`, `BENCH_CGROUP`

- [ ] **Step 1: Take only `bench/` from the stale branch**

The branch predates the #70/#71 housekeeping and would delete tests `main` has, so never merge it. Take the directory alone:

```bash
git checkout bench-harness -- bench/
git status --short bench/
```

- [ ] **Step 2: Write the failing test**

Append to `bench/tests.R` (create it):

```r
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
```

- [ ] **Step 3: Run it**

Run: `Rscript bench/tests.R`
Expected: PASS. If `bench_make_data()` names its factor column something other than `grp`, fix the test to the actual name rather than the function, since `sweep-ops.R` already builds evidence from `levels(X$grp)`.

- [ ] **Step 4: Fix the ignore rules**

The README claims `bench/results/` is gitignored, but `bench/.gitignore` only holds `*.html` and `*.png`. Replace its contents:

```
*.html
*.png
results/
logs/
lib/
*.log
```

- [ ] **Step 5: Verify nothing unexpected is tracked or dirty**

Run: `git status --short bench/`
Expected: the nine harness files staged as additions, plus `tests.R`; no `results/`, `logs/`, or `*.log` entries. The pre-existing untracked experiment scripts (`ux.R`, `psock-direct.R`) stay untracked and visible, which is deliberate: they are Lukas's to keep or delete.

- [ ] **Step 6: Hand off**

Suggested message: `bench: land harness from bench-harness branch, fix ignore rules`

---

### Task 2: Clean-tree precondition and ref parsing

**Files:**
- Create: `bench/refs.R`
- Test: `bench/tests.R`

**Interfaces:**
- Consumes: nothing
- Produces: `bench_parse_ref(spec)` returning `list(kind = "git"|"cran", value = character(1), label = character(1))`; `bench_tree_dirty()` returning `character()` of modified tracked paths; `bench_assert_clean_tree()` returning invisibly or stopping

- [ ] **Step 1: Write the failing test**

Append to `bench/tests.R`:

```r
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

test_that("untracked files do not count as a dirty tree", {
  f <- file.path("bench", paste0("untracked-", Sys.getpid(), ".tmp"))
  on.exit(unlink(f), add = TRUE)
  writeLines("x", f)
  expect_false(f %in% bench_tree_dirty())
  expect_silent(bench_assert_clean_tree())
})
```

- [ ] **Step 2: Run it to verify it fails**

Run: `Rscript bench/tests.R`
Expected: FAIL, cannot open `bench/refs.R`.

- [ ] **Step 3: Implement**

Create `bench/refs.R`:

```r
# Ref resolution for the version-comparing benchmark: see bench/DESIGN.md.

# "git:<ref>", "cran:<version>", or a bare git ref such as "HEAD" or "main".
bench_parse_ref <- function(spec) {
  if (!nzchar(spec)) {
    stop("empty ref spec")
  }
  if (!grepl(":", spec, fixed = TRUE)) {
    return(list(kind = "git", value = spec, label = spec))
  }
  kind <- sub(":.*$", "", spec)
  value <- sub("^[^:]*:", "", spec)
  if (!kind %in% c("git", "cran")) {
    stop("unknown ref kind: ", kind)
  }
  if (!nzchar(value)) {
    stop(kind, " ref needs a version or revision")
  }
  list(
    kind = kind,
    value = value,
    label = if (kind == "cran") paste0("cran-", value) else value
  )
}

# Modified TRACKED paths only. Untracked files are expected here: results/,
# logs/ and local experiment scripts all live under bench/.
bench_tree_dirty <- function() {
  out <- system2("git", c("status", "--porcelain", "--untracked-files=no"), stdout = TRUE)
  if (!length(out)) {
    return(character())
  }
  trimws(substring(out, 4))
}

# A number from a tree nobody can reconstruct cannot be cited, so refuse it
# rather than label it (bench/DESIGN.md, "A dirty tree refuses to run").
bench_assert_clean_tree <- function() {
  dirty <- bench_tree_dirty()
  if (length(dirty)) {
    stop(
      "working tree has modified tracked files, so results would not be reproducible:\n  ",
      paste(dirty, collapse = "\n  "),
      "\nCommit them, including to a throwaway branch, then rerun.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}
```

- [ ] **Step 4: Run it to verify it passes**

Run: `Rscript bench/tests.R`
Expected: PASS. Note the clean-tree test only passes on a clean tree; if it fails because of your own edits, that is the function working.

- [ ] **Step 5: Format and hand off**

Run: `air format bench/`
Suggested message: `bench: add ref parsing and clean-tree precondition`

---

### Task 3: Install refs into per-ref libraries

**Files:**
- Modify: `bench/refs.R`
- Test: `bench/tests.R`

**Interfaces:**
- Consumes: `bench_parse_ref()` from Task 2
- Produces: `bench_assert_deps()` returning invisibly or stopping; `bench_install_ref(spec, root)` returning `list(label, lib, arf_version, commit)`; `bench_assert_args(lib, calls)` returning invisibly or stopping

- [ ] **Step 1: Write the failing test**

Append to `bench/tests.R`:

```r
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
```

- [ ] **Step 2: Run it to verify it fails**

Run: `Rscript bench/tests.R`
Expected: FAIL, `bench_install_ref` not found.

- [ ] **Step 3: Implement**

Append to `bench/refs.R`:

```r
# The shared dependency layer is simply the session's own libpaths: every dep
# already resolves there, and installing them per ref would let their versions
# diverge between refs and confound the comparison with somebody else's
# performance change. A per-ref library therefore holds only arf, and children
# prepend it with `.libPaths(c(ref_lib, .libPaths()))`.
bench_assert_deps <- function() {
  deps <- c("data.table", "ranger", "stringr", "truncnorm", "foreach",
            "doParallel", "mirai", "mori", "digest", "callr", "batchtools")
  missing <- deps[!vapply(deps, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("missing benchmark dependencies: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

# Identical flags for every ref. Never add --no-byte-compile: it would change
# the thing being timed.
.bench_install_flags <- c("--no-docs", "--no-help", "--no-multiarch")

.bench_install_tarball <- function(tarball, lib) {
  args <- c("CMD", "INSTALL", .bench_install_flags, "-l", shQuote(lib), shQuote(tarball))
  status <- system2(file.path(R.home("bin"), "R"), args, stdout = FALSE, stderr = FALSE)
  if (status != 0L) {
    stop("R CMD INSTALL failed for ", basename(tarball), call. = FALSE)
  }
  invisible(TRUE)
}

.bench_install_git <- function(rev, lib) {
  wt <- file.path(tempdir(), paste0("arf-wt-", substr(rev, 1, 12), "-", Sys.getpid()))
  unlink(wt, recursive = TRUE)
  status <- system2("git", c("worktree", "add", "--detach", shQuote(wt), shQuote(rev)),
                    stdout = FALSE, stderr = FALSE)
  if (status != 0L) {
    stop("git worktree add failed for ", rev, call. = FALSE)
  }
  on.exit({
    system2("git", c("worktree", "remove", "--force", shQuote(wt)), stdout = FALSE, stderr = FALSE)
  }, add = TRUE)
  commit <- system2("git", c("-C", shQuote(wt), "rev-parse", "--short", "HEAD"), stdout = TRUE)
  .bench_install_tarball(wt, lib)
  commit
}

.bench_install_cran <- function(version, lib) {
  # The tarball CRAN ships is the faithful anchor, and a current release lives
  # at a different URL than an archived one, so try both.
  urls <- sprintf(
    c("https://cran.r-project.org/src/contrib/arf_%s.tar.gz",
      "https://cran.r-project.org/src/contrib/Archive/arf/arf_%s.tar.gz"),
    version
  )
  dest <- file.path(tempdir(), sprintf("arf_%s.tar.gz", version))
  ok <- FALSE
  for (u in urls) {
    ok <- isTRUE(tryCatch(
      utils::download.file(u, dest, quiet = TRUE, mode = "wb") == 0L,
      error = function(e) FALSE,
      warning = function(w) FALSE
    ))
    if (ok) break
  }
  if (!ok) {
    stop("could not download arf ", version, " from CRAN or its archive", call. = FALSE)
  }
  .bench_install_tarball(dest, lib)
  paste0("cran-", version)
}

bench_install_ref <- function(spec, root) {
  ref <- bench_parse_ref(spec)
  lib <- file.path(root, ref$label)
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  commit <- if (ref$kind == "git") {
    .bench_install_git(ref$value, lib)
  } else {
    .bench_install_cran(ref$value, lib)
  }
  list(
    label = ref$label,
    lib = lib,
    arf_version = as.character(utils::packageVersion("arf", lib.loc = lib)),
    commit = commit
  )
}

# R drops or partially matches an argument a ref does not have, so the call
# would succeed against that version's defaults and produce a plausible number
# for a different computation. Check up front instead.
bench_assert_args <- function(lib, calls) {
  # In a clean child, because this session may already have arf attached from a
  # different library and loadNamespace() would hand back that cached copy.
  have <- callr::r(
    function(lib, fns) {
      .libPaths(c(lib, .libPaths()))
      loadNamespace("arf", lib.loc = lib)
      stats::setNames(lapply(fns, function(f) names(formals(getExportedValue("arf", f)))), fns)
    },
    args = list(lib = lib, fns = names(calls))
  )
  for (fn in names(calls)) {
    missing <- setdiff(calls[[fn]], have[[fn]])
    if (length(missing)) {
      stop("arf in ", lib, " has no argument(s) ", paste(missing, collapse = ", "),
           " for ", fn, "()", call. = FALSE)
    }
  }
  invisible(TRUE)
}
```

- [ ] **Step 4: Run it to verify it passes**

Run: `Rscript bench/tests.R`
Expected: PASS. The git install takes a minute or two because it compiles nothing but does run `R CMD INSTALL`. If `loadNamespace()` in `bench_assert_args()` collides with an already-attached `arf`, run the check in a `callr::r()` child instead, passing `lib`.

- [ ] **Step 5: Format and hand off**

Run: `air format bench/`
Suggested message: `bench: install refs into per-ref libraries over shared deps`

---

### Task 4: The grid as data

**Files:**
- Create: `bench/cells.R`
- Create: `bench/anchors.csv`
- Test: `bench/tests.R`

**Interfaces:**
- Consumes: nothing
- Produces: `bench_cells(tier)` returning a `data.frame` with columns `op, n, p, trees, backend, workers, iters, n_evidence, n_synth, n_folds, rowmode`; `bench_anchors(path)` returning a character vector of ref specs

- [ ] **Step 1: Write the failing test**

Append to `bench/tests.R`:

```r
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
```

- [ ] **Step 2: Run it to verify it fails**

Run: `Rscript bench/tests.R`
Expected: FAIL, cannot open `bench/cells.R`.

- [ ] **Step 3: Implement**

Create `bench/anchors.csv`:

```csv
spec,note
cran:0.2.5,last CRAN release before the mirai/mori work
```

Create `bench/cells.R`:

```r
# The benchmark grid as data. Tiers and sizes come from bench/DESIGN.md.

.bench_op_knobs <- function(cells) {
  cells$n_evidence <- ifelse(cells$op %in% c("forge", "expct"), 100L, NA_integer_)
  cells$n_synth <- ifelse(cells$op == "forge", 1L, NA_integer_)
  cells$n_folds <- ifelse(cells$op == "lik", 8L, NA_integer_)
  cells$rowmode <- ifelse(cells$op %in% c("forge", "expct"), "separate", NA_character_)
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
  .bench_op_knobs(cells)
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
```

- [ ] **Step 4: Run it to verify it passes**

Run: `Rscript bench/tests.R`
Expected: PASS.

- [ ] **Step 5: Format and hand off**

Run: `air format bench/`
Suggested message: `bench: define the grid and anchor set as data`

---

### Task 5: Measure an installed ref and digest its result

**Files:**
- Modify: `bench/bench-helpers.R`
- Test: `bench/tests.R`

**Interfaces:**
- Consumes: `bench_install_ref()` from Task 3
- Produces: `bench_digest(x, kind)` returning `character(1)`; `.bench_cell_fn(lib, data_path, backend, n_workers, dt_threads, ranger_threads, iters, op, op_args)` returning `list(seconds = numeric(iters), digest = character(1), digest_kind = character(1), arf_version = character(1))`; `bench_measure_cell(...)` with a `lib` argument in place of `pkgdir`, returning `list(seconds, peak_mb, digest, digest_kind, arf_version)`

- [ ] **Step 1: Write the failing test**

Append to `bench/tests.R`:

```r
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
  set.seed(1)
  X <- bench_make_data(300, 10)
  data_path <- file.path(tempdir(), paste0("benchdata-", Sys.getpid(), ".rds"))
  on.exit(unlink(data_path), add = TRUE)
  arf <- adversarial_rf(X, num_trees = 5, verbose = FALSE, parallel = FALSE)
  psi <- forde(arf, X, parallel = FALSE)
  saveRDS(list(arf = arf, X = X, psi = psi, evidence = NULL), data_path)

  m <- bench_measure_cell("sequential", data_path, NA_integer_, 1L, ref$lib,
                          iters = 2L, op = "forde")
  expect_length(m$seconds, 2L)
  expect_true(m$peak_mb >= 0)
  expect_equal(m$arf_version, ref$arf_version)
  expect_equal(m$digest_kind, "exact")
  m2 <- bench_measure_cell("sequential", data_path, NA_integer_, 1L, ref$lib,
                           iters = 1L, op = "forde")
  expect_equal(m$digest, m2$digest)
})
```

- [ ] **Step 2: Run it to verify it fails**

Run: `Rscript bench/tests.R`
Expected: FAIL, `bench_digest` not found, and `bench_measure_cell()` rejects `lib`.

- [ ] **Step 3: Implement the digest**

Append to `bench/bench-helpers.R`:

```r
## ---- correctness digest -----------------------------------------------------
# A performance number from a ref that computes something different is worse
# than no number, so every cell fingerprints its own result. "exact" hashes the
# rounded values (forde parameters, lik values). "summary" hashes rounded column
# means and standard deviations, for ops whose output is a random sample
# (forge, expct): an exact hash there also moves when a refactor merely
# consumes RNG draws in a different order.
# adversarial_rf is deliberately "none": its result is a fitted ranger forest,
# where an exact hash would flag benign RNG-consumption differences and the
# summary path has no as.data.frame() method to stand on. The fit's correctness
# belongs to the test suite, not to a benchmark fingerprint.
BENCH_DIGEST_KIND <- c(forde = "exact", lik = "exact", adversarial_rf = "none",
                       forge = "summary", expct = "summary")

.bench_round_num <- function(x, digits) {
  if (is.numeric(x)) round(x, digits) else x
}

bench_digest <- function(x, kind = c("exact", "summary", "none")) {
  kind <- match.arg(kind)
  if (kind == "none") {
    return(NA_character_)
  }
  payload <- if (kind == "exact") {
    .bench_digest_exact(x)
  } else {
    .bench_digest_summary(x)
  }
  digest::digest(payload, algo = "xxhash64")
}

.bench_digest_exact <- function(x) {
  if (is.null(x)) return("NULL")
  if (is.list(x) && !is.data.frame(x)) {
    return(lapply(x[order(names(x))], .bench_digest_exact))
  }
  if (is.data.frame(x)) {
    return(lapply(as.list(x)[order(names(x))], .bench_round_num, digits = 8))
  }
  .bench_round_num(x, 8)
}

# Order-free on purpose: row order is not part of the distribution, and sample
# order is not stable across refactors.
.bench_digest_summary <- function(x) {
  if (is.null(x)) return("NULL")
  x <- as.data.frame(x)
  if (!nrow(x) || !ncol(x)) return("empty")
  cols <- as.list(x)[order(names(x))]
  lapply(cols, function(col) {
    if (is.numeric(col)) {
      round(c(mean(col), stats::sd(col)), 4)
    } else {
      tab <- table(as.character(col))
      round(as.numeric(tab[order(names(tab))]) / length(col), 4)
    }
  })
}
```

- [ ] **Step 4: Switch the child from a source dir to an installed library**

In `bench/bench-helpers.R`, change the first line of `.bench_cell_fn()` body and its signature. Replace:

```r
.bench_cell_fn <- function(pkgdir, data_path, backend, n_workers, dt_threads,
                           ranger_threads, iters, op = "forde",
                           op_args = list()) {
  suppressWarnings(suppressMessages(pkgload::load_all(pkgdir, quiet = TRUE)))
```

with:

```r
.bench_cell_fn <- function(lib, data_path, backend, n_workers, dt_threads,
                           ranger_threads, iters, op = "forde",
                           op_args = list()) {
  # Prepend, never replace: the ref's library holds only arf, and the session's
  # own libpaths are the shared dependency layer, so a comparison cannot be
  # confounded by a different data.table version.
  .libPaths(c(lib, .libPaths()))
  suppressWarnings(suppressMessages(library(arf)))
```

In the same function, replace the two `pkgload::load_all(pkgdir, ...)` calls used for cluster and daemon setup with the installed-library equivalent:

```r
    parallel::clusterCall(cl, function(l, t) {
      .libPaths(c(l, .libPaths()))
      suppressMessages(library(arf))
      data.table::setDTthreads(t)
    }, lib, dt_threads)
```

```r
    if (op %in% c("forge", "expct", "lik")) {
      mirai::everywhere({
        .libPaths(c(lib, .libPaths()))
        suppressMessages(library(arf))
      }, lib = lib)
    }
```

Then replace the timing tail of the function:

```r
  secs <- vapply(seq_len(iters),
                 function(i) system.time(run())[["elapsed"]], numeric(1))
  secs
}
```

with a version that also fingerprints the result:

```r
  kind <- unname(BENCH_DIGEST_KIND[[op]])
  # Digest the first iteration only: one hash per cell is enough, and a fixed
  # seed makes the stochastic ops comparable across refs.
  set.seed(1)
  first <- NULL
  secs <- vapply(seq_len(iters), function(i) {
    t <- system.time(res <- run())[["elapsed"]]
    if (i == 1L) first <<- res
    t
  }, numeric(1))
  list(
    seconds = secs,
    digest = bench_digest(first, kind),
    digest_kind = kind,
    arf_version = as.character(utils::packageVersion("arf"))
  )
}
```

Note `BENCH_DIGEST_KIND` and `bench_digest()` must be passed into the child, since `callr::r_bg()` only ships the function and its arguments. Add them to the `args` list in Step 5 rather than relying on the parent's globals.

- [ ] **Step 5: Update the parent-side measurement**

In `bench_measure_cell()`, rename `pkgdir` to `lib`, pass the digest helpers through, and return the full timing vector. Replace the `callr::r_bg()` call and the return value:

```r
  proc <- callr::r_bg(
    function(lib, data_path, backend, n_workers, dt_threads, ranger_threads,
             iters, op, op_args, helpers) {
      for (nm in names(helpers)) assign(nm, helpers[[nm]], envir = globalenv())
      .bench_cell_fn(lib, data_path, backend, n_workers, dt_threads,
                     ranger_threads, iters, op, op_args)
    },
    args = list(lib, data_path, backend, n_workers, dt_threads, ranger_threads,
                iters, op, op_args,
                helpers = list(.bench_cell_fn = .bench_cell_fn,
                               bench_digest = bench_digest,
                               .bench_digest_exact = .bench_digest_exact,
                               .bench_digest_summary = .bench_digest_summary,
                               .bench_round_num = .bench_round_num,
                               BENCH_DIGEST_KIND = BENCH_DIGEST_KIND,
                               `%||%` = `%||%`))
  )
```

```r
  res <- tryCatch(proc$get_result(),
                  error = function(e) list(seconds = NA_real_, digest = NA_character_,
                                           digest_kind = NA_character_,
                                           arf_version = NA_character_))
  list(seconds = res$seconds, peak_mb = peak_kb / 1024, digest = res$digest,
       digest_kind = res$digest_kind, arf_version = res$arf_version)
```

- [ ] **Step 6: Update the existing callers**

`bench/sweep.R`, `bench/sweep-ops.R` and `bench/mem-backends.R` pass `pkgdir` and read `m$seconds` as a scalar. In each, replace the `pkgdir <- normalizePath(".")` line with an installed library of the working tree and wrap the scalar use:

```r
bench_root <- file.path(tempdir(), "benchlib")
lib <- bench_install_ref("git:HEAD", bench_root)$lib
```

and change `m$seconds` to `stats::median(m$seconds)` at each use site.

- [ ] **Step 7: Run it to verify it passes**

Run: `Rscript bench/tests.R`
Expected: PASS.

- [ ] **Step 8: Format and hand off**

Run: `air format bench/`
Suggested message: `bench: measure installed refs and fingerprint results`

---

### Task 6: Run one cell against all refs

**Files:**
- Create: `bench/run-cell.R`
- Test: `bench/tests.R`

**Interfaces:**
- Consumes: `bench_cells()` (Task 4), `bench_install_ref()` (Task 3), `bench_measure_cell()` (Task 5)
- Produces: `bench_schema()` returning the ordered column names; `bench_run_cell(cell, refs, data_path)` returning one `data.frame` row per ref with exactly those columns

- [ ] **Step 1: Write the failing test**

Append to `bench/tests.R`:

```r
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
  bad <- list(label = "broken", lib = file.path(root, "nonexistent"),
              arf_version = "0.0.0", commit = "0000000")

  cell <- bench_cells("quick")[1, ]
  cell$n <- 300
  cell$trees <- 5L
  cell$iters <- 1L
  rows <- bench_run_cell(cell, list(good, bad))

  expect_equal(nrow(rows), 2L)
  expect_false(is.na(rows$time_median[rows$ref == good$label]))
  expect_true(is.na(rows$time_median[rows$ref == "broken"]))
})
```

- [ ] **Step 2: Run it to verify it fails**

Run: `Rscript bench/tests.R`
Expected: FAIL, cannot open `bench/run-cell.R`.

- [ ] **Step 3: Implement**

Create `bench/run-cell.R`:

```r
# One cell measured against every ref, in one process, on one node.
# The ref loop lives HERE rather than in the batchtools grid on purpose: if
# `ref` were a job dimension, refs would scatter across nodes and reintroduce
# the drift that same-node A/B exists to remove (bench/DESIGN.md).

bench_schema <- function() {
  c("ref", "arf_version", "commit", "op", "backend", "workers", "n", "p",
    "trees", "iters", "dt_threads", "ranger_threads", "metric", "peak_mb",
    "time_median", "time_min", "time_max", "digest", "digest_kind", "tier",
    "n_evidence", "n_synth", "n_folds", "rowmode", "host", "kernel",
    "r_version", "job_id", "timestamp")
}

# Build the fixture once per cell and let every ref read the same file, so a
# comparison never differs by its input data.
.bench_cell_data <- function(cell) {
  set.seed(1)
  X <- bench_make_data(cell$n, cell$p)
  arf <- adversarial_rf(X, num_trees = cell$trees, verbose = FALSE, parallel = FALSE)
  psi <- forde(arf, X, parallel = FALSE)
  evidence <- if (!is.na(cell$n_evidence)) {
    data.frame(grp = sample(levels(X$grp), cell$n_evidence, replace = TRUE))
  } else {
    NULL
  }
  path <- tempfile(fileext = ".rds")
  saveRDS(list(arf = arf, X = X, psi = psi, evidence = evidence), path)
  path
}

.bench_op_args <- function(cell) {
  switch(cell$op,
    forge = list(n_synth = cell$n_synth, stepsize = 0L, rowmode = cell$rowmode),
    expct = list(stepsize = 0L, rowmode = cell$rowmode),
    lik = list(batch = ceiling(cell$n / cell$n_folds)),
    adversarial_rf = list(trees = cell$trees),
    list())
}

bench_run_cell <- function(cell, refs, data_path = NULL,
                           dt_threads = 1L, ranger_threads = 1L,
                           job_id = NA_character_) {
  stopifnot(nrow(cell) == 1L)
  if (is.null(data_path)) {
    data_path <- .bench_cell_data(cell)
    on.exit(unlink(data_path), add = TRUE)
  }
  uname <- tryCatch(system2("uname", "-r", stdout = TRUE), error = function(e) NA_character_)
  rows <- lapply(refs, function(ref) {
    # A ref that dies, from an OOM kill or a missing library, must not cost the
    # other refs their measurements.
    m <- tryCatch(
      bench_measure_cell(cell$backend, data_path,
                         if (is.na(cell$workers)) NA_integer_ else cell$workers,
                         dt_threads, ref$lib, ranger_threads = ranger_threads,
                         iters = cell$iters, op = cell$op,
                         op_args = .bench_op_args(cell)),
      error = function(e) {
        message("  ref ", ref$label, " failed: ", conditionMessage(e))
        list(seconds = NA_real_, peak_mb = NA_real_, digest = NA_character_,
             digest_kind = NA_character_, arf_version = ref$arf_version)
      }
    )
    secs <- m$seconds[!is.na(m$seconds)]
    data.frame(
      ref = ref$label, arf_version = ref$arf_version, commit = ref$commit,
      op = cell$op, backend = cell$backend, workers = cell$workers,
      n = cell$n, p = cell$p, trees = cell$trees, iters = cell$iters,
      dt_threads = dt_threads, ranger_threads = ranger_threads,
      metric = BENCH_METRIC, peak_mb = round(m$peak_mb, 1),
      time_median = if (length(secs)) round(stats::median(secs), 3) else NA_real_,
      time_min = if (length(secs)) round(min(secs), 3) else NA_real_,
      time_max = if (length(secs)) round(max(secs), 3) else NA_real_,
      digest = m$digest, digest_kind = m$digest_kind, tier = cell$tier,
      n_evidence = cell$n_evidence, n_synth = cell$n_synth,
      n_folds = cell$n_folds, rowmode = cell$rowmode,
      host = Sys.info()[["nodename"]], kernel = uname[1],
      r_version = paste0(R.version$major, ".", R.version$minor),
      job_id = job_id, timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[, bench_schema()]
}
```

- [ ] **Step 4: Run it to verify it passes**

Run: `Rscript bench/tests.R`
Expected: PASS.

- [ ] **Step 5: Format and hand off**

Run: `air format bench/`
Suggested message: `bench: run one cell against all refs, one row per ref`

---

### Task 7: batchtools orchestration

**Files:**
- Create: `bench/registry.R`
- Create: `bench/run.R`
- Test: `bench/tests.R`

**Interfaces:**
- Consumes: `bench_cells()` (Task 4), `bench_install_ref()` (Task 3), `bench_run_cell()` (Task 6), `bench_assert_clean_tree()` (Task 2)
- Produces: `bench_make_registry(dir, cluster)` returning a batchtools registry; `bench_submit(tier, refs, cluster, dir)` returning the registry after submitting; `bench_collect(reg)` returning a `data.frame` in `bench_schema()` order

- [ ] **Step 1: Write the failing test**

Append to `bench/tests.R`:

```r
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

  reg <- bench_make_registry(dir, "local")
  reg <- bench_submit_cells(reg, cells, c("git:HEAD"))
  batchtools::waitForJobs(reg = reg)
  rows <- bench_collect(reg)

  expect_equal(nrow(rows), 1L)
  expect_equal(names(rows), bench_schema())
  expect_false(is.na(rows$time_median))
})
```

- [ ] **Step 2: Run it to verify it fails**

Run: `Rscript bench/tests.R`
Expected: FAIL, cannot open `bench/registry.R`.

- [ ] **Step 3: Implement the registry**

Create `bench/registry.R`:

```r
# batchtools orchestrates; the cgroup sampler in bench-helpers.R measures.
# Slurm's own accounting (MaxRSS) is unusable for this comparison, which is why
# the metric never comes from the scheduler (bench/DESIGN.md).

bench_make_registry <- function(dir, cluster = c("local", "slurm"),
                                slurm_template = Sys.getenv("ARF_BENCH_SLURM_TMPL", "")) {
  cluster <- match.arg(cluster)
  reg <- batchtools::makeRegistry(
    file.dir = dir,
    source = c("bench/bench-helpers.R", "bench/refs.R", "bench/cells.R",
               "bench/run-cell.R"),
    packages = character()
  )
  reg$cluster.functions <- if (cluster == "slurm") {
    if (!nzchar(slurm_template)) {
      stop("set ARF_BENCH_SLURM_TMPL to the BIPS cluster slurm template", call. = FALSE)
    }
    batchtools::makeClusterFunctionsSlurm(template = slurm_template)
  } else {
    batchtools::makeClusterFunctionsInteractive()
  }
  reg
}

# One job per (op, cell) with every ref inside it. Refs are installed INSIDE
# the job so each comparison uses libraries built on the node that measures it.
.bench_job <- function(cell, ref_specs, lib_root) {
  refs <- lapply(ref_specs, bench_install_ref, root = lib_root)
  calls <- list(forde = "parallel", lik = c("arf", "batch", "parallel"),
                forge = c("n_synth", "evidence", "stepsize", "evidence_row_mode"),
                expct = c("evidence", "stepsize", "evidence_row_mode"),
                adversarial_rf = c("num_trees", "parallel"))
  for (ref in refs) {
    bench_assert_args(ref$lib, calls[cell$op])
  }
  bench_run_cell(cell, refs, job_id = Sys.getenv("SLURM_JOB_ID", NA_character_))
}

bench_submit_cells <- function(reg, cells, ref_specs,
                               lib_root = file.path(tempdir(), "arf-bench-lib"),
                               resources = list()) {
  batchtools::batchMap(
    function(i, ref_specs, lib_root) .bench_job(CELLS[i, ], ref_specs, lib_root),
    i = seq_len(nrow(cells)),
    more.args = list(ref_specs = ref_specs, lib_root = lib_root),
    reg = reg
  )
  # CELLS travels to the workers through the registry's export, keeping the
  # job parameter a single index so the registry stays readable.
  batchtools::batchExport(list(CELLS = cells), reg = reg)
  batchtools::submitJobs(resources = resources, reg = reg)
  reg
}

bench_collect <- function(reg) {
  done <- batchtools::findDone(reg = reg)
  if (!nrow(done)) {
    stop("no jobs finished; see batchtools::getErrorMessages()", call. = FALSE)
  }
  rows <- batchtools::reduceResultsList(ids = done, reg = reg)
  out <- do.call(rbind, rows)
  out[, bench_schema()]
}
```

- [ ] **Step 4: Implement the entry point**

Create `bench/run.R`:

```r
#!/usr/bin/env Rscript
# Entry point for the version-comparing benchmark. See bench/DESIGN.md.
#
#   Rscript bench/run.R                       # quick tier, local
#   ARF_BENCH_TIER=full ARF_BENCH_CLUSTER=slurm Rscript bench/run.R
#
# Refs default to HEAD, main and the curated anchors.

source("bench/bench-helpers.R")
source("bench/refs.R")
source("bench/cells.R")
source("bench/run-cell.R")
source("bench/registry.R")

bench_assert_clean_tree()

tier <- Sys.getenv("ARF_BENCH_TIER", "quick")
cluster <- Sys.getenv("ARF_BENCH_CLUSTER", "local")
refs <- c("HEAD", "main", bench_anchors())
refs <- unique(refs[nzchar(refs)])

cells <- bench_cells(tier)
message(sprintf("arf bench | tier %s | %d cells | refs %s | metric %s",
                tier, nrow(cells), paste(refs, collapse = ", "), BENCH_METRIC))

dir <- file.path("bench", "registry", format(Sys.time(), "%Y%m%d-%H%M%S"))
reg <- bench_make_registry(dir, cluster)
reg <- bench_submit_cells(reg, cells, refs)
batchtools::waitForJobs(reg = reg)

rows <- bench_collect(reg)
dir.create("bench/results", showWarnings = FALSE, recursive = TRUE)
out <- file.path("bench/results",
                 sprintf("bench-%s-%s.csv", tier, format(Sys.time(), "%Y%m%d-%H%M%S")))
write.csv(rows, out, row.names = FALSE)
message("Written to ", out)

source("bench/collate.R")
cat(bench_report(rows), sep = "\n")
writeLines(bench_report(rows), "bench/results/report.md")
message("Report at bench/results/report.md")
```

- [ ] **Step 5: Add `registry/` to the ignore rules**

In `bench/.gitignore`, add:

```
registry/
```

- [ ] **Step 6: Run it to verify it passes**

Run: `Rscript bench/tests.R`
Expected: PASS. `bench/run.R` itself is exercised in Task 8 once `bench_report()` exists.

- [ ] **Step 7: Format and hand off**

Run: `air format bench/`
Suggested message: `bench: orchestrate cells with batchtools, local or slurm`

---

### Task 8: Collate, report, and append history

**Files:**
- Create: `bench/collate.R`
- Create: `bench/history.csv`
- Test: `bench/tests.R`

**Interfaces:**
- Consumes: `bench_schema()` (Task 6)
- Produces: `bench_deltas(rows, baseline)` returning a `data.frame` with `delta_mem_pct`, `delta_time_pct`, `mem_verdict`, `time_verdict`, `digest_differs`; `bench_report(rows, baseline)` returning `character()` lines; `bench_history_append(rows, path)` returning the appended `data.frame` invisibly

- [ ] **Step 1: Write the failing test**

Append to `bench/tests.R`:

```r
source("bench/collate.R")

.bench_fake <- function(...) {
  base <- data.frame(
    ref = NA_character_, arf_version = "0.2.5.9000", commit = "abc1234",
    op = "lik", backend = "sequential", workers = NA_integer_, n = 1e4, p = 10L,
    trees = 50L, iters = 5L, dt_threads = 1L, ranger_threads = 1L,
    metric = "cgroup-anon+shmem", peak_mb = NA_real_, time_median = NA_real_,
    time_min = NA_real_, time_max = NA_real_, digest = NA_character_,
    digest_kind = "exact", tier = "quick", n_evidence = NA_integer_,
    n_synth = NA_integer_, n_folds = 8L, rowmode = NA_character_,
    host = "h", kernel = "k", r_version = "4.5", job_id = NA_character_,
    timestamp = "2026-10-05T00:00:00", stringsAsFactors = FALSE
  )
  rows <- lapply(list(...), function(o) {
    r <- base
    for (nm in names(o)) r[[nm]] <- o[[nm]]
    r
  })
  do.call(rbind, rows)
}

test_that("a memory win above the threshold reads as real and a small time move does not", {
  rows <- .bench_fake(
    list(ref = "main", peak_mb = 412, time_median = 1.84, digest = "a3f1"),
    list(ref = "HEAD", peak_mb = 371, time_median = 1.79, digest = "a3f1")
  )
  d <- bench_deltas(rows, baseline = "main")
  head_row <- d[d$ref == "HEAD", ]
  expect_equal(round(head_row$delta_mem_pct, 1), -10.0)
  expect_equal(head_row$mem_verdict, "real")
  expect_equal(head_row$time_verdict, "inconclusive")
  expect_false(head_row$digest_differs)
})

test_that("a digest mismatch is flagged in the report", {
  rows <- .bench_fake(
    list(ref = "main", peak_mb = 412, time_median = 1.84, digest = "a3f1"),
    list(ref = "HEAD", peak_mb = 371, time_median = 1.79, digest = "9c2e")
  )
  txt <- paste(bench_report(rows, baseline = "main"), collapse = "\n")
  expect_match(txt, "OUTPUT DIFFERS")
  expect_match(txt, "not apples-to-apples")
})

test_that("a ref identical to the baseline yields no finding and no division by zero", {
  rows <- .bench_fake(
    list(ref = "main", peak_mb = 400, time_median = 2, digest = "a"),
    list(ref = "HEAD", peak_mb = 400, time_median = 2, digest = "a")
  )
  d <- bench_deltas(rows, baseline = "main")
  expect_equal(d$delta_mem_pct[d$ref == "HEAD"], 0)
  expect_equal(d$mem_verdict[d$ref == "HEAD"], "unchanged")
  expect_false(any(is.nan(d$delta_mem_pct)))
})

test_that("a failed baseline measurement does not render as Inf or NaN", {
  rows <- .bench_fake(
    list(ref = "main", peak_mb = 0, time_median = NA_real_, digest = NA_character_),
    list(ref = "HEAD", peak_mb = 371, time_median = 1.79, digest = "a3f1")
  )
  d <- bench_deltas(rows, baseline = "main")
  expect_true(is.na(d$delta_mem_pct[d$ref == "HEAD"]))
  expect_equal(d$mem_verdict[d$ref == "HEAD"], "no baseline")
  txt <- paste(bench_report(rows, baseline = "main"), collapse = "\n")
  expect_false(grepl("Inf|NaN", txt))
})

test_that("history keeps anchor rows only and is append-only", {
  path <- file.path(tempdir(), paste0("hist-", Sys.getpid(), ".csv"))
  on.exit(unlink(path), add = TRUE)
  rows <- .bench_fake(
    list(ref = "HEAD", peak_mb = 371, time_median = 1.79, digest = "a"),
    list(ref = "cran-0.2.5", peak_mb = 412, time_median = 1.84, digest = "b")
  )
  h1 <- bench_history_append(rows, path)
  expect_equal(nrow(h1), 1L)
  expect_equal(h1$ref, "cran-0.2.5")
  h2 <- bench_history_append(rows, path)
  expect_equal(nrow(h2), 2L)
  expect_equal(names(read.csv(path)), bench_schema())
})
```

- [ ] **Step 2: Run it to verify it fails**

Run: `Rscript bench/tests.R`
Expected: FAIL, cannot open `bench/collate.R`.

- [ ] **Step 3: Implement**

Create `bench/collate.R`:

```r
# Deltas, the PR report, and the anchor history. Thresholds come from
# bench/DESIGN.md: memory above 3% is real, time below 10% is inconclusive.

BENCH_MEM_THRESHOLD <- 3
BENCH_TIME_THRESHOLD <- 10

.bench_cell_key <- function(rows) {
  paste(rows$op, rows$backend, rows$workers, rows$n, rows$p, rows$trees,
        rows$rowmode, sep = "|")
}

.bench_pct <- function(new, old) {
  # A baseline of 0 or NA means the baseline measurement failed, which is not a
  # 100% improvement.
  ifelse(is.na(new) | is.na(old) | old <= 0, NA_real_, (new - old) / old * 100)
}

.bench_verdict <- function(pct, threshold) {
  ifelse(is.na(pct), "no baseline",
    ifelse(pct == 0, "unchanged",
      ifelse(abs(pct) >= threshold, "real", "inconclusive")))
}

bench_deltas <- function(rows, baseline = "main") {
  rows$key <- .bench_cell_key(rows)
  base <- rows[rows$ref == baseline, ]
  if (!nrow(base)) {
    stop("baseline ref '", baseline, "' is not in the results", call. = FALSE)
  }
  idx <- match(rows$key, base$key)
  rows$delta_mem_pct <- .bench_pct(rows$peak_mb, base$peak_mb[idx])
  rows$delta_time_pct <- .bench_pct(rows$time_median, base$time_median[idx])
  rows$mem_verdict <- .bench_verdict(rows$delta_mem_pct, BENCH_MEM_THRESHOLD)
  rows$time_verdict <- .bench_verdict(rows$delta_time_pct, BENCH_TIME_THRESHOLD)
  rows$digest_differs <- !is.na(rows$digest) & !is.na(base$digest[idx]) &
    rows$digest != base$digest[idx]
  rows$key <- NULL
  rows
}

.bench_fmt_pct <- function(x) ifelse(is.na(x), "   n/a", sprintf("%+6.1f%%", x))

bench_report <- function(rows, baseline = "main") {
  d <- bench_deltas(rows, baseline)
  d <- d[d$ref != baseline, ]
  metrics <- unique(rows$metric)
  out <- c(
    "# arf benchmark report",
    "",
    sprintf("Baseline: `%s`. Metric: %s. Host: %s, kernel %s, R %s.",
            baseline, paste(metrics, collapse = " / "), rows$host[1],
            rows$kernel[1], rows$r_version[1]),
    sprintf("Memory deltas at or above %g%% are real; time deltas below %g%% are inconclusive.",
            BENCH_MEM_THRESHOLD, BENCH_TIME_THRESHOLD),
    ""
  )
  for (i in seq_len(nrow(d))) {
    r <- d[i, ]
    out <- c(out, sprintf(
      "- %-14s %-10s n=%-6g trees=%-4g %-10s mem %s (%s)  time %s (%s)",
      r$op, r$backend, r$n, r$trees, r$ref,
      .bench_fmt_pct(r$delta_mem_pct), r$mem_verdict,
      .bench_fmt_pct(r$delta_time_pct), r$time_verdict))
    if (isTRUE(r$digest_differs)) {
      out <- c(out, sprintf(
        "    !! OUTPUT DIFFERS (%s digest) against %s: this delta is not apples-to-apples, confirm the change is intended",
        r$digest_kind, baseline))
    }
  }
  c(out, "")
}

# Anchors only: HEAD and main move, so recording them would be noise. Committed
# by hand afterwards, for the runs worth keeping.
bench_history_append <- function(rows, path = "bench/history.csv") {
  keep <- rows[!rows$ref %in% c("HEAD", "main"), bench_schema()]
  if (!nrow(keep)) {
    return(invisible(keep))
  }
  if (file.exists(path)) {
    old <- utils::read.csv(path, stringsAsFactors = FALSE)
    keep <- rbind(old, keep)
  }
  utils::write.csv(keep, path, row.names = FALSE)
  invisible(keep)
}
```

- [ ] **Step 4: Seed the history file**

Create `bench/history.csv` with the header only, so the first append has a schema to match:

```bash
Rscript -e 'source("bench/run-cell.R"); write.csv(setNames(data.frame(matrix(nrow = 0, ncol = length(bench_schema()))), bench_schema()), "bench/history.csv", row.names = FALSE)'
```

- [ ] **Step 5: Run it to verify it passes**

Run: `Rscript bench/tests.R`
Expected: PASS.

- [ ] **Step 6: Smoke-test the whole pipeline**

Run: `ARF_BENCH_TIER=quick Rscript -e 'source("bench/cells.R"); cells <- bench_cells("quick"); message(nrow(cells), " cells")'`
Then run the real thing on a trimmed grid to confirm `bench/run.R` works end to end:

```bash
Rscript -e '
source("bench/bench-helpers.R"); source("bench/refs.R"); source("bench/cells.R")
source("bench/run-cell.R"); source("bench/registry.R"); source("bench/collate.R")
bench_assert_clean_tree()
cells <- bench_cells("quick"); cells <- cells[cells$op == "forde", ][1, ]
cells$n <- 500; cells$trees <- 5L; cells$iters <- 2L
dir <- file.path(tempdir(), "smoke-reg")
reg <- bench_make_registry(dir, "local")
reg <- bench_submit_cells(reg, cells, c("HEAD", "main"))
batchtools::waitForJobs(reg = reg)
rows <- bench_collect(reg)
cat(bench_report(rows), sep = "\n")
'
```

Expected: two rows, one report line, no `Inf` or `NaN`. On a branch where `HEAD` equals `main` the delta is 0% and reads `unchanged`, which is the Review Focus case.

- [ ] **Step 7: Format and hand off**

Run: `air format bench/`
Suggested message: `bench: collate deltas, render the PR report, append anchor history`

---

### Task 9: make target, README, and the trend panel

**Files:**
- Modify: `Makefile`
- Modify: `bench/README.md`
- Modify: `bench/viz.qmd`

**Interfaces:**
- Consumes: `bench/run.R` (Task 7), `bench/history.csv` (Task 8)
- Produces: `make bench` running the quick tier

- [ ] **Step 1: Add the make target**

In `Makefile`, after the `test` target:

```make
.PHONY: bench
bench:
	Rscript bench/run.R
```

- [ ] **Step 2: Verify it runs**

Run: `make bench`
Expected: the quick tier runs locally against `HEAD`, `main` and `cran:0.2.5`, writes `bench/results/bench-quick-*.csv` and `bench/results/report.md`, and prints the report. First run is slower because it installs three refs.

- [ ] **Step 3: Document it in the README**

Add to `bench/README.md`, after the intro:

```markdown
## Version comparison

`bench/run.R` compares refs rather than backends: it installs each ref into its
own library, measures the same grid against all of them in one allocation, and
reports deltas against `main`. See `DESIGN.md` for why the job unit is one cell
against every ref, and why a dirty tree refuses to run.

    make bench                                        # quick tier, local
    ARF_BENCH_TIER=full ARF_BENCH_CLUSTER=slurm \
      ARF_BENCH_SLURM_TMPL=/path/to/slurm.tmpl \
      Rscript bench/run.R                             # full tier, cluster

Refs default to `HEAD`, `main`, and the anchors in `anchors.csv`. Anchor rows
are appended to `history.csv` with `bench_history_append()` and committed by
hand for runs worth keeping.
```

- [ ] **Step 4: Add the trend panel**

Append to `bench/viz.qmd`:

````markdown
## Release trend

```{r}
#| label: history-trend
#| fig-height: 8
hist <- read.csv("history.csv")
if (nrow(hist)) {
  library(ggplot2)
  long <- rbind(
    transform(hist[, c("arf_version", "op", "n", "peak_mb")],
              metric = "peak memory (MB)", value = peak_mb, peak_mb = NULL),
    transform(hist[, c("arf_version", "op", "n", "time_median")],
              metric = "time (s)", value = time_median, time_median = NULL)
  )
  ggplot(long, aes(arf_version, value, group = interaction(op, n))) +
    geom_line() +
    geom_point() +
    facet_grid(metric ~ op, scales = "free_y") +
    labs(x = NULL, y = NULL,
         title = "arf performance across anchors",
         caption = "Anchors are re-measured on one node, so this trend is not stitched across machines.") +
    theme_minimal()
} else {
  cat("No anchor history yet; run bench/run.R and append it.\n")
}
```
````

- [ ] **Step 5: Render the report to check the panel**

Run: `quarto render bench/viz.qmd`
Expected: renders, and with an empty `history.csv` prints the placeholder rather than erroring.

- [ ] **Step 6: Hand off**

Suggested message: `bench: add make bench, document version comparison, plot the anchor trend`
