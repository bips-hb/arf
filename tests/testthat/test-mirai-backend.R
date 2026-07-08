# Correctness gate for the mirai+mori backend:
# equality vs the sequential/foreach paths on one dataset.

test_that("mirai backend produces equal forde output on iris", {
  skip_if_not_installed("mirai")
  skip_if_not_installed("mori")

  arf <- adversarial_rf(iris, verbose = FALSE, parallel = FALSE)

  old <- options(arf.backend = "foreach")
  on.exit(options(old), add = TRUE)
  psi_foreach <- forde(arf, iris, parallel = FALSE)

  options(arf.backend = "mirai")
  setup_mirai_daemons(2)
  on.exit(mirai::daemons(0), add = TRUE)
  psi_mirai <- forde(arf, iris, parallel = TRUE)

  expect_equal(psi_mirai$cnt, psi_foreach$cnt, ignore_attr = TRUE)
  expect_equal(psi_mirai$cat, psi_foreach$cat, ignore_attr = TRUE)
  expect_equal(psi_mirai$forest, psi_foreach$forest, ignore_attr = TRUE)
})

test_that("mirai backend errors clearly when daemons are not set", {
  skip_if_not_installed("mirai")
  skip_if_not_installed("mori")

  arf <- adversarial_rf(iris, verbose = FALSE, parallel = FALSE)
  mirai::daemons(0)
  old <- options(arf.backend = "mirai")
  on.exit(options(old), add = TRUE)

  expect_error(forde(arf, iris, parallel = TRUE), "daemons")
})

test_that("mirai backend produces structurally consistent forge() output", {
  skip_if_not_installed("mirai")
  skip_if_not_installed("mori")

  # forge() is stochastic, so compare structure (not values) across backends.
  arf <- adversarial_rf(iris, verbose = FALSE, parallel = FALSE)
  psi <- forde(arf, iris, parallel = FALSE)
  # 21 separate exact conditions (mixed species) -> multi-step
  evi <- iris[c(1:7, 51:57, 101:107), "Species", drop = FALSE]

  old <- options(arf.backend = "foreach")
  on.exit(options(old), add = TRUE)
  x_foreach <- forge(psi, n_synth = 3, evidence = evi, parallel = FALSE,
                     stepsize = 5, verbose = FALSE)

  options(arf.backend = "mirai")
  setup_mirai_daemons(2)
  on.exit(mirai::daemons(0), add = TRUE)
  x_mirai <- forge(psi, n_synth = 3, evidence = evi, parallel = TRUE,
                   stepsize = 5, verbose = FALSE)

  expect_equal(nrow(x_mirai), nrow(x_foreach))
  expect_equal(colnames(x_mirai), colnames(x_foreach))
  expect_equal(sapply(x_mirai, class), sapply(x_foreach, class))
  # exact conditions must map onto their output rows in evidence order:
  # catches scrambled step-to-row assembly that structural checks miss
  expect_equal(as.character(x_mirai$Species),
               as.character(rep(evi$Species, each = 3)))
})

test_that("mirai backend gives identical lik() (deterministic)", {
  skip_if_not_installed("mirai")
  skip_if_not_installed("mori")

  arf <- adversarial_rf(iris, verbose = FALSE, parallel = FALSE)
  psi <- forde(arf, iris, parallel = FALSE)

  old <- options(arf.backend = "foreach")
  on.exit(options(old), add = TRUE)
  ll_foreach <- lik(psi, iris, arf = arf, batch = 30, parallel = FALSE)

  options(arf.backend = "mirai")
  setup_mirai_daemons(2)
  on.exit(mirai::daemons(0), add = TRUE)
  ll_mirai <- lik(psi, iris, arf = arf, batch = 30, parallel = TRUE)

  expect_equal(ll_mirai, ll_foreach)  # lik is deterministic
})

test_that("mirai backend gives structurally consistent expct()", {
  skip_if_not_installed("mirai")
  skip_if_not_installed("mori")

  arf <- adversarial_rf(iris, verbose = FALSE, parallel = FALSE)
  psi <- forde(arf, iris, parallel = FALSE)
  evi <- iris[1:20, "Species", drop = FALSE]  # 20 separate conditions -> multi-step

  old <- options(arf.backend = "foreach")
  on.exit(options(old), add = TRUE)
  x_foreach <- expct(psi, evidence = evi, parallel = FALSE, stepsize = 5, verbose = FALSE)

  options(arf.backend = "mirai")
  setup_mirai_daemons(2)
  on.exit(mirai::daemons(0), add = TRUE)
  x_mirai <- expct(psi, evidence = evi, parallel = TRUE, stepsize = 5, verbose = FALSE)

  expect_equal(dim(x_mirai), dim(x_foreach))
  expect_equal(colnames(x_mirai), colnames(x_foreach))
})

test_that("arf_n_workers reflects active mirai daemons (stepsize sizing)", {
  skip_if_not_installed("mirai")
  skip_on_cran()

  mirai::daemons(0)
  n_idle <- arf_n_workers()  # no mirai, no foreach -> 1
  expect_equal(n_idle, 1L)

  mirai::daemons(3)
  on.exit(mirai::daemons(0), add = TRUE)
  expect_gte(arf_n_workers(), 3L)  # must see daemons, else step_no==1 kills mirai
})

test_that("mirai backend gives identical cforde() (deterministic)", {
  skip_if_not_installed("mirai")
  skip_if_not_installed("mori")

  arf <- adversarial_rf(iris, verbose = FALSE, parallel = FALSE)
  psi <- forde(arf, iris, parallel = FALSE)
  set.seed(1)
  evi <- data.frame(Sepal.Length = runif(20, 4.5, 7))  # 20 conditions -> multi-step

  old <- options(arf.backend = "foreach")
  on.exit(options(old), add = TRUE)
  cf_foreach <- arf:::cforde(psi, evi, stepsize = 5, parallel = FALSE, verbose = FALSE)

  options(arf.backend = "mirai")
  setup_mirai_daemons(2)
  on.exit(mirai::daemons(0), add = TRUE)
  cf_mirai <- arf:::cforde(psi, evi, stepsize = 5, parallel = TRUE, verbose = FALSE)

  expect_equal(cf_mirai$forest, cf_foreach$forest)  # cforde is deterministic
  expect_equal(cf_mirai$cnt, cf_foreach$cnt)
  expect_equal(cf_mirai$cat, cf_foreach$cat)
})

test_that("mirai backend preserves row order and class in expct() (regression)", {
  skip_if_not_installed("mirai")
  skip_if_not_installed("mori")

  # step_no (8) > n_workers (2): guards against interleaved-chunk row scrambling
  # and data.table-vs-data.frame class drift in arf_mirai_tree_map.
  arf <- adversarial_rf(iris, verbose = FALSE, parallel = FALSE)
  psi <- forde(arf, iris, parallel = FALSE)
  evi <- data.frame(Sepal.Length = seq(4.5, 7, length.out = 8))

  old <- options(arf.backend = "foreach")
  on.exit(options(old), add = TRUE)
  x_foreach <- expct(psi, query = "Petal.Length", evidence = evi,
                     parallel = FALSE, stepsize = 1, verbose = FALSE)

  options(arf.backend = "mirai")
  setup_mirai_daemons(2)
  on.exit(mirai::daemons(0), add = TRUE)
  x_mirai <- expct(psi, query = "Petal.Length", evidence = evi,
                   parallel = TRUE, stepsize = 1, verbose = FALSE)

  expect_equal(x_mirai, x_foreach)  # exact: order + values + class + row.names
})

test_that("mirai backend gives identical adversarial_rf() pruning (deterministic)", {
  skip_if_not_installed("mirai")
  skip_if_not_installed("mori")

  # ranger draws per-tree seeds from the R RNG, so training is reproducible
  # under set.seed() regardless of threading; prune is deterministic given a
  # forest. That makes the full run comparable end to end, exercising the
  # actual mirai prune dispatch in adversarial_rf().
  set.seed(7)
  a_serial <- adversarial_rf(iris, num_trees = 50, parallel = FALSE,
                             verbose = FALSE)

  old <- options(arf.backend = "mirai")
  on.exit(options(old), add = TRUE)
  setup_mirai_daemons(2)
  on.exit(mirai::daemons(0), add = TRUE)
  set.seed(7)
  a_mirai <- adversarial_rf(iris, num_trees = 50, parallel = TRUE,
                            verbose = FALSE)

  expect_identical(a_mirai$forest$child.nodeIDs, a_serial$forest$child.nodeIDs)
})

test_that("mirai backend runs adversarial_rf() end to end", {
  skip_if_not_installed("mirai")
  skip_if_not_installed("mori")

  old <- options(arf.backend = "mirai")
  on.exit(options(old), add = TRUE)
  setup_mirai_daemons(2)
  on.exit(mirai::daemons(0), add = TRUE)

  a <- adversarial_rf(iris, num_trees = 30, parallel = TRUE, verbose = FALSE)
  expect_s3_class(a, "ranger")
  expect_length(a$forest$child.nodeIDs, a$num.trees)
})

test_that("arf.block_rows caps expct blocks without changing results", {
  a <- adversarial_rf(iris, num_trees = 10, parallel = FALSE, verbose = FALSE)
  psi <- forde(a, iris, parallel = FALSE)
  evi <- data.frame(Species = sample(levels(iris$Species), 6, replace = TRUE))
  evi1 <- evi[1, , drop = FALSE]

  set.seed(1)
  ref <- expct(psi, evidence = evi, parallel = FALSE)
  set.seed(2)
  ref1 <- expct(psi, evidence = evi1, parallel = FALSE)

  old <- options(arf.block_rows = 1)  # force one condition per block
  on.exit(options(old), add = TRUE)
  set.seed(1)
  blocked <- expct(psi, evidence = evi, parallel = FALSE)
  expect_identical(blocked, ref)

  # single condition: cannot be split, must take the unblocked path
  set.seed(2)
  blocked1 <- expct(psi, evidence = evi1, parallel = FALSE)
  expect_identical(blocked1, ref1)
})

test_that("arf_tree_chunks yields contiguous in-order blocks covering all trees", {
  for (nt in c(1L, 5L, 7L)) {
    for (nw in c(1L, 3L, 10L)) {
      ch <- arf_tree_chunks(nt, nw)
      # concatenating chunks in chunk order must reproduce 1..nt exactly:
      # this is the contiguity invariant positional ops rely on
      expect_identical(unlist(ch, use.names = FALSE), seq_len(nt))
      expect_lte(length(ch), max(1L, min(nw, nt)))
    }
  }
})

test_that("arf_select_backend applies the documented precedence", {
  skip_if_not_installed("mirai")
  skip_if_not_installed("mori")
  skip_on_cran()

  old <- options(arf.backend = NULL, arf.verbose = FALSE)
  on.exit(options(old), add = TRUE)

  expect_identical(arf_select_backend(FALSE), "sequential")

  mirai::daemons(0)
  expect_identical(arf_select_backend(TRUE), "foreach")  # no daemons, option unset

  setup_mirai_daemons(2)
  on.exit(mirai::daemons(0), add = TRUE)
  expect_identical(arf_select_backend(TRUE), "mirai")  # daemons auto-detected

  options(arf.backend = "foreach")
  expect_identical(arf_select_backend(TRUE), "foreach")  # explicit option wins

  options(arf.backend = "bogus")
  expect_error(arf_select_backend(TRUE))
})

test_that("arf_load_on_daemons caches per daemon pool and self-invalidates", {
  skip_if_not_installed("mirai")

  setup_mirai_daemons(2)
  on.exit(mirai::daemons(0), add = TRUE)
  assign("arf_loaded_key", NULL, envir = arf:::.arf_env)

  expect_true(arf_load_on_daemons())    # first call loads
  expect_false(arf_load_on_daemons())   # same pool -> cached

  mirai::daemons(0)
  setup_mirai_daemons(2)                # rebuilt pool mints a new key
  expect_true(arf_load_on_daemons())
})

test_that("mirai worker errors propagate instead of corrupting results", {
  skip_if_not_installed("mirai")
  skip_if_not_installed("mori")

  setup_mirai_daemons(2)
  on.exit(mirai::daemons(0), add = TRUE)

  boom <- function(tree) stop("boom")
  environment(boom) <- globalenv()
  expect_error(arf_mirai_tree_map(4, boom, list()), "boom")
})

test_that("foreach doParallel backend gives equal forde output", {
  skip_if_not_installed("doParallel")
  skip_on_cran()

  a <- adversarial_rf(iris, verbose = FALSE, parallel = FALSE)
  psi_seq <- forde(a, iris, parallel = FALSE)

  # PSOCK, not fork: forking after mirai/nanonext threads exist is unsafe
  cl <- parallel::makeCluster(2)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  if (requireNamespace("pkgload", quietly = TRUE) &&
      isTRUE(tryCatch(pkgload::is_dev_package("arf"), error = function(e) FALSE))) {
    pdir <- pkgload::pkg_path()
    parallel::clusterCall(cl, function(p) {
      suppressMessages(pkgload::load_all(p, quiet = TRUE))
    }, pdir)
  }
  doParallel::registerDoParallel(cl)
  on.exit(foreach::registerDoSEQ(), add = TRUE)
  old <- options(arf.backend = "foreach")
  on.exit(options(old), add = TRUE)

  psi_par <- forde(a, iris, parallel = TRUE)

  expect_equal(psi_par$cnt, psi_seq$cnt, ignore_attr = TRUE)
  expect_equal(psi_par$cat, psi_seq$cat, ignore_attr = TRUE)
  expect_equal(psi_par$forest, psi_seq$forest, ignore_attr = TRUE)
})
