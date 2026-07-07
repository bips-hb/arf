# Correctness gate for the experimental mirai+mori backend.
# Prototype-scope: one dataset, equality vs the foreach path.

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
  evi <- iris[1:20, "Species", drop = FALSE]  # 20 separate conditions -> multi-step

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
  expect_true(all(as.character(x_mirai$Species) %in% levels(iris$Species)))
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

  # Prune is deterministic given a forest. Build one unpruned forest, then prune
  # it both ways and compare (can't compare full adversarial_rf across backends:
  # ranger threading makes the forest itself non-reproducible).
  set.seed(7)
  a0 <- adversarial_rf(iris, num_trees = 50, prune = FALSE, parallel = FALSE,
                       verbose = FALSE)
  pred <- stats::predict(a0, prep_x(iris), type = "terminalNodes")$predictions + 1L
  nt <- 50L
  serial <- lapply(seq_len(nt), arf_prune_tree,
                   child_nodeIDs = a0$forest$child.nodeIDs, pred = pred,
                   min_node_size = 2L)

  setup_mirai_daemons(2)
  on.exit(mirai::daemons(0), add = TRUE)
  child_shared <- mori::share(a0$forest$child.nodeIDs)
  pred_shared <- mori::share(pred)
  n_chunks <- max(1L, min(as.integer(mirai::status()$connections), nt))
  chunks <- split(seq_len(nt), sort(rep(seq_len(n_chunks), length.out = nt)))
  res <- mirai::mirai_map(
    chunks,
    function(trees, worker, child_nodeIDs, pred, min_node_size) {
      lapply(trees, worker, child_nodeIDs = child_nodeIDs, pred = pred,
             min_node_size = min_node_size)
    },
    .args = list(worker = arf_prune_tree, child_nodeIDs = child_shared,
                 pred = pred_shared, min_node_size = 2L))[]
  mirai_out <- unname(do.call(c, res))

  expect_identical(mirai_out, serial)
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

test_that("arf.chunk_factor > 1 leaves results unchanged", {
  skip_if_not_installed("mirai")
  skip_if_not_installed("mori")

  a <- adversarial_rf(iris, num_trees = 20, parallel = FALSE, verbose = FALSE)
  serial <- forde(a, iris, parallel = FALSE)

  old <- options(arf.backend = "mirai", arf.chunk_factor = 4)
  on.exit(options(old), add = TRUE)
  on.exit(options(arf.chunk_factor = NULL), add = TRUE)
  setup_mirai_daemons(2)
  on.exit(mirai::daemons(0), add = TRUE)

  # finer chunks split trees 8 ways on 2 daemons; output must not change
  psi <- forde(a, iris, parallel = TRUE)
  expect_equal(psi, serial)

  a_mirai <- adversarial_rf(iris, num_trees = 20, parallel = TRUE,
                            verbose = FALSE)
  expect_length(a_mirai$forest$child.nodeIDs, 20L)
})

test_that("arf.block_rows caps expct blocks without changing results", {
  a <- adversarial_rf(iris, num_trees = 10, parallel = FALSE, verbose = FALSE)
  psi <- forde(a, iris, parallel = FALSE)
  evi <- data.frame(Species = sample(levels(iris$Species), 6, replace = TRUE))

  set.seed(1)
  ref <- expct(psi, evidence = evi, parallel = FALSE)

  old <- options(arf.block_rows = 1)  # force one condition per block
  on.exit(options(old), add = TRUE)
  set.seed(1)
  blocked <- expct(psi, evidence = evi, parallel = FALSE)

  expect_identical(blocked, ref)
})
