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
  mirai::daemons(2)
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
  mirai::daemons(2)
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
  mirai::daemons(2)
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
  mirai::daemons(2)
  on.exit(mirai::daemons(0), add = TRUE)
  x_mirai <- expct(psi, evidence = evi, parallel = TRUE, stepsize = 5, verbose = FALSE)

  expect_equal(dim(x_mirai), dim(x_foreach))
  expect_equal(colnames(x_mirai), colnames(x_foreach))
})
