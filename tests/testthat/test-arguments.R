test_that("FORDE works with alpha>0", {
  arf <- adversarial_rf(iris, parallel = FALSE)
  expect_silent(forde(arf, iris, parallel = FALSE, alpha = 0.01))
})
test_that("mtry default is at least 2 for few features, capped at p", {
  arf2 <- adversarial_rf(iris[, 1:2], num_trees = 2, verbose = FALSE, parallel = FALSE)
  expect_equal(arf2$mtry, 2)
  arf1 <- adversarial_rf(iris[, 1, drop = FALSE], num_trees = 2, verbose = FALSE, parallel = FALSE)
  expect_equal(arf1$mtry, 1)
  arf5 <- adversarial_rf(iris, num_trees = 2, verbose = FALSE, parallel = FALSE)
  expect_equal(arf5$mtry, 2)
  arf_user <- adversarial_rf(iris, mtry = 3, num_trees = 2, verbose = FALSE, parallel = FALSE)
  expect_equal(arf_user$mtry, 3)
})

test_that("FORDE rejects invalid arguments", {
  arf <- adversarial_rf(iris, num_trees = 2, verbose = FALSE, parallel = FALSE)
  expect_error(forde(arf, iris, family = "normal", parallel = FALSE), "family not recognized")
  expect_error(forde(arf, iris, alpha = -1, parallel = FALSE), "alpha must be nonnegative")
  expect_error(forde(arf, iris, epsilon = -1, parallel = FALSE), "epsilon must be nonnegative")
  expect_error(forde(arf, iris[1:10, ], oob = TRUE, parallel = FALSE), "trained on x when oob = TRUE")
  iris_inf <- iris
  iris_inf$Sepal.Length[1] <- Inf
  expect_error(forde(arf, iris_inf, parallel = FALSE), "infinite values")
})

test_that("FORDE with uniform family resets missing finite bounds with a warning", {
  arf <- adversarial_rf(iris, num_trees = 2, verbose = FALSE, parallel = FALSE)
  expect_warning(
    psi <- forde(arf, iris, family = "unif", finite_bounds = "no", parallel = FALSE),
    "Resetting finite_bounds"
  )
  expect_true(all(is.finite(psi$cnt$min)))
})

test_that("LIK rejects unknown query columns", {
  arf <- adversarial_rf(iris, num_trees = 2, verbose = FALSE, parallel = FALSE)
  psi <- forde(arf, iris, parallel = FALSE)
  q <- iris[1:5, ]
  q$bogus <- 1
  expect_error(lik(psi, q, parallel = FALSE), "Unrecognized feature")
})

test_that("FORDE handles purely categorical data", {
  x <- data.frame(
    a = factor(sample(letters[1:3], 100, TRUE)),
    b = factor(sample(c("x", "y"), 100, TRUE))
  )
  arf <- adversarial_rf(x, num_trees = 3, verbose = FALSE, parallel = FALSE)
  psi <- forde(arf, x, parallel = FALSE)
  expect_equal(nrow(psi$cnt), 0)
  expect_true(nrow(psi$cat) > 0)
  expect_equal(nrow(forge(psi, n_synth = 5, parallel = FALSE)), 5)
  expect_length(lik(psi, x[1:10, ], arf = arf, parallel = FALSE), 10)
})
