test_that("ARF returns ranger object", {
  arf <- adversarial_rf(iris, num_trees = 2, verbose = FALSE, parallel = FALSE)
  expect_s3_class(arf, "ranger")
})

test_that("FORDE returns correct list object", {
  arf <- adversarial_rf(iris, num_trees = 2, verbose = FALSE, parallel = FALSE)
  psi <- forde(arf, iris, parallel = FALSE)
  expect_type(psi, "list")
  expect_named(psi, c("cnt", "cat", "forest", "meta", "input_class"))
  expect_s3_class(psi$cnt, "data.table")
  expect_s3_class(psi$cat, "data.table")
  expect_s3_class(psi$forest, "data.table")
  expect_s3_class(psi$meta, "data.table")
  expect_type(psi$input_class, "character")
})

test_that("FORDE categories sum to unity", {
  arf <- adversarial_rf(iris, num_trees = 2, verbose = FALSE, parallel = FALSE)
  psi <- forde(arf, iris, parallel = FALSE)
  tmp <- psi$cat[, sum(prob), by = f_idx]
  expect_true(all.equal(tmp$V1, rep(1, times = tmp[, .N])))
})

test_that("Likelihood calculation returns vector of log-likelihoods", {
  arf <- adversarial_rf(iris, num_trees = 2, verbose = FALSE, parallel = FALSE)
  psi <- forde(arf, iris, parallel = FALSE)
  loglik <- lik(psi, iris, arf = arf, parallel = FALSE)
  loglik2 <- lik(psi, iris, parallel = FALSE)
  expect_type(loglik, "double")
  expect_length(loglik, nrow(iris))
  expect_true(all(!is.na(loglik)))
  expect_true(all.equal(loglik, loglik2))
})

test_that("FORGE returns data frame when called with data frame", {
  arf <- adversarial_rf(iris, num_trees = 2, verbose = FALSE, parallel = FALSE)
  psi <- forde(arf, iris, parallel = FALSE)
  x_synth <- forge(psi, n_synth = 20, parallel = FALSE)
  expect_s3_class(x_synth, "data.frame")
})

test_that("FORGE returns data table when called with data table", {
  arf <- adversarial_rf(data.table::as.data.table(iris), num_trees = 2, verbose = FALSE, parallel = FALSE)
  psi <- forde(arf, data.table::as.data.table(iris), parallel = FALSE)
  x_synth <- forge(psi, n_synth = 20, parallel = FALSE)
  expect_s3_class(x_synth, "data.table")
})

test_that("FORGE returns matrix when called with matrix", {
  n <- 50
  p <- 4
  x <- matrix(runif(n * p), ncol = p)
  
  arf <- adversarial_rf(x, num_trees = 2, verbose = FALSE, parallel = FALSE)
  psi <- forde(arf, x, parallel = FALSE)
  x_synth <- forge(psi, n_synth = 20, parallel = FALSE)
  expect_type(x_synth, "double")
  expect_true(is.matrix(x_synth))
})

test_that("FORGE returns same column types", {
  n <- 50
  dat <- data.frame(numeric = rnorm(n), 
                    integer = sample(1L:5L, n, replace = TRUE), 
                    character = sample(letters[1:5], n, replace = TRUE), 
                    factor = factor(sample(letters[1:5], n, replace = TRUE)), 
                    logical = (sample(0:1, n, replace = TRUE) == 1))
  
  expect_warning(arf <- adversarial_rf(dat, num_trees = 2, verbose = FALSE, parallel = FALSE))
  expect_warning(psi <- forde(arf, dat, parallel = FALSE))
  x_synth <- forge(psi, n_synth = 20, parallel = FALSE)
  
  # No NAs
  expect_true(all(!is.na(x_synth)))

  # Keeps column types
  types <- sapply(dat, typeof)
  types_synth <- sapply(x_synth, typeof)
  expect_equal(types, types_synth)
})
