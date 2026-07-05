test_that("sample_from_leaves returns a data.table without params", {
  arf <- adversarial_rf(iris, num_trees = 2, verbose = FALSE, parallel = FALSE)
  x_synth <- sample_from_leaves(arf, iris)
  expect_s3_class(x_synth, "data.table")
  expect_equal(nrow(x_synth), nrow(iris))
  expect_equal(colnames(x_synth), colnames(iris))
})

test_that("sample_from_leaves restores input class and column types with params", {
  # Guard copied from the analogous forge column-types test; 0.16.1 introduced
  # the vector-valued min.bucket adversarial_rf relies on. Left as a pointer in
  # case this can be relaxed for leaf sampling in the future.
  if (utils::packageVersion("ranger") < "0.16.1") {
    skip("can only test this with recent ranger version.")
  }

  n <- 50
  dat <- data.frame(numeric = rnorm(n),
                    integer_factor = sample(1L:5L, n, replace = TRUE),
                    integer_numeric = sample(1L:50L, n, replace = FALSE),
                    character = sample(letters[1:5], n, replace = TRUE),
                    factor = factor(sample(letters[1:5], n, replace = TRUE)),
                    logical = (sample(0:1, n, replace = TRUE) == 1))

  arf <- adversarial_rf(dat, num_trees = 2, verbose = FALSE, parallel = FALSE)
  psi <- forde(arf, dat, parallel = FALSE)
  x_synth <- sample_from_leaves(arf, dat, params = psi)

  # data.frame in, data.frame out
  expect_s3_class(x_synth, "data.frame")
  expect_equal(nrow(x_synth), nrow(dat))

  # No NAs and preserved column types
  expect_true(all(!is.na(x_synth)))
  classes <- sapply(dat, class)
  classes_synth <- sapply(x_synth, class)
  expect_equal(classes, classes_synth)
})

test_that("sample_from_leaves returns a data.table when called with a data.table", {
  dt <- data.table::as.data.table(iris)
  arf <- adversarial_rf(dt, num_trees = 2, verbose = FALSE, parallel = FALSE)
  psi <- forde(arf, dt, parallel = FALSE)
  x_synth <- sample_from_leaves(arf, dt, params = psi)
  expect_s3_class(x_synth, "data.table")
})

test_that("sample_from_leaves returns factors with same levels (and order) with params", {
  arf <- adversarial_rf(iris, num_trees = 2, verbose = FALSE, parallel = FALSE)
  psi <- forde(arf, iris, parallel = FALSE)
  x_synth <- sample_from_leaves(arf, iris, params = psi)
  expect_s3_class(x_synth$Species, "factor")
  expect_equal(levels(x_synth$Species), levels(iris$Species))
})

test_that("sample_from_leaves forwards round to post-processing", {
  # See note on the guard above: 0.16.1 introduced vector-valued min.bucket.
  if (utils::packageVersion("ranger") < "0.16.1") {
    skip("can only test this with recent ranger version.")
  }
  # Leaf sampling reuses observed values, which already sit at the real data's
  # precision, so rounding only shows up where a type is coerced: an
  # integer-valued numeric column stays numeric with round = FALSE and becomes
  # integer with round = TRUE (mirrors forge()).
  n <- 50
  dat <- data.frame(numeric = rnorm(n),
                    integer_numeric = sample(1L:50L, n, replace = FALSE))
  arf <- adversarial_rf(dat, num_trees = 2, verbose = FALSE, parallel = FALSE)
  psi <- forde(arf, dat, parallel = FALSE)

  x_round <- sample_from_leaves(arf, dat, params = psi, round = TRUE)
  x_noround <- sample_from_leaves(arf, dat, params = psi, round = FALSE)

  expect_equal(class(x_round$integer_numeric), "integer")
  expect_equal(class(x_noround$integer_numeric), "numeric")
})

test_that("sample_from_leaves only draws values present in the real data", {
  # Marginal intra-leaf sampling reuses observed values, so every synthetic
  # value must appear in the corresponding real column.
  arf <- adversarial_rf(iris, num_trees = 2, verbose = FALSE, parallel = FALSE)
  psi <- forde(arf, iris, parallel = FALSE)
  x_synth <- sample_from_leaves(arf, iris, params = psi, round = FALSE)
  for (j in colnames(iris)) {
    expect_true(all(x_synth[[j]] %in% iris[[j]]))
  }
})
