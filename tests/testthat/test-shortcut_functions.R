# Shortcut wrappers: train + forde + lik/forge/expct in one call.
iris_small <- iris[c(1:20, 51:70, 101:120), ]

test_that("darf returns one log-likelihood per query row", {
  ll <- darf(iris_small, num_trees = 5, parallel = FALSE)
  expect_type(ll, "double")
  expect_length(ll, nrow(iris_small))
  ll_q <- darf(iris_small, query = iris_small[1:7, ], num_trees = 5, parallel = FALSE)
  expect_length(ll_q, 7)
})

test_that("rarf defaults to nrow(x) draws, or one per evidence row", {
  x <- rarf(iris_small, num_trees = 5, parallel = FALSE)
  expect_s3_class(x, "data.frame")
  expect_equal(dim(x), dim(iris_small))
  expect_identical(names(x), names(iris_small))
  evi <- data.frame(Species = c("setosa", "virginica"))
  x_evi <- rarf(iris_small, evidence = evi, num_trees = 5, parallel = FALSE)
  expect_equal(nrow(x_evi), 2)
  expect_identical(as.character(x_evi$Species), c("setosa", "virginica"))
  x3 <- rarf(iris_small, n_synth = 3, num_trees = 5, parallel = FALSE)
  expect_equal(nrow(x3), 3)
})

test_that("earf returns expectations, unconditional and conditional", {
  e <- earf(iris_small, num_trees = 5, parallel = FALSE)
  expect_equal(nrow(e), 1)
  expect_identical(names(e), names(iris_small))
  e_evi <- earf(
    iris_small,
    evidence = data.frame(Species = "setosa"),
    query = "Sepal.Length",
    num_trees = 5,
    parallel = FALSE
  )
  expect_identical(names(e_evi), "Sepal.Length")
})

test_that("shortcuts accept pre-computed params and skip training", {
  arf <- adversarial_rf(iris_small, num_trees = 5, verbose = FALSE, parallel = FALSE)
  psi <- forde(arf, iris_small, parallel = FALSE)
  expect_equal(nrow(rarf(iris_small, n_synth = 4, params = psi, parallel = FALSE)), 4)
  expect_warning(
    ll <- darf(iris_small, params = psi, parallel = FALSE),
    "faster to include the pre-trained arf"
  )
  expect_length(ll, nrow(iris_small))
  expect_length(darf(iris_small, arf = arf, params = psi, parallel = FALSE), nrow(iris_small))
  expect_equal(nrow(earf(iris_small, params = psi, parallel = FALSE)), 1)
  # a pre-trained arf alone skips training but still runs forde
  expect_equal(nrow(rarf(iris_small, n_synth = 2, arf = arf, parallel = FALSE)), 2)
})
