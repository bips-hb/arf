test_that("nomatch = 'force' fills impossible evidence rows in separate mode", {
  set.seed(42)
  n <- 200
  dat <- data.frame(
    x1 = runif(n, 0, 10),
    x2 = runif(n, 0, 10),
    z = factor(sample(c("a", "b"), n, TRUE))
  )
  arf <- adversarial_rf(dat, num_trees = 20, verbose = FALSE, parallel = FALSE)
  psi <- forde(arf, dat, finite_bounds = "local", parallel = FALSE)
  # Row 1 matches leaves; rows 2/3 lie outside all leaf bounds
  evidence <- data.frame(x1 = c(5, 999, -5), z = factor(c(NA, "a", NA), levels = c("a", "b")))
  evidence$x2 <- c(NA, NA, 999)

  expect_warning(
    x_synth <- forge(psi, n_synth = 4, evidence = evidence, evidence_row_mode = "separate", parallel = FALSE),
    "no matching leaves"
  )
  expect_equal(nrow(x_synth), 12)
  expect_false(anyNA(x_synth))
  expect_s3_class(x_synth$z, "factor")
  expect_equal(x_synth$x1[5:8], rep(999, 4))
  expect_equal(as.character(x_synth$z[5:8]), rep("a", 4))
  expect_equal(x_synth$x2[9:12], rep(999, 4))

  expect_warning(
    e <- expct(psi, evidence = evidence, evidence_row_mode = "separate", parallel = FALSE),
    "no matching leaves"
  )
  expect_equal(nrow(e), 3)
  expect_false(anyNA(e))
})

test_that("nomatch = 'force' works with data.table input", {
  set.seed(42)
  n <- 200
  dat <- data.table::data.table(x1 = runif(n, 0, 10), x2 = runif(n, 0, 10))
  arf <- adversarial_rf(dat, num_trees = 20, verbose = FALSE, parallel = FALSE)
  psi <- forde(arf, dat, finite_bounds = "local", parallel = FALSE)
  evidence <- data.table::data.table(x1 = c(5, 999))

  expect_warning(
    x_synth <- forge(psi, n_synth = 3, evidence = evidence, evidence_row_mode = "separate", parallel = FALSE),
    "no matching leaves"
  )
  expect_s3_class(x_synth, "data.table")
  expect_equal(nrow(x_synth), 6)
  expect_false(anyNA(x_synth))

  expect_warning(
    e <- expct(psi, evidence = evidence, evidence_row_mode = "separate", parallel = FALSE),
    "no matching leaves"
  )
  expect_equal(nrow(e), 2)
  expect_false(anyNA(e))
  expect_named(e, "x2")
})
