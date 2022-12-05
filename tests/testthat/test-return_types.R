test_that("FORGE returns data frame when called with data frame", {
  arf <- adversarial_rf(iris, num_trees = 2)
  psi <- forde(arf, iris)
  x_synth <- forge(psi, n_synth = 20)
  expect_s3_class(x_synth, "data.frame")
})
