
test_that("FORGE works if y is a column name", {
  # Categorical column
  dat <- iris
  colnames(dat)[5] <- "y"
  arf <- adversarial_rf(dat, parallel = FALSE)
  psi <- forde(arf, dat, parallel = FALSE)
  x_synth <- forge(psi, n_synth = 10, parallel = FALSE)
  expect_in("y", psi$cat$variable)
  expect_equal(sum(is.na(x_synth)), 0) 
  
  # Continuous column
  dat <- iris
  colnames(dat)[1] <- "y"
  arf <- adversarial_rf(dat, parallel = FALSE)
  psi <- forde(arf, dat, parallel = FALSE)
  x_synth <- forge(psi, n_synth = 10, parallel = FALSE)
  expect_in("y", psi$cnt$variable)
  expect_equal(sum(is.na(x_synth)), 0) 
  
  # Categorical column and alpha>0
  dat <- iris
  colnames(dat)[5] <- "y"
  arf <- adversarial_rf(dat, parallel = FALSE)
  psi <- forde(arf, dat, parallel = FALSE, alpha = 1e-5)
  x_synth <- forge(psi, n_synth = 10, parallel = FALSE)
  expect_in("y", psi$cat$variable)
  expect_equal(sum(is.na(x_synth)), 0) 
})


