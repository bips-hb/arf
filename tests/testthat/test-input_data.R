
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

test_that("Imputation works if integers recoded to factors", {
  if (utils::packageVersion("ranger") >= "0.16.1") {
    skip("can only test this with recent ranger version.")
  }
  
  dat <- iris
  dat$int <- as.integer(sample(1:4, nrow(iris), replace = TRUE))
  dat[1:50, "int"] <- NA_integer_
  dat[30:40, "Sepal.Length"] <- NA_real_
  
  expect_warning(arf <- adversarial_rf(dat, parallel = FALSE))
  psi <- forde(arf, dat, parallel = FALSE)
  x_synth <- forge(psi, n_synth = 10, parallel = FALSE)
  expect_in("int", psi$cat$variable)
  expect_equal(sum(is.na(x_synth)), 0) 
  
  # y as column name
  dat <- iris
  dat$y <- as.integer(sample(1:4, nrow(iris), replace = TRUE))
  dat[1:50, "y"] <- NA_integer_
  dat[30:40, "Sepal.Length"] <- NA_real_
  
  expect_warning(arf <- adversarial_rf(dat, parallel = FALSE))
  psi <- forde(arf, dat, parallel = FALSE)
  x_synth <- forge(psi, n_synth = 10, parallel = FALSE)
  expect_in("y", psi$cat$variable)
  expect_equal(sum(is.na(x_synth)), 0) 
})


