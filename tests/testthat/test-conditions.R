
arf <- adversarial_rf(iris, parallel = FALSE)
psi <- forde(arf, iris, parallel = FALSE)

test_that("numeric conditions work", {
  synth <- forge(psi, evidence = data.frame(Sepal.Length = "5"), n_synth = 10, parallel = FALSE)
  expect_equal(synth$Sepal.Length, rep(5, 10))
})

test_that("same result with different conditioning syntax", {
  set.seed(10)
  synth1 <- forge(psi, evidence = data.frame(Sepal.Length = ">5", Species = "setosa"), n_synth = 10, parallel = FALSE)
  set.seed(10)
  synth2 <- forge(psi, evidence = data.frame(Sepal.Length = "(5,Inf)", Species = "setosa"), n_synth = 10, parallel = FALSE)
  expect_equal(synth1, synth2) #all(synth1 == synth2)
})

test_that("categorical conditions work", {
  synth <- forge(psi, evidence = data.frame(Species = "setosa"), n_synth = 10, parallel = FALSE)
  expect_equal(as.character(synth$Species), rep("setosa", 10))
})

test_that("categorical conditions work with not", {
  synth <- forge(psi, evidence = data.frame(Species = "!setosa"), n_synth = 10, parallel = FALSE)
  expect_in(as.character(synth$Species), c("versicolor", "virginica"))
})

test_that("if nomatch='force' and verbose=TRUE, run through with a warning", {
  # Zero likelihood case (no finite bounds)
  psi_no <- forde(arf, iris, finite_bounds = "no", parallel = FALSE)
  expect_warning(x_synth <- forge(psi_no, evidence = data.frame(Sepal.Length = 100), 
                                  nomatch = "force", verbose = TRUE,
                                  n_synth = 10, parallel = FALSE), 
      "All leaves have zero likelihood for some entered evidence rows\\. This is probably because evidence contains an \\(almost\\) impossible combination\\. Sampling from all possible leaves with equal probability \\(can be changed with 'nomatch' argument\\)\\.")
  expect_true(all(!is.na(x_synth) & x_synth$Sepal.Length == 100))
  
  # No matching leaf case (finite bounds)
  psi_global <- forde(arf, iris, finite_bounds = "global", parallel = FALSE)
  expect_warning(x_synth <- forge(psi_global, evidence = data.frame(Sepal.Length = 100), 
                                  nomatch = "force", verbose = TRUE,
                                  n_synth = 10, parallel = FALSE), 
                 "For some entered evidence rows, no matching leaves could be found\\. This is probably because evidence lies outside of the distribution calculated by FORDE\\. For continuous data, consider setting epsilon>0 or finite_bounds='no' in forde\\(\\)\\. For categorical data, consider setting alpha>0 in forde\\(\\)\\. Sampling from all leaves with equal probability \\(can be changed with 'nomatch' argument\\)\\.")
  expect_true(all(!is.na(x_synth) & x_synth$Sepal.Length == 100))
})

test_that("if nomatch='force' and verbose=FALSE, run through without a warning", {
  # Zero likelihood case (no finite bounds)
  psi_no <- forde(arf, iris, finite_bounds = "no", parallel = FALSE)
  expect_silent(x_synth <- forge(psi_no, evidence = data.frame(Sepal.Length = 100), 
                                 nomatch = "force", verbose = FALSE,
                                 n_synth = 10, parallel = FALSE))
  expect_true(all(!is.na(x_synth) & x_synth$Sepal.Length == 100))
  
  # No matching leaf case (finite bounds)
  psi_global <- forde(arf, iris, finite_bounds = "global", parallel = FALSE)
  expect_silent(x_synth <- forge(psi_global, evidence = data.frame(Sepal.Length = 100), 
                                 nomatch = "force", verbose = FALSE,
                                 n_synth = 10, parallel = FALSE))
  expect_true(all(!is.na(x_synth) & x_synth$Sepal.Length == 100))
})

test_that("if nomatch='na' and verbose=TRUE, run through with a warning and return NA", {
  # Zero likelihood case (no finite bounds)
  psi_no <- forde(arf, iris, finite_bounds = "no", parallel = FALSE)
  expect_warning(x_synth <- forge(psi_no, evidence = data.frame(Sepal.Length = 100), 
                                  nomatch = "na", verbose = TRUE, 
                                  n_synth = 10, parallel = FALSE), 
                 "All leaves have zero likelihood for some entered evidence rows\\. This is probably because evidence contains an \\(almost\\) impossible combination\\. Returning NA for those rows \\(can be changed with 'nomatch' argument\\)\\.")
  expect_true(all(is.na(x_synth[, -1]) & x_synth$Sepal.Length == 100))
  
  # No matching leaf case (finite bounds)
  psi_global <- forde(arf, iris, finite_bounds = "global", parallel = FALSE)
  expect_warning(x_synth <- forge(psi_global, evidence = data.frame(Sepal.Length = 100), 
                                  nomatch = "na", verbose = TRUE, 
                                  n_synth = 10, parallel = FALSE), 
                 "For some entered evidence rows, no matching leaves could be found\\. This is probably because evidence lies outside of the distribution calculated by FORDE\\. For continuous data, consider setting epsilon>0 or finite_bounds='no' in forde\\(\\)\\. For categorical data, consider setting alpha>0 in forde\\(\\)\\. Returning NA for those rows \\(can be changed with 'nomatch' argument\\)\\.")
  expect_true(all(is.na(x_synth[, -1]) & x_synth$Sepal.Length == 100))
})

test_that("if nomatch='na' and verbose=FALSE, run through without a warning and return NA", {
  # Zero likelihood case (no finite bounds)
  psi_no <- forde(arf, iris, finite_bounds = "no", parallel = FALSE)
  expect_silent(x_synth <- forge(psi_no, evidence = data.frame(Sepal.Length = 100), 
                                 nomatch = "na", verbose = FALSE, 
                                 n_synth = 10, parallel = FALSE))
  expect_true(all(is.na(x_synth[, -1]) & x_synth$Sepal.Length == 100))
  
  # No matching leaf case (finite bounds)
  psi_global <- forde(arf, iris, finite_bounds = "global", parallel = FALSE)
  expect_silent(x_synth <- forge(psi_global, evidence = data.frame(Sepal.Length = 100), 
                                 nomatch = "na", verbose = FALSE, 
                                 n_synth = 10, parallel = FALSE))
  expect_true(all(is.na(x_synth[, -1]) & x_synth$Sepal.Length == 100))
})

