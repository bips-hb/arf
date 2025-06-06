
arf <- adversarial_rf(iris, parallel = FALSE)
psi <- forde(arf, iris, parallel = FALSE)

test_that("expct returns correct values", {
  # For all classes
  res <- expct(psi, query = "Sepal.Length", parallel = FALSE)
  expect_equal(colnames(res), "Sepal.Length")
  expect_equal(res$Sepal.Length, mean(iris$Sepal.Length))
  
  # Only for setosa
  res <- expct(psi, query = "Sepal.Length", evidence = data.frame(Species = "setosa"), parallel = FALSE)
  expect_equal(res$Sepal.Length, mean(iris[iris$Species == "setosa", "Sepal.Length"]), tolerance = .1)
})

test_that("expct works for vectorized evidence", {
  evi <- data.frame(Species = c("setosa", "versicolor", "virginica"))
  res <- expct(psi, query = "Petal.Width", evidence = evi, parallel = FALSE)
  expect_equal(res$Petal.Width, c(0.2, 1.3, 2.0), tolerance = .2)
})

test_that("expct works with NAs", {
  evi <- iris[1:3,]
  evi[1, 1] <- NA_real_
  evi[1, 5] <- NA_character_
  evi[2, 2] <- NA_real_
  res <- expct(psi, evidence = evi, parallel = FALSE)
  
  # Result should be the same but with expected values for NAs
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), nrow(evi))
  expect_equal(ncol(res), ncol(evi))
  expect_equal(colnames(res), colnames(evi))
  expect_equal(res[1, 2], evi[1, 2])
})

test_that("expct works with partial sample", {
  # Petal.* missing
  evi <- iris[1:3, c(1:2, 5)]
  res <- expct(psi, evidence = evi, parallel = FALSE)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), nrow(evi))
  expect_equal(ncol(res), ncol(iris) - ncol(evi))
  expect_equal(colnames(res), colnames(iris)[3:4])
  
  # Petal.* and Species missing
  evi <- iris[1:3, 1:2]
  res <- expct(psi, evidence = evi, parallel = FALSE)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), nrow(evi))
  expect_equal(ncol(res), ncol(iris) - ncol(evi))
  expect_equal(colnames(res), colnames(iris)[3:5])
})

test_that("if nomatch='force', do not return NA", {
  # Zero likelihood case (no finite bounds)
  psi_no <- forde(arf, iris, finite_bounds = "no", parallel = FALSE)
  x_synth <- expct(psi_no, evidence = data.frame(Sepal.Length = 100), 
                   nomatch = "force", verbose = FALSE, 
                   parallel = FALSE)
  expect_true(all(!is.na(x_synth) & x_synth$Sepal.Length == 100))
  
  # No matching leaf case (finite bounds)
  psi_global <- forde(arf, iris, finite_bounds = "global", parallel = FALSE)
  x_synth <- expct(psi_global, evidence = data.frame(Sepal.Length = 100), 
                   nomatch = "force", verbose = FALSE, 
                   parallel = FALSE)
  expect_true(all(!is.na(x_synth) & x_synth$Sepal.Length == 100))
})

test_that("if nomatch='na', return NA", {
  # Zero likelihood case (no finite bounds)
  psi_no <- forde(arf, iris, finite_bounds = "no", parallel = FALSE)
  x_synth <- expct(psi_no, evidence = data.frame(Sepal.Length = 100), 
                   nomatch = "na", verbose = FALSE, 
                   parallel = FALSE)
  expect_true(all(is.na(x_synth[, -1]) & x_synth$Sepal.Length == 100))
  
  # No matching leaf case (finite bounds)
  psi_global <- forde(arf, iris, finite_bounds = "global", parallel = FALSE)
  x_synth <- expct(psi_global, evidence = data.frame(Sepal.Length = 100), 
                   nomatch = "na", verbose = FALSE, 
                   parallel = FALSE)
  expect_true(all(is.na(x_synth[, -1]) & x_synth$Sepal.Length == 100))
})
