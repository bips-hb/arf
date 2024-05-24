
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

