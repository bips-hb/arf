test_that("FORDE works with alpha>0", {
  arf <- adversarial_rf(iris, parallel = FALSE)
  expect_silent(forde(arf, iris, parallel = FALSE, alpha = 0.01))
})
test_that("mtry default is at least 2 for few features, capped at p", {
  arf2 <- adversarial_rf(iris[, 1:2], num_trees = 2, verbose = FALSE, parallel = FALSE)
  expect_equal(arf2$mtry, 2)
  arf1 <- adversarial_rf(iris[, 1, drop = FALSE], num_trees = 2, verbose = FALSE, parallel = FALSE)
  expect_equal(arf1$mtry, 1)
  arf5 <- adversarial_rf(iris, num_trees = 2, verbose = FALSE, parallel = FALSE)
  expect_equal(arf5$mtry, 2)
  arf_user <- adversarial_rf(iris, mtry = 3, num_trees = 2, verbose = FALSE, parallel = FALSE)
  expect_equal(arf_user$mtry, 3)
})
