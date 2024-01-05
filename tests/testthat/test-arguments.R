test_that("FORDE works with alpha>0", {
  arf <- adversarial_rf(iris)
  expect_silent(forde(arf, iris, alpha = 0.01))
})