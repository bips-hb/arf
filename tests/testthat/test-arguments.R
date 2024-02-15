test_that("FORDE works with alpha>0", {
  arf <- adversarial_rf(iris, parallel = FALSE)
  expect_silent(forde(arf, iris, parallel = FALSE, alpha = 0.01))
})