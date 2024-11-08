
# Generate some missings
iris_na <- iris
for (j in 1:ncol(iris)) {
  iris_na[sample(1:nrow(iris), 5), j] <- NA
}
  
test_that("impute returns same data with message if no missings", {
  expect_message(iris_imputed <- arf::impute(iris, parallel = FALSE), "No missing values found\\. Returning input data\\.")
  expect_equal(iris, iris_imputed)
})

test_that("Imputation fills missing values", {
  if (utils::packageVersion("ranger") < "0.16.4") {
    skip("can only test this with recent ranger version.")
  }
  
  # Single imputation
  iris_imputed <- arf::impute(iris_na, m = 1, parallel = FALSE)
  expect_s3_class(iris_imputed, "data.frame")
  expect_true(!anyNA(iris_imputed))
  
  # Multiple imputation
  iris_imputed <- arf::impute(iris_na, parallel = FALSE)
  expect_type(iris_imputed, "list")
  expect_length(iris_imputed, 20)
  expect_true(all(sapply(iris_imputed, function(x) !anyNA(x))))  
})
