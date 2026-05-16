# Correctness gate for the experimental mirai+mori backend.
# Prototype-scope: one dataset, equality vs the foreach path.

test_that("mirai backend produces equal forde output on iris", {
  skip_if_not_installed("mirai")
  skip_if_not_installed("mori")

  arf <- adversarial_rf(iris, verbose = FALSE, parallel = FALSE)

  old <- options(arf.backend = "foreach")
  on.exit(options(old), add = TRUE)
  psi_foreach <- forde(arf, iris, parallel = FALSE)

  options(arf.backend = "mirai")
  mirai::daemons(2)
  on.exit(mirai::daemons(0), add = TRUE)
  psi_mirai <- forde(arf, iris, parallel = TRUE)

  expect_equal(psi_mirai$cnt, psi_foreach$cnt, ignore_attr = TRUE)
  expect_equal(psi_mirai$cat, psi_foreach$cat, ignore_attr = TRUE)
  expect_equal(psi_mirai$forest, psi_foreach$forest, ignore_attr = TRUE)
})

test_that("mirai backend errors clearly when daemons are not set", {
  skip_if_not_installed("mirai")
  skip_if_not_installed("mori")

  arf <- adversarial_rf(iris, verbose = FALSE, parallel = FALSE)
  mirai::daemons(0)
  old <- options(arf.backend = "mirai")
  on.exit(options(old), add = TRUE)

  expect_error(forde(arf, iris, parallel = TRUE), "daemons")
})
