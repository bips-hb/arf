# Gate for every daemon-backed test. Call it as the FIRST line of the test:
# the arf training and forde() calls that set these tests up are expensive, and
# a skip further down would pay for them on CRAN before bailing.
# skip_on_cran() is required, not just belt-and-braces: mirai/mori are in
# Suggests, so skip_if_not_installed() does NOT skip on CRAN, and daemon
# processes brush against CRAN's 2-core policy (a known source of flakiness).
skip_if_no_daemons <- function() {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("mirai")
  testthat::skip_if_not_installed("mori")
}

# Set up mirai daemons for backend tests. In CI the package is installed, so the
# package's own arf_load_on_daemons() (loadNamespace) suffices. Under load_all
# (devtools::test()) daemons can't see the source build, so dev-load arf on them
# here; the package's loadNamespace call is then a no-op.
setup_mirai_daemons <- function(n = 2, seed = NULL) {
  if (is.null(seed)) {
    mirai::daemons(n)
  } else {
    mirai::daemons(n, seed = seed)
  }
  if (requireNamespace("pkgload", quietly = TRUE) &&
      isTRUE(tryCatch(pkgload::is_dev_package("arf"), error = function(e) FALSE))) {
    pdir <- pkgload::pkg_path()
    mirai::everywhere(suppressMessages(pkgload::load_all(pdir, quiet = TRUE)), pdir = pdir)
  }
}
