# Set up mirai daemons for backend tests. In CI the package is installed, so the
# package's own arf_load_on_daemons() (loadNamespace) suffices. Under load_all
# (devtools::test()) daemons can't see the source build, so dev-load arf on them
# here; the package's loadNamespace call is then a no-op.
setup_mirai_daemons <- function(n = 2, seed = NULL) {
  # mirai/mori are in Suggests, so skip_if_not_installed() does NOT skip on
  # CRAN; skip here instead -- daemon processes would brush against CRAN's
  # 2-core policy and are a known source of flakiness on its builders.
  testthat::skip_on_cran()
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
