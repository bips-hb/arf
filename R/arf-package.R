#' @seealso
#' \code{\link{adversarial_rf}}, \code{\link{forde}}, \code{\link{forge}}, \code{\link{expct}}, \code{\link{lik}}
#' 
#' Useful links:
#' \itemize{
#'   \item \url{https://github.com/bips-hb/arf}
#'   \item \url{https://bips-hb.github.io/arf/}
#'   \item Report bugs at \url{https://github.com/bips-hb/arf/issues}
#' }
#' @examples
#' # Train ARF and estimate leaf parameters
#' arf <- adversarial_rf(iris)
#' psi <- forde(arf, iris)
#' 
#' # Generate 100 synthetic samples from the iris dataset
#' x_synth <- forge(psi, n_synth = 100)
#'
#' # Condition on Species = "setosa" and Sepal.Length > 6
#' evi <- data.frame(Species = "setosa",
#'                   Sepal.Length = "(6, Inf)")
#' x_synth <- forge(psi, n_synth = 100, evidence = evi)
#' 
#' # Estimate average log-likelihood
#' ll <- lik(psi, iris, arf = arf, log = TRUE)
#' mean(ll)
#' 
#' # Expectation of Sepal.Length for class setosa
#' evi <- data.frame(Species = "setosa")
#' expct(psi, query = "Sepal.Length", evidence = evi)
#' 
#' \dontrun{
#' # Parallelization with doParallel
#' doParallel::registerDoParallel(cores = 4)
#'
#' # ... or with doFuture
#' doFuture::registerDoFuture()
#' future::plan("multisession", workers = 4)
#' }
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
