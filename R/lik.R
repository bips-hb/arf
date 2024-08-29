#' Likelihood Estimation
#' 
#' Compute the likelihood of input data, optionally conditioned on some event(s).
#' 
#' @param params Circuit parameters learned via \code{\link{forde}}. 
#' @param query Data frame of samples, optionally comprising just a subset of 
#'   training features. Likelihoods will be computed for each sample. Missing
#'   features will be marginalized out. See Details.
#' @param evidence Optional set of conditioning events. This can take one of 
#'   three forms: (1) a partial sample, i.e. a single row of data with
#'   some but not all columns; (2) a data frame of conditioning events, 
#'   which allows for inequalities and intervals, where multiple rows are 
#'   combined by logical OR; or (3) a posterior distribution over leaves;
#'   see Details and Examples.
#' @param arf Pre-trained \code{\link{adversarial_rf}} or other object of class 
#'   \code{ranger}. This is not required but speeds up computation considerably
#'   for total evidence queries. (Ignored for partial evidence queries.)
#' @param log Return likelihoods on log scale?
#' @param stepsize Stepsize defining number of query and evidence (if provided)
#'  rows handled in one for each step.
#'  Defaults to nrow(evidence)/num_registered_workers for \code{parallel == TRUE}.
#' @param parallel Compute in parallel? Must register backend beforehand, e.g. 
#'   via \code{doParallel}.
#'   
#'   
#' @details 
#' This function computes the likelihood of input data, optionally conditioned
#' on some event(s). Queries may be partial, i.e. covering some but not all
#' features, in which case excluded variables will be marginalized out. 
#' 
#' There are three methods for (optionally) encoding conditioning events via the 
#' \code{evidence} argument. The first is to provide a partial sample, where
#' some but not all columns from the training data are present. The second is to 
#' provide a data frame with three columns: \code{variable}, \code{relation}, 
#' and \code{value}. This supports inequalities via \code{relation}. 
#' Alternatively, users may directly input a pre-calculated posterior 
#' distribution over leaves, with columns \code{f_idx} and \code{wt}. This may 
#' be preferable for complex constraints. See Examples.
#' 
#' 
#' @return 
#' A vector of likelihoods, optionally on the log scale. 
#' 
#' 
#' @references 
#' Watson, D., Blesch, K., Kapar, J., & Wright, M. (2023). Adversarial random 
#' forests for density estimation and generative modeling. In \emph{Proceedings 
#' of the 26th International Conference on Artificial Intelligence and 
#' Statistics}, pp. 5357-5375.
#' 
#' 
#' @examples
#' # Train ARF and estimate leaf parameters
#' arf <- adversarial_rf(iris)
#' psi <- forde(arf, iris)
#' 
#' # Estimate average log-likelihood
#' ll <- lik(psi, iris, arf = arf, log = TRUE)
#' mean(ll)
#' 
#' # Identical but slower
#' ll <- lik(psi, iris, log = TRUE)
#' mean(ll)
#' 
#' # Partial evidence query
#' lik(psi, query = iris[1, 1:3])
#' 
#' # Condition on Species = "setosa"
#' evi <- data.frame(Species = "setosa")
#' lik(psi, query = iris[1, 1:3], evidence = evi)
#' 
#' # Condition on Species = "setosa" and Petal.Width > 0.3
#' evi <- data.frame(Species = "setosa", 
#'                   Petal.Width = ">0.3")
#' lik(psi, query = iris[1, 1:3], evidence = evi)
#' 
#' @seealso
#' \code{\link{arf}}, \code{\link{adversarial_rf}}, \code{\link{forde}}, \code{\link{forge}}, \code{\link{expct}}
#' 
#'
#' @export
#' @import data.table
#' 

lik <- function(
    params, 
    query,
    evidence = NULL,
    arf = NULL,
    log = TRUE, 
    stepsize = 0, 
    parallel = TRUE) {
  
  # To avoid data.table check issues
  c_idx <- c_idx_extract <- cvg <- wt <- NULL
  
  # Check query
  x <- as.data.frame(query)
  colnames_x <- colnames(x)
  n <- nrow(x)
  d <- ncol(x)
  if (d == params$meta[, .N] & is.null(arf)) {
    warning('For total evidence queries, it is faster to include the ', 
            'pre-trained arf.')
  }
  if (any(!colnames(x) %in% params$meta$variable)) {
    err <- setdiff(colnames(x), params$meta$variable)
    stop('Unrecognized feature(s) among colnames: ', err)
  }
  
  # Prepare evidence
  if (!is.null(evidence)) {
    evidence <- as.data.table(evidence)
    if (ncol(evidence) == 2 && all(colnames(evidence) == c("f_idx", "wt"))) {
      cparams_full <- copy(params)
      cparams_full$forest <- copy(evidence)[, wt := wt/sum(wt)]
    } else {
      cparams_full <- cforde(params = params, evidence = evidence, row_mode = "or", output = "params", nomatch = "na_warning", stepsize = stepsize, parallel = parallel)
    }
  } else {
    cparams_full <- copy(params)
  }
  
  if (nrow(cparams_full$forest) == 0) {
    lik <- NA
  } else {
    lik <- cforde(params = cparams_full, evidence = query, row_mode = "separate", output = "p", stepsize = stepsize, parallel = parallel)
  
    # Export
    if (log) {
      lik <- log(lik)
    }
  }
  
  lik
}
