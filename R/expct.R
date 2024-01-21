#' Expected Value
#' 
#' Compute the expectation of some query variable(s), optionally conditioned
#' on some event(s).
#' 
#' @param params Circuit parameters learned via \code{\link{forde}}. 
#' @param query Optional character vector of variable names. Estimates will be
#'   computed for each. If \code{NULL}, all variables other than those in 
#'   \code{evidence} will be estimated.
#' @param evidence Optional set of conditioning events. This can take one of 
#'   three forms: (1) a partial sample, i.e. a single row of data with some but
#'   not all columns; (2) a data frame of conditioning events, which allows for 
#'   inequalities; or (3) a posterior distribution over leaves. See Details.
#'   
#'   
#' @details 
#' This function computes expected values for any subset of features, optionally 
#' conditioned on some event(s). 
#' 
#' 
#' @return 
#' A one row data frame with values for all query variables.
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
#' # Train ARF and corresponding circuit
#' arf <- adversarial_rf(iris)
#' psi <- forde(arf, iris)
#' 
#' # What is the expected value Sepal.Length?
#' expct(psi, query = "Sepal.Length")
#' 
#' # What if we condition on Species = "setosa"?
#' evi <- data.frame(Species = "setosa")
#' expct(psi, query = "Sepal.Length", evidence = evi)
#' 
#' # Compute expectations for all features other than Species
#' expct(psi, evidence = evi)
#' 
#' 
#' @seealso
#' \code{\link{adversarial_rf}}, \code{\link{forde}}, \code{\link{lik}}
#' 
#'
#' @export
#' @import data.table
#' @importFrom truncnorm etruncnorm
#' 

expct <- function(
    params, 
    query = NULL, 
    evidence = NULL) {
  
  # To avoid data.table check issues
  variable <- tree <- f_idx <- cvg <- wt <- V1 <- value <- val <- family <-
    mu <- sigma <- obs <- prob <- . <- NULL
  
  # Prep evidence
  conj <- FALSE
  if (!is.null(evidence)) {
    evidence <- prep_evi(params, evidence)
    if (!all(c('f_idx', 'wt') %in% colnames(evidence))) {
      conj <- TRUE
    }
  } 
  
  # Check query
  if (is.null(query)) {
    if (isTRUE(conj)) {
      query <- setdiff(params$meta$variable, evidence$variable)
    } else {
      query <- params$meta$variable
      if (!is.null(evidence)) {
        warning('Computing expectations for all variables. To avoid this ',
                'for conditioning variables, consider passing evidence in the ',
                'form of a partial sample or data frame of events.')
      }
    }
  } else if (any(!query %in% params$meta$variable)) {
    err <- setdiff(query, params$meta$variable)
    stop('Unrecognized feature(s) in query: ', err)
  }
  factor_cols <- params$meta[variable %in% query, family == 'multinom']
  
  # PMF over leaves
  if (is.null(evidence)) {
    num_trees <- params$forest[, max(tree)]
    omega <- params$forest[, .(f_idx, cvg)]
    omega[, wt := cvg / num_trees]
    omega[, cvg := NULL]
  } else if (conj) {
    omega <- leaf_posterior(params, evidence)
  } else {
    omega <- evidence
  }
  omega <- omega[wt > 0]
  
  psi_cnt <- psi_cat <- NULL
  # Continuous data
  if (any(!factor_cols)) {
    tmp <- merge(params$cnt[variable %in% query], omega, by = 'f_idx', sort = FALSE)
    # tmp[, expct := truncnorm::etruncnorm(min, max, mu, sigma)]
    # psi_cnt <- tmp[, crossprod(wt, expct), by = variable]
    psi_cnt <- tmp[, crossprod(wt, mu), by = variable]
    psi_cnt <- dcast(psi_cnt, . ~ variable, value.var = 'V1')[, . := NULL]
  }
  
  # Categorical data
  if (any(factor_cols)) {
    tmp <- merge(params$cat[variable %in% query], omega, by = 'f_idx', sort = FALSE)
    tmp <- tmp[, crossprod(prob, wt), by = .(variable, val)]
    tmp <- tmp[order(match(variable, query[factor_cols]))]
    vals <- tmp[tmp[, .I[which.max(V1)], by = variable]$V1]$val
    psi_cat <- setDT(lapply(seq_along(vals), function(j) vals[j]))
    setnames(psi_cat, query[factor_cols])
  }
  
  # Clean up, export
  out <- cbind(psi_cnt, psi_cat)
  out <- post_x(out, params)
  return(out)
}




