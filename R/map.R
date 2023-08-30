#' Maximum A Posteriori Estimation
#' 
#' Compute the most likely value of some query variable(s), optionally 
#' conditioned on some event(s).
#' 
#' @param params Circuit parameters learned via \code{\link{forde}}. 
#' @param query Optional character vector of variable names. MAP estimates will 
#'   be computed for each. If \code{NULL}, all variables other than those 
#'   in \code{evidence} will be estimated.
#' @param evidence Optional set of conditioning events. This can take one of 
#'   three forms: (1) a partial sample, i.e. a single row of data with some but
#'   not all columns; (2) a data frame of conditioning events, which allows for 
#'   inequalities; or (3) a posterior distribution over leaves. See Details.
#' @param n_eval Number of points to use for grid search.
#'   
#'   
#' @details 
#' This function computes maximum a posteriori (MAP) values for any subset of
#' features, optionally conditioned on some event(s). This technically covers
#' a range of related query types, including most probable explanations (MPE) 
#' and marginal MAP, which are special instances of MAP. 
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
#' # What is the most likely Sepal.Length?
#' map(psi, query = "Sepal.Length")
#' 
#' # What if we condition on Species = "setosa"?
#' evi <- data.frame(Species = "setosa")
#' map(psi, query = "Sepal.Length", evidence = evi)
#' 
#' # Compute MAP estimates for all features other than Species
#' map(psi, evidence = evi)
#' 
#' 
#' @seealso
#' \code{\link{adversarial_rf}}, \code{\link{forde}}, \code{\link{lik}}
#' 
#'
#' @export
#' @import data.table
#' @importFrom truncnorm dtruncnorm 
#' 

map <- function(
    params, 
    query = NULL, 
    evidence = NULL, 
    n_eval = 100) {
  
  # To avoid data.table check issues
  variable <- family <- tree <- f_idx <- cvg <- wt <- V1 <- value <- val <- 
    mu <- sigma <- obs <- prob <- fold <- . <- NULL
  
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
        warning('Computing MAP estimates for all variables. To avoid this ',
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
  # Continuous data...
  if (any(!factor_cols)) {
    # Grid search for continuous features...?
    # Take the union of the bounds of all weighted likelihood-maximizing leaves
    # for each tree-variable combo
    tmp <- merge(params$cnt[variable %in% query], omega, by = 'f_idx', sort = FALSE)
    tmp <- merge(tmp, params$forest[, .(f_idx, tree)], by = 'f_idx', sort = FALSE)
    tmp[, wt_lik := wt * truncnorm::dtruncnorm(mu, min, max, mu, sigma)]
    min_j <- tmp[is.finite(min), min(min), by = variable]
    max_j <- tmp[is.finite(max), max(max), by = variable]
    tmp[, new_min := min][, new_max := max]
    for (j in query[!factor_cols]) {
      tmp[!is.finite(min) & variable == j, new_min := min_j[variable == j, V1]]
      tmp[!is.finite(max) & variable == j, new_max := max_j[variable == j, V1]]
    }
    tmp <- tmp[tmp[, .I[which.max(wt_lik)], by = .(tree, variable)]$V1]
    min_j <- tmp[, min(new_min), by = variable]
    max_j <- tmp[, max(new_max), by = variable]
    x <- rbindlist(lapply(which(!factor_cols), function(j) {
      data.table(
        'obs' = seq_len(n_eval), 'variable' = query[j], 
        'value' = seq(min_j[variable == query[j], V1], 
                      max_j[variable == query[j], V1], length.out = n_eval)
      )
    }))
    tmp <- merge(tmp, x, by = 'variable', sort = FALSE, allow.cartesian = TRUE)
    tmp[, lik := truncnorm::dtruncnorm(value, min, max, mu, sigma)]
    tmp[value == min, lik := 0]
    tmp <- tmp[, crossprod(lik, wt), by = .(obs, variable)]
    tmp <- tmp[order(match(variable, query[!factor_cols]))]
    mle_dt <- tmp[tmp[, .I[which.max(V1)], by = variable]$V1]
    psi_cnt <- setDT(lapply(which(!factor_cols), function(j) {
      mle_idx <- mle_dt[variable == query[j], obs]
      x[variable == query[j] & obs == mle_idx, value]
    }))
    setnames(psi_cnt, query[!factor_cols])
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


# Allow evidence to be a data frame, and return one answer per row
# Need to consider the possibility of shared omega's across instances


