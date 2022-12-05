#' Forests for Generative Modeling
#' 
#' Uses pre-trained FORDE model to simulate synthetic data.
#' 
#' @param params Parameters learned via \code{forde}. 
#' @param n_synth Number of synthetic samples to generate.
#' @param parallel Compute in parallel? Must register backend beforehand, e.g. 
#'   via \code{doParallel}.
#'   
#' @details  
#' \code{forge} simulates a synthetic dataset of \code{n_synth} samples. First,
#' leaves are sampled in proportion to their coverage. Then, each feature is
#' sampled independently within each leaf according to the probability mass or
#' density function learned by \code{\link{forde}}. This will create realistic
#' data so long as the random forest used in the previous step satisfies the 
#' local independence criterion. 
#' 
#' 
#' @return  
#' A dataset of \code{n_synth} synthetic samples. Because continuous and 
#' categorical features are treated separately, columns may not be in the same
#' order as the input data passed to \code{adversarial_rf}.
#' 
#' 
#' @references 
#' Watson, D., Blesch, K., Kapar, J., & Wright, M. (2022). Adversarial random 
#' forests for density estimation and generative modeling. \emph{arXiv} preprint,
#' 2205.09435.
#'
#'
#' @examples
#' arf <- adversarial_rf(iris)
#' psi <- forde(arf, iris)
#' x_synth <- forge(psi, n_synth = 100)
#'
#' @seealso
#' \code{\link{adversarial_rf}}, \code{\link{forde}}
#' 
#' 
#' @export
#' @import data.table
#' @importFrom foreach foreach %do% %dopar%
#' @importFrom truncnorm rtruncnorm 
#' 

forge <- function(
    params, 
    n_synth, 
    parallel = TRUE) {
  
  # To avoid data.table check issues
  tree <- cvg <- leaf <- idx <- family <- mu <- sigma <- prob <- dat <- variable <- j <- . <- NULL
  
  # Draw random leaves with probability proportional to coverage
  num_trees <- params[, max(tree)]
  omega <- unique(params[, .(tree, leaf, cvg)])[, idx := .I]
  draws <- data.table('idx' = sample(omega$idx, size = n_synth, replace = TRUE, 
                                     prob = omega$cvg / num_trees))
  omega <- merge(draws, omega, sort = FALSE)[, idx := .I]
  params_idx <- merge(omega, params, by = c('tree', 'leaf', 'cvg'), sort = FALSE, 
                      allow.cartesian = TRUE)
  
  # Simulate data
  synth_cnt <- synth_cat <- NULL
  if (params[family != 'multinom', .N > 0L]) {  # Continuous
    params_cnt <- params_idx[family != 'multinom']
    fams <- params_cnt[, unique(family)]
    if ('truncnorm' %in% fams) {
      params_cnt[family == 'truncnorm', 
                 dat := truncnorm::rtruncnorm(.N, a = min, b = max, mean = mu, sd = sigma)]
    } 
    if ('unif' %in% fams) {
      params_cnt[family == 'unif', dat := stats::runif(.N, min = min, max = max)]
    }
    synth_cnt <- dcast(params_cnt, idx ~ variable, value.var = 'dat')
    synth_cnt[, idx := NULL]
  }
  if (params[family == 'multinom', .N > 0L]) {  # Categorical
    params_idx[prob == 1, dat := cat]
    synth_cat_fn <- function(j) {
      params_j <- params_idx[variable == j]
      params_j[prob < 1, dat := sample(cat, 1, prob = prob), by = idx]
      out <- data.table(unique(params_j[, .(idx, dat)])[, dat])
      colnames(out) <- j
      return(out)
    }
    cat_vars <- params[family == 'multinom', unique(variable)]
    if (isTRUE(parallel) & length(cat_vars) > 1) {
      synth_cat <- foreach(j = cat_vars, .combine = cbind) %dopar% synth_cat_fn(j)
    } else {
      synth_cat <- foreach(j = cat_vars, .combine = cbind) %do% synth_cat_fn(j)
    }
  }
  
  # Export
  x_synth <- cbind(synth_cnt, synth_cat)
  return(x_synth)
}


