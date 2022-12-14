#' Forests for Generative Modeling
#' 
#' Uses pre-trained FORDE model to simulate synthetic data.
#' 
#' @param params Parameters learned via \code{\link{forde}}. 
#' @param n_synth Number of synthetic samples to generate.
#' @param parallel Compute in parallel? Must register backend beforehand, e.g. 
#'   via \code{doParallel}.
#'   
#' @details  
#' \code{forge} simulates a synthetic dataset of \code{n_synth} samples. First,
#' leaves are sampled in proportion to their coverage. Then, each feature is
#' sampled independently within each leaf according to the probability mass or
#' density function learned by \code{\link{forde}}. This will create realistic
#' data so long as the adversarial RF used in the previous step satisfies the 
#' local independence criterion. See Watson et al. (2022).
#' 
#' 
#' @return  
#' A dataset of \code{n_synth} synthetic samples. 
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
  omega <- params$forest
  draws <- data.table(
    'f_idx' = omega[, sample(f_idx, size = n_synth, replace = TRUE, prob = cvg)]
  )
  omega <- merge(draws, omega, sort = FALSE)[, idx := .I]
  
  # Simulate data
  synth_cnt <- synth_cat <- NULL
  if (!is.null(params$cnt)) {
    fam <- params$meta[family != 'multinom', unique(family)]
    psi <- merge(omega, params$cnt, by = 'f_idx', sort = FALSE, allow.cartesian = TRUE)
    if (fam == 'truncnorm') {
      psi[, dat := truncnorm::rtruncnorm(.N, a = min, b = max, mean = mu, sd = sigma)]
    } else if (fam == 'unif') {
      psi[, dat := stats::runif(.N, min = min, max = max)]
    }
    synth_cnt <- dcast(psi, idx ~ variable, value.var = 'dat')[, idx := NULL]
  }
  if (!is.null(params$cat)) {
    sim_cat <- function(j) {
      psi <- merge(omega, params$cat[variable == j], by = 'f_idx', sort = FALSE, 
                   allow.cartesian = TRUE)
      psi[prob == 1, dat := val]
      psi[prob < 1, dat := sample(val, 1, prob = prob), by = idx]
      setnames(data.table(unique(psi[, .(idx, dat)])[, dat]), j)
    }
    cat_vars <- params$meta[family == 'multinom', variable]
    if (isTRUE(parallel)) {
      synth_cat <- foreach(j = cat_vars, .combine = cbind) %dopar% sim_cat(j)
    } else {
      synth_cat <- foreach(j = cat_vars, .combine = cbind) %do% sim_cat(j)
    }
  }
  
  # Clean up, export
  x_synth <- cbind(synth_cnt, synth_cat)
  setcolorder(x_synth, params$meta$variable)
  setDF(x_synth)
  idx_factor <- params$meta[, which(class == 'factor')]
  idx_logical <- params$meta[, which(class == 'logical')]
  idx_integer <- params$meta[, which(class == 'integer')]
  if (sum(idx_factor) > 0L) {
    x_synth[, idx_factor] <- as.data.frame(
      lapply(x_synth[, idx_factor, drop = FALSE], as.factor)
    )
  }
  if (sum(idx_logical) > 0L) {
    x_synth[, idx_logical] <- as.data.frame(
      lapply(x_synth[, idx_logical, drop = FALSE], function(x) {x == "TRUE"})
    )
  }
  if (sum(idx_integer) > 0L) {
    x_synth[, idx_integer] <- as.data.frame(
      lapply(x_synth[, idx_integer, drop = FALSE], function(x) as.integer(levels(x))[x]) 
    )
  }
  if ("data.table" %in% params$input_class) {
    x_synth <- as.data.table(x_synth)
  } else if ("matrix" %in% params$input_class) {
    x_synth <- as.matrix(x_synth)
  }
  return(x_synth)
}


