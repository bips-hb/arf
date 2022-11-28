#' Forests for generative modeling
#' 
#' Uses pre-trained FORDE model to simulate synthetic data.
#' 
#' @param psi Parameters learned via \code{forde}. 
#' @param n_synth Number of synthetic samples to generate.
#' @param family Distribution to use for random sampling. Current options 
#'   include truncated normal (the default \code{family = "truncnorm"}) and 
#'   uniform (\code{family = "unif"}). See Details.
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
#' Currently, \code{forge} only provides support for truncated normal or uniform
#' densities when features are continuous. Future releases will accommodate 
#' a larger class of distributional families.
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
#' fd <- forde(arf, iris)
#' x_synth <- forge(fd$psi, n_synth = 100)
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
    psi, 
    n_synth, 
    family = 'truncnorm',
    parallel = TRUE) {
  
  # Draw random leaves with probability proportional to coverage
  omega <- unique(psi[, .(tree, leaf, cvg)])
  omega[, pr := cvg / max(tree)][, idx := .I]
  draws <- sample(omega$idx, size = n_synth, replace = TRUE, prob = omega$pr)
  psi_fn <- function(i) {
    id <- draws[i]
    out <- psi[tree == omega[idx == id, tree] & leaf == omega[idx == id, leaf]]
    out[, idx := i]
  }
  if (isTRUE(parallel)) {
    psi_idx <- foreach(i = 1:n_synth, .combine = rbind) %dopar% psi_fn(i)
  } else {
    psi_idx <- foreach(i = 1:n_synth, .combine = rbind) %do% psi_fn(i)
  }
  
  # Simulate data
  synth_cnt <- synth_cat <- NULL
  if (nrow(psi[type == 'cnt']) > 0L) {  # Continuous
    psi_cnt <- psi_idx[type == 'cnt']
    if (family == 'truncnorm') {
      psi_cnt[, dat := rtruncnorm(nrow(psi_cnt), a = min, b = max,
                                  mean = mu, sd = sigma)]
    } else if (family == 'unif') {
      psi_cnt[, dat := runif(nrow(psi_cnt), min = min, max = max)]
    }
    synth_cnt <- dcast(psi_cnt, idx ~ variable, value.var = 'dat')
    synth_cnt[, idx := NULL]
  }
  if (nrow(psi[type == 'cat']) > 0L) { # Categorical
    psi_cat <- psi_idx[type == 'cat']
    psi_cat[, dat := sample(cat, 1, prob = prob), by = .(idx, variable)]
    synth_cat <- dcast(unique(psi_cat[, .(idx, variable, dat)]), 
                       idx ~ variable, value.var = 'dat')
    synth_cat[, idx := NULL]
  }
  
  # Export
  x_synth <- cbind(synth_cnt, synth_cat)
  return(x_synth)
}
