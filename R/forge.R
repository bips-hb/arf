#' Forests for Generative Modeling
#' 
#' Uses pre-trained FORDE model to simulate synthetic data.
#' 
#' @param params Circuit parameters learned via \code{\link{forde}}. 
#' @param n_synth Number of synthetic samples to generate.
#' @param evidence Optional set of conditioning events. This can take one of 
#'   three forms: (1) a partial sample, i.e. a single row of data with
#'   some but not all columns; (2) a data frame of conditioning events, 
#'   which allows for inequalities; or (3) a posterior distribution over leaves.
#'   See Details.
#'   
#' @details  
#' \code{forge} simulates a synthetic dataset of \code{n_synth} samples. First,
#' leaves are sampled in proportion to either their coverage (if 
#' \code{evidence = NULL}) or their posterior probability. Then, each feature is 
#' sampled independently within each leaf according to the probability mass or 
#' density function learned by \code{\link{forde}}. This will create realistic 
#' data so long as the adversarial RF used in the previous step satisfies the 
#' local independence criterion. See Watson et al. (2023).
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
#' A dataset of \code{n_synth} synthetic samples. 
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
#' arf <- adversarial_rf(iris)
#' psi <- forde(arf, iris)
#' x_synth <- forge(psi, n_synth = 100)
#'
#' # Condition on Species = "setosa"
#' evi <- data.frame(Species = "setosa")
#' x_synth <- forge(psi, n_synth = 100, evidence = evi)
#' 
#' # Condition in Species = "setosa" and Sepal.Length > 6
#' evi <- data.frame(variable = c("Species", "Sepal.Length"),
#'                   relation = c("==", ">"), 
#'                   value = c("setosa", 6))
#' x_synth <- forge(psi, n_synth = 100, evidence = evi)
#' 
#' # Or just input some distribution on leaves
#' # (Weights that do not sum to unity are automatically scaled)
#' n_leaves <- nrow(psi$forest)
#' evi <- data.frame(f_idx = psi$forest$f_idx, wt = rexp(n_leaves))
#' x_synth <- forge(psi, n_synth = 100, evidence = evi)
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
    evidence = NULL) {
  
  # To avoid data.table check issues
  tree <- cvg <- leaf <- idx <- family <- mu <- sigma <- prob <- dat <- 
    variable <- relation <- wt <- j <- f_idx <- val <- . <- NULL
  
  # Prep evidence
  conj <- FALSE
  to_sim <- params$meta$variable
  if (!is.null(evidence)) {
    evidence <- prep_evi(params, evidence)
    if (!all(c('f_idx', 'wt') %in% colnames(evidence))) {
      conj <- TRUE
      to_sim <- setdiff(params$meta$variable, evidence[relation == '==', variable])
    }
  }
  factor_cols <- params$meta[variable %in% to_sim, family == 'multinom']
  
  # Prepare the event space
  if (is.null(evidence)) {
    num_trees <- params$forest[, max(tree)]
    omega <- params$forest[, .(f_idx, cvg)]
    omega[, wt := cvg / num_trees]
    omega[, cvg := NULL]
  } else if (isTRUE(conj)) {
    omega <- leaf_posterior(params, evidence)
  } else {
    omega <- evidence
  }
  omega <- omega[wt > 0]
  
  if (nrow(omega) == 0) {
    stop("All leaves have zero likelihood.")
  }
  
  # Draw random leaves with probability proportional to weight
  draws <- data.table(
    'f_idx' = omega[, sample(f_idx, size = n_synth, replace = TRUE, prob = wt)]
  )
  omega <- merge(draws, omega, sort = FALSE)[, idx := .I]
  
  # Simulate continuous data
  synth_cnt <- synth_cat <- NULL
  if (any(!factor_cols)) {
    fam <- params$meta[family != 'multinom', unique(family)]
    psi <- merge(omega, params$cnt[variable %in% to_sim], by = 'f_idx', 
                 sort = FALSE, allow.cartesian = TRUE)
    if (isTRUE(conj)) {
      if (any(evidence$relation %in% c('<', '<=', '>', '>='))) {
        for (k in evidence[, which(grepl('<', relation))]) {
          j <- evidence$variable[k]
          value <- as.numeric(evidence$value[k])
          psi[variable == j & max > value, max := value]
        }
        for (k in evidence[, which(grepl('>', relation))]) {
          j <- evidence$variable[k]
          value <- as.numeric(evidence$value[k])
          psi[variable == j & min < value, min := value]
        }
      }
    }
    if (fam == 'truncnorm') {
      psi[, dat := truncnorm::rtruncnorm(.N, a = min, b = max, mean = mu, sd = sigma)]
    } else if (fam == 'unif') {
      psi[, dat := stats::runif(.N, min = min, max = max)]
    }
    synth_cnt <- dcast(psi, idx ~ variable, value.var = 'dat')[, idx := NULL]
  }
  
  # Simulate categorical data
  if (any(factor_cols)) {
    psi <- merge(omega, params$cat[variable %in% to_sim], by = 'f_idx',
                 sort = FALSE, allow.cartesian = TRUE)
    psi[prob == 1, dat := val] 
    if (isTRUE(conj)) {
      if (any(evidence[, relation == '!='])) {
        tmp <- evidence[relation == '!=', .(variable, value)]
        psi <- merge(psi, tmp, by = 'variable', sort = FALSE)
        psi <- psi[val != value]
        psi[, value := NULL]
      }
    }
    psi[prob < 1, dat := sample(val, 1, prob = prob), by = .(variable, idx)]
    psi <- unique(psi[, .(idx, variable, dat)])
    synth_cat <- dcast(psi, idx ~ variable, value.var = 'dat')[, idx := NULL]
  }
  
  # Combine, optionally impose constraint(s)
  x_synth <- cbind(synth_cnt, synth_cat)
  if (length(to_sim) != params$meta[, .N]) {
    tmp <- evidence[relation == '==']
    add_on <- setDT(lapply(tmp[, .I], function(k) rep(tmp[k, value], n_synth)))
    setnames(add_on, tmp[, variable])
    x_synth <- cbind(x_synth, add_on)
  }
  
  # Clean up, export
  x_synth <- post_x(x_synth, params)
  return(x_synth)
}


