#' Forests for Generative Modeling
#' 
#' Uses pre-trained FORDE model to simulate synthetic data.
#' 
#' @param params Parameters learned via \code{\link{forde}}. 
#' @param n_synth Number of synthetic samples to generate.
#' @param evidence Optional data frame of conditioning event(s) or posterior 
#'   distribution over leaves. See Details.
#' @param parallel Compute in parallel? Must register backend beforehand, e.g. 
#'   via \code{doParallel}.
#'   
#' @details  
#' \code{forge} simulates a synthetic dataset of \code{n_synth} samples. First,
#' leaves are sampled in proportion to their coverage. Then, each feature is
#' sampled independently within each leaf according to the probability mass or
#' density function learned by \code{\link{forde}}. This will create realistic
#' data so long as the adversarial RF used in the previous step satisfies the 
#' local independence criterion. See Watson et al. (2023).
#' 
#' There are two methods for (optionally) encoding conditioning events via the 
#' \code{evidence} argument. The first is to provide a data frame with three 
#' columns: \code{variable}, \code{operator}, and \code{value}. Each row will be 
#' treated as a separate conjunct. Alternatively, users may directly input a 
#' pre-calculated posterior distribution over leaves, with columns \code{f_idx} 
#' and \code{wt}. This may be preferable for complex constraints. See Examples.
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
#' evi <- data.frame(variable = "Species", operator = "==", value = "setosa")
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
#' @importFrom tibble as_tibble
#' 

forge <- function(
    params, 
    n_synth, 
    evidence = NULL,
    parallel = TRUE) {
  
  # To avoid data.table check issues
  tree <- cvg <- leaf <- idx <- family <- mu <- sigma <- prob <- dat <- 
    variable <- j <- f_idx <- val <- . <- NULL
  
  # Check evidence
  if (!is.null(evidence)) {
    evidence <- as.data.table(evidence)
    conj <- all(c('variable', 'operator', 'value') %in% colnames(evidence))
    post <- all(c('f_idx', 'wt') %in% colnames(evidence))
    if (conj + post != 1L) {
      stop('evidence must either be a data frame of conjuncts or a posterior
           distribution on leaves.')
    }
  }
  
  # Prepare the event space
  if (is.null(evidence)) {
    omega <- params$forest
    omega[, wt := cvg / max(tree)]
    omega <- omega[, .(f_idx, wt)]
  } else if (any(grepl('variable', colnames(evidence)))) {
    omega <- leaf_posterior(params, evidence)
  } else {
    omega <- evidence
  }
  
  # Draw random leaves with probability proportional to weight
  draws <- data.table(
    'f_idx' = omega[, sample(f_idx, size = n_synth, replace = TRUE, prob = wt)]
  )
  omega <- merge(draws, omega, sort = FALSE)[, idx := .I]
  
  # Simulate data
  synth_cnt <- synth_cat <- NULL
  if (!is.null(params$cnt)) {
    fam <- params$meta[family != 'multinom', unique(family)]
    psi <- merge(omega, params$cnt, by = 'f_idx', sort = FALSE, allow.cartesian = TRUE)
    if (!is.null(evidence) & any(evidence$operator %in% c('<', '<=', '>', '>='))) {
      for (k in evidence[, which(grepl('<', operator))]) {
        j <- evidence$variable[k]
        value <- as.numeric(evidence$value[k])
        psi[variable == j & max > value, max := value]
      }
      for (k in evidence[, which(grepl('>', operator))]) {
        j <- evidence$variable[k]
        value <- as.numeric(evidence$value[k])
        psi[variable == j & min < value, min := value]
      }
    }
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
  
  # Combine, optionally impose constraint(s)
  x_synth <- cbind(synth_cnt, synth_cat)
  if (!is.null(evidence) & any(grepl('==', evidence$operator))) {
    for (k in evidence[, which(operator == '==')]) {
      j <- evidence$variable[k]
      value <- evidence$value[k]
      if (params$meta[variable == j, class == 'numeric']) {
        value <- as.numeric(value)
      }
      x_synth[[j]] <- value
    }
  }
  
  # Clean up, export
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
      lapply(x_synth[, idx_logical, drop = FALSE], function(x) {x == 'TRUE'})
    )
  }
  if (sum(idx_integer) > 0L) {
    x_synth[, idx_integer] <- as.data.frame(
      lapply(x_synth[, idx_integer, drop = FALSE], function(x) as.integer(x)) 
    )
  }
  if ('data.table' %in% params$input_class) {
    x_synth <- as.data.table(x_synth)
  } else if ('tbl_df' %in% params$input_class) {
    x_synth <- tibble::as_tibble(x_synth)
  } else if ('matrix' %in% params$input_class) {
    x_synth <- as.matrix(x_synth)
  }
  return(x_synth)
}


