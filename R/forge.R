#' Forests for Generative Modeling
#' 
#' Uses pre-trained FORDE model to simulate synthetic data.
#' 
#' @param pc Probabilistic circuit learned via \code{\link{forde}}. 
#' @param n_synth Number of synthetic samples to generate.
#' @param evidence Optional set of conditioning events. This can take one of 
#'   three forms: (1) a partial sample, i.e. a single row of data with
#'   some but not all columns; (2) a data frame of conditioning events, 
#'   which allows for inequalities; or (3) a posterior distribution over leaves.
#'   See Details.
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
#' There are three methods for (optionally) encoding conditioning events via the 
#' \code{evidence} argument. The first is to provide a partial sample, where
#' some but not all columns from the training data are present. The second is to 
#' provide a data frame with three columns: \code{variable}, \code{operator}, 
#' and \code{value}. This supports inequalities via \code{operator}. 
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
#'                   operator = c("==", ">"), 
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
#' @importFrom tibble as_tibble
#' 

forge <- function(
    pc, 
    n_synth, 
    evidence = NULL,
    parallel = TRUE) {
  
  # To avoid data.table check issues
  tree <- cvg <- leaf <- idx <- family <- mu <- sigma <- prob <- dat <- 
    variable <- j <- f_idx <- val <- . <- NULL
  
  # Check evidence
  if (!is.null(evidence)) {
    evidence <- as.data.table(evidence)
    part <- all(colnames(evidence) %in% pc$meta$variable)
    conj <- all(c('variable', 'operator', 'value') %in% colnames(evidence))
    post <- all(c('f_idx', 'wt') %in% colnames(evidence))
    if (part + conj + post != 1L) {
      stop('evidence must either be a partial sample, a data frame of conjuncts, ', 
           'or a posterior distribution over leaves.')
    }
    if (part) {
      if (!all(colnames(evidence) %in% pc$meta$variable)) {
        err <- setdiff(colnames(evidence), pc$meta$variable)
        stop('Unrecognized feature(s) among colnames: ', err)
      }
      evidence <- melt(evidence, measure.vars = colnames(evidence))
      evidence[, operator := '==']
      conj <- TRUE
    }
    if (conj) {
      if (max(evidence[, .N, by = variable]$N > 1L)) {
        stop('Only one constraint per variable allowed when using conjuncts.')
      }
      evi <- merge(pc$meta, evidence, by = 'variable', sort = FALSE)
      if (evi[class == 'numeric' & operator == '!=', .N] > 0) {
        evidence <- evidence[!(class == 'numeric' & operator == '!=')]
        warning('With continuous features, "!=" is not a valid operator. ', 
                'This constraint has been removed.')
      }
      if (evi[class != 'numeric' & !operator %in% c('==', '!='), .N] > 0) {
        stop('With categorical features, the only valid operators are ',
             '"==" or "!=".')
      }
    }
  } else {
    conj <- FALSE
  }
  
  # Prepare the event space
  if (is.null(evidence)) {
    omega <- pc$forest
    omega[, wt := cvg / max(tree)]
    omega <- omega[, .(f_idx, wt)]
  } else if (conj) {
    omega <- leaf_posterior(pc, evidence, parallel)
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
  if (!is.null(pc$cnt)) {
    fam <- pc$meta[family != 'multinom', unique(family)]
    psi <- merge(omega, pc$cnt, by = 'f_idx', sort = FALSE, allow.cartesian = TRUE)
    if (conj) {
      if (any(evidence$operator %in% c('<', '<=', '>', '>='))) {
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
    }
    if (fam == 'truncnorm') {
      psi[, dat := truncnorm::rtruncnorm(.N, a = min, b = max, mean = mu, sd = sigma)]
    } else if (fam == 'unif') {
      psi[, dat := stats::runif(.N, min = min, max = max)]
    }
    synth_cnt <- dcast(psi, idx ~ variable, value.var = 'dat')[, idx := NULL]
  }
  if (!is.null(pc$cat)) {
    sim_cat <- function(j) {
      psi <- merge(omega, pc$cat[variable == j], by = 'f_idx', sort = FALSE, 
                   allow.cartesian = TRUE)
      psi[prob == 1, dat := val]
      if (conj) {
        if (evidence[variable == j & operator == '!=', .N] == 1L) {
          value <- evidence[variable == j, value]
          psi <- psi[val != value]
        }
      }
      psi[prob < 1, dat := sample(val, 1, prob = prob), by = idx] 
      setnames(data.table(unique(psi[, .(idx, dat)])[, dat]), j)
    }
    cat_vars <- pc$meta[family == 'multinom', variable]
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
      if (pc$meta[variable == j, class == 'numeric']) {
        value <- as.numeric(value)
      }
      x_synth[[j]] <- value
    }
  }
  
  # Clean up, export
  setcolorder(x_synth, pc$meta$variable)
  setDF(x_synth)
  idx_factor <- pc$meta[, which(class == 'factor')]
  idx_logical <- pc$meta[, which(class == 'logical')]
  idx_integer <- pc$meta[, which(class == 'integer')]
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
  if ('data.table' %in% pc$input_class) {
    x_synth <- as.data.table(x_synth)
  } else if ('tbl_df' %in% pc$input_class) {
    x_synth <- tibble::as_tibble(x_synth)
  } else if ('matrix' %in% pc$input_class) {
    x_synth <- as.matrix(x_synth)
  }
  return(x_synth)
}


