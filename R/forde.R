#' Forests for density estimation
#' 
#' Uses a pre-trained ARF model to estimate leaf and distribution parameters.
#' 
#' @param rf Pre-trained adversarial random forest. Alternatively, any object
#'   of class \code{ranger}.
#' @param x_trn Training data for estimating parameters.
#' @param x_tst Optional test data. If supplied, the function computes 
#'   log-likelihoods on test data (measured in nats).
#' @param family Distribution to use for density estimation of continuous 
#'   features. Current options include truncated normal (the default
#'   \code{family = "truncnorm"}) and uniform (\code{family = "unif"}). See 
#'   Details.
#' @param epsilon Slack parameter on empirical bounds when \code{family = "unif"}.
#'   This avoids zero-density points when test data fall outside the support
#'   of training data. The gap between lower and upper bounds is expanded by 
#'   a factor of \code{epsilon}. Only used when a variable is never selected for
#'   splitting.
#' @param loglik Return log-likelihood of training or test data? If \code{FALSE},
#'   function only returns leaf and distribution parameters \code{psi}. This
#'   can save time when the ultimate goal is data synthesis.
#' @param batch Batch size. The default is to compute parameters for the full 
#'   dataset in one round, which is always the fastest option if memory allows. 
#'   However, with large samples or many trees, it can be more memory efficient 
#'   to split the data into batches. This has no impact on results.
#' @param parallel Compute in parallel? Must register backend beforehand, e.g. 
#'   via \code{doParallel}.
#'   
#'   
#' @details 
#' \code{forde} extracts leaf parameters from a pretrained forest and learns
#' distribution parameters for data within each leaf. The former includes 
#' coverage (proportion of data falling into the leaf) and split criteria. The 
#' latter includes proportions for categorical features and mean/variance for
#' continuous features. These values are stored in a \code{data.table} called
#' \code{psi}, which can be passed onto \code{\link{forge}} for data synthesis. 
#' They are also used to compute log-likelihoods for \code{x_trn} or \code{x_tst} 
#' if \code{loglik = TRUE}.
#' 
#' Currently, \code{forde} only provides support for truncated normal or uniform
#' densities when features are continuous. Future releases will accommodate 
#' a larger class of distributional families.
#' 
#' Though \code{forde} was designed to take an adversarial random forest as 
#' input, \code{rf} can in principle be any object of class \code{ranger}. This
#' allows users to test performance with alternative pipelines (e.g., with 
#' supervised forest input). There is also no requirement that \code{x_trn} be 
#' the data used to fit \code{rf}.
#' 
#' 
#' @return 
#' A list with two elements: 
#' \itemize{
#'   \item \code{psi}, a \code{data.table} with leaf and distribution parameters
#'   for each feature. These are used for computing log-likelihoods or 
#'   generating data with \code{\link{forge}}.
#'   \item \code{loglik}, a vector of log-likelihoods, one for each sample in 
#'   \code{x_trn} or \code{x_tst}, if the latter is provided. To skip these 
#'   calculations, set \code{loglik = FALSE}, in which case this item is 
#'   \code{NULL}.
#' }
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
#' head(fd$psi)
#' head(fd$loglik)
#' 
#' 
#' @seealso
#' \code{\link{adversarial_rf}}, \code{\link{forge}}
#' 
#'
#' @export
#' @import ranger 
#' @import data.table
#' @importFrom foreach foreach %do% %dopar%
#' @importFrom truncnorm dtruncnorm 
#' 

forde <- function(
    rf, 
    x_trn, 
    x_tst = NULL, 
    family = 'truncnorm', 
    epsilon = 0.1, 
    loglik = TRUE, 
    batch = NULL, 
    parallel = TRUE) {

  # Prelimz
  if (!is.null(x_tst)) {
    if (ncol(x_tst) != ncol(x_trn)) {
      stop('x_trn and x_tst must be the same dimensionality.')
    }
    if (colnames(x_tst) != colnames(x_trn)) {
      stop('x_trn and x_tst must have identical colnames.')
    }
  }
  if (!family %in% c('truncnorm', 'unif')) {
    stop('family not recognized.')
  }
  
  # Prep data
  x <- as.data.frame(x_trn)
  n <- nrow(x)
  d <- ncol(x)
  idx_char <- sapply(x, is.character)
  if (any(idx_char)) {
    x[, idx_char] <- as.data.frame(
      lapply(x[, idx_char, drop = FALSE], as.factor)
    )
  }
  idx_logical <- sapply(x, is.logical)
  if (any(idx_logical)) {
    x[, idx_logical] <- as.data.frame(
      lapply(x[, idx_logical, drop = FALSE], as.factor)
    )
  }
  factor_cols <- sapply(x, is.factor)
  
  # Compute leaf bounds and coverage
  num_trees <- rf$num.trees
  pred <- predict(rf, x, type = 'terminalNodes')$predictions + 1
  bnd_fn <- function(tree) {
    num_nodes <- length(rf$forest$split.varIDs[[tree]])
    lb <- matrix(-Inf, nrow = num_nodes, ncol = d)
    ub <- matrix(Inf, nrow = num_nodes, ncol = d)
    if (family == 'unif') {
      for (j in seq_len(d)) {
        if (!isTRUE(factor_cols[j])) {
          gap <- max(x[[j]]) - min(x[[j]])
          lb[, j] <- min(x[[j]]) - epsilon/2 * gap
          ub[, j] <- max(x[[j]]) + epsilon/2 * gap
        }
      }
    }
    for (i in 1:num_nodes) {
      left_child <- rf$forest$child.nodeIDs[[tree]][[1]][i] + 1
      right_child <- rf$forest$child.nodeIDs[[tree]][[2]][i] + 1
      splitvarID <- rf$forest$split.varIDs[[tree]][i] + 1
      splitval <- rf$forest$split.value[[tree]][i]
      if (left_child > 1 & left_child != right_child) {
        ub[left_child, ] <- ub[right_child, ] <- ub[i, ]
        lb[left_child, ] <- lb[right_child, ] <- lb[i, ]
        ub[left_child, splitvarID] <- lb[right_child, splitvarID] <- splitval
      }
    }
    leaves <- which(rf$forest$child.nodeIDs[[tree]][[1]] == 0) 
    colnames(lb) <- rf$forest$independent.variable.names
    colnames(ub) <- rf$forest$independent.variable.names
    merge(melt(data.table(tree = tree, leaf = leaves, lb[leaves, ]), 
               id.vars = c('tree', 'leaf'), value.name = 'min'), 
          melt(data.table(tree = tree, leaf = leaves, ub[leaves, ]), 
               id.vars = c('tree', 'leaf'), value.name = 'max'), 
          by = c('tree', 'leaf', 'variable'))
  }
  if (isTRUE(parallel)) {
    bnds <- foreach(tree = 1:num_trees, .combine = rbind) %dopar% bnd_fn(tree)
  } else {
    bnds <- foreach(tree = 1:num_trees, .combine = rbind) %do% bnd_fn(tree)
  }
  bnds[, cvg := sum(pred[, tree] == leaf) / n, by = .(tree, leaf)]
  # Can't do anything with coverage 1/n
  if (any(!factor_cols)) {
    bnds[cvg == 1/n, cvg := 0]
  }
  
  # Calculate distribution parameters
  psi_cnt <- psi_cat <- NULL
  # Continuous case
  if (any(!factor_cols)) {
    psi_cnt_fn <- function(tree) {
      dt <- data.table(tree = tree, x[, !factor_cols, drop = FALSE], leaf = pred[, tree])
      long <- melt(dt, id.vars = c('tree', 'leaf'))
      if (family == 'truncnorm') {
        long[, list(cat = NA_character_, prob = NA_real_, 
                    mu = mean(value), sigma = sd(value), type = 'cnt'), 
             by = .(tree, leaf, variable)]
      } else if (family == 'unif') {
        long[, list(cat = NA_character_, prob = NA_real_, 
                    mu = NA_real_, sigma = NA_real_, type = 'cnt'), 
             by = .(tree, leaf, variable)]
      }
    }
    if (isTRUE(parallel)) {
      psi_cnt <- foreach(tree = 1:num_trees, .combine = rbind) %dopar% psi_cnt_fn(tree)
    } else {
      psi_cnt <- foreach(tree = 1:num_trees, .combine = rbind) %do% psi_cnt_fn(tree)
    }
  } 
  # Categorical case
  if (any(factor_cols)) {
    psi_cat_fn <- function(tree) {
      dt <- data.table(tree = tree, x[, factor_cols, drop = FALSE], leaf = pred[, tree])
      long <- melt(dt, id.vars = c('tree', 'leaf'), value.factor = FALSE, value.name = 'cat')
      long[, count := .N, by = .(tree, leaf, variable)]
      unique(setDT(long)[, list(prob = .N/count, mu = NA_real_, sigma = NA_real_, type = 'cat'), 
                         by = .(tree, leaf, variable, cat)])
    }
    if (isTRUE(parallel)) {
      psi_cat <- foreach(tree = 1:num_trees, .combine = rbind) %dopar% psi_cat_fn(tree)
    } else {
      psi_cat <- foreach(tree = 1:num_trees, .combine = rbind) %do% psi_cat_fn(tree)
    }
  } 
  psi_tmp <- rbind(psi_cnt, psi_cat)
  psi <- merge(psi_tmp, bnds, by = c('tree', 'leaf', 'variable'))
  rm(psi_cnt, psi_cat, psi_tmp)
  
  # Log-likelihood calculation
  if (isTRUE(loglik)) {
    # Optionally prep test data
    if (!is.null(x_tst)) {
      x <- as.data.frame(x_tst)
      n <- nrow(x_tst)
      idx_char <- sapply(x, is.character)
      if (any(idx_char)) {
        x[, idx_char] <- as.data.frame(
          lapply(x[, idx_char, drop = FALSE], as.factor)
        )
      }
      idx_logical <- sapply(x, is.logical)
      if (any(idx_logical)) {
        x[, idx_logical] <- as.data.frame(
          lapply(x[, idx_logical, drop = FALSE], as.factor)
        )
      }
      factor_cols <- sapply(x, is.factor)
      pred <- predict(rf, x, type = 'terminalNodes')$predictions + 1
    }
    
    # Optional batch index
    if (!is.null(batch)) {
      k <- round(n / batch)
      batch_idx <- suppressWarnings(split(1:n, seq_len(k)))
    } else {
      k <- 1L
      batch_idx <- list(1:n)
    }
    # Compute per-feature likelihoods
    loglik_fn <- function(fold) {
      psi_x_cnt <- psi_x_cat <- NULL
      # Predictions
      preds <- rbindlist(lapply(1:ncol(pred), function(b) {
        data.table(tree = b, leaf = pred[batch_idx[[fold]], b], obs = batch_idx[[fold]])
      }))
      # Continuous data
      if (any(!factor_cols)) {
        x_long_cnt <- melt(
          data.table(obs = batch_idx[[fold]], 
                     x[batch_idx[[fold]], !factor_cols, drop = FALSE]), 
          id.vars = 'obs'
        )
        preds_x_cnt <- merge(preds, x_long_cnt, by = 'obs', allow.cartesian = TRUE)
        psi_x_cnt <- merge(psi[type == 'cnt', .(tree, leaf, cvg, variable, min, max, mu, sigma)], 
                           preds_x_cnt, by = c('tree', 'leaf', 'variable'))
        if (family == 'truncnorm') {
          psi_x_cnt[, lik := dtruncnorm(value, a = min, b = max, mean = mu, sd = sigma)]
        } else if (family == 'unif') {
          psi_x_cnt[, lik := dunif(value, min = min, max = max)]
        } 
        psi_x_cnt <- psi_x_cnt[, .(tree, obs, cvg, lik)]
        rm(x_long_cnt, preds_x_cnt)
      }
      # Categorical data
      if (any(factor_cols)) {
        x_long_cat <- melt(
          data.table(obs = batch_idx[[fold]], 
                     x[batch_idx[[fold]], factor_cols, drop = FALSE]), 
          id.vars = 'obs', value.name = 'cat'
        )
        preds_x_cat <- merge(preds, x_long_cat, by = 'obs', allow.cartesian = TRUE)
        psi_x_cat <- merge(psi[type == 'cat', .(tree, leaf, cvg, variable, cat, prob)], 
                           preds_x_cat, by = c('tree', 'leaf', 'variable', 'cat'), 
                           allow.cartesian = TRUE)
        psi_x_cat[, lik := prob]
        psi_x_cat <- psi_x_cat[, .(tree, obs, cvg, lik)]
        rm(x_long_cat, preds_x_cat)
      } 
      rm(preds)
      # Put it together
      psi_x <- rbind(psi_x_cnt, psi_x_cat)
      rm(psi_x_cnt, psi_x_cat)
      # Compute per-sample log-likelihoods
      loglik <- unique(psi_x[, prod(lik) * cvg, by = .(obs, tree)])
      loglik[is.na(V1), V1 := 0]
      loglik <- loglik[, log(mean(V1)), by = obs]
      return(loglik)
    }
    if (k == 1L) {
      ll <- loglik_fn(1)
    } else {
      if (isTRUE(parallel)) {
        ll <- foreach(fold = 1:k, .combine = rbind) %dopar% loglik_fn(fold)
      } else {
        ll <- foreach(fold = 1:k, .combine = rbind) %do% loglik_fn(fold)
      }
    }
    loglik <- ll[order(obs), V1]
  } else {
    loglik <- NULL
  }
  
  # Export
  out <- list('psi' = psi, 'loglik' = loglik)
  return(out)
}
