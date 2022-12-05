#' Likelihood Estimation
#' 
#' Compute the density of input data.
#' 
#' @param arf Pre-trained adversarial random forest. Alternatively, any object
#'   of class \code{ranger}.
#' @param params Parameters learned via \code{forde}. 
#' @param x Input data. Function outputs densities for each sample.
#' @param oob Only use out-of-bag leaves for likelihood estimation? If 
#'   \code{TRUE}, \code{x} must be the same dataset used to train \code{rf}.
#' @param log Return likelihoods on log scale?
#' @param batch Batch size. The default is to compute densities for all of 
#'   \code{x} in one round, which is always the fastest option if memory allows. 
#'   However, with large samples or many trees, it can be more memory efficient 
#'   to split the data into batches. This has no impact on results.
#' @param parallel Compute in parallel? Must register backend beforehand, e.g. 
#'   via \code{doParallel}.
#'   
#'   
#' @details 
#' This function computes the density of input data according to a FORDE model
#' using a pre-trained ARF. Each sample's likelihood is a weighted average of 
#' its likelihood in all leaves whose split criteria it satisfies. Intra-leaf
#' densities are fully factorized, since ARFs satisfy the local independence
#' criterion by construction.
#' 
#' 
#' @return 
#' A vector of likelihoods, optionally on the log scale. 
#' 
#' 
#' @references 
#' Watson, D., Blesch, K., Kapar, J., & Wright, M. (2022). Adversarial random 
#' forests for density estimation and generative modeling. \emph{arXiv} preprint,
#' 2205.09435.
#' 
#' 
#' @examples
#' # Estimate average log-likelihood
#' arf <- adversarial_rf(iris)
#' psi <- forde(arf, iris)
#' ll <- lik(arf, psi, iris, log = TRUE)
#' mean(ll)
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

lik <- function(
    arf, 
    params, 
    x, 
    oob = FALSE, 
    log = TRUE, 
    batch = NULL, 
    parallel = TRUE) {
  
  # To avoid data.table check issues
  tree <- cvg <- leaf <- variable <- mu <- sigma <- value <- obs <- prob <- V1 <- loglik <- family <- fold <- . <- NULL
  
  # Prelimz
  x <- as.data.frame(x)
  n <- nrow(x)
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
  idx_intgr <- sapply(x, is.integer)
  if (any(idx_intgr)) {
    for (j in which(idx_intgr)) {
      lvls <- sort(unique(x[, j]))
      x[, j] <- factor(x[, j], levels = lvls, ordered = TRUE)
    }
  }
  factor_cols <- sapply(x, is.factor)
  pred <- predict(arf, x, type = 'terminalNodes')$predictions + 1L
  
  # Optional batch index
  if (is.null(batch)) {
    batch <- n
  }
  k <- round(n/batch)
  if (k < 1) {
    k <- 1
  }
  batch_idx <- suppressWarnings(split(1:n, seq_len(k)))
  
  # Likelihood function
  num_trees <- arf$num.trees
  lik_fn <- function(fold) {
    params_x_cnt <- params_x_cat <- NULL
    # Predictions
    preds <- rbindlist(lapply(1:num_trees, function(b) {
      data.table(tree = b, leaf = pred[batch_idx[[fold]], b], obs = batch_idx[[fold]])
    }))
    if (isTRUE(oob)) {
      preds <- preds[!is.na(leaf)]
    }
    # Continuous data
    if (any(!factor_cols)) {
      x_long_cnt <- melt(
        data.table(obs = batch_idx[[fold]], 
                   x[batch_idx[[fold]], !factor_cols, drop = FALSE]), 
        id.vars = 'obs', variable.factor = FALSE
      )
      preds_x_cnt <- merge(preds, x_long_cnt, by = 'obs', sort = FALSE, allow.cartesian = TRUE)
      params_x_cnt <- merge(params[family != 'multinom', .(tree, leaf, cvg, variable, min, max, mu, sigma)], 
                            preds_x_cnt, by = c('tree', 'leaf', 'variable'), sort = FALSE)
      family <- unique(params[family != "multinom", family])
      if (family == 'truncnorm') {
        params_x_cnt[, lik := truncnorm::dtruncnorm(value, a = min, b = max, mean = mu, sd = sigma)]
      } else if (family == 'unif') {
        params_x_cnt[, lik := stats::dunif(value, min = min, max = max)]
      } 
      params_x_cnt <- params_x_cnt[, .(tree, obs, cvg, lik)]
      rm(x_long_cnt, preds_x_cnt)
    }
    # Categorical data
    if (any(factor_cols)) {
      x_long_cat <- melt(
        data.table(obs = batch_idx[[fold]], 
                   x[batch_idx[[fold]], factor_cols, drop = FALSE]), 
        id.vars = 'obs', value.name = 'cat', variable.factor = FALSE
      )
      preds_x_cat <- merge(preds, x_long_cat, by = 'obs', sort = FALSE, allow.cartesian = TRUE)
      params_x_cat <- merge(params[family == 'multinom', .(tree, leaf, cvg, variable, cat, prob)], 
                            preds_x_cat, by = c('tree', 'leaf', 'variable', 'cat'), 
                            sort = FALSE, allow.cartesian = TRUE)
      params_x_cat[, lik := prob]
      params_x_cat <- params_x_cat[, .(tree, obs, cvg, lik)]
      rm(x_long_cat, preds_x_cat)
    } 
    rm(preds)
    # Put it together
    params_x <- rbind(params_x_cnt, params_x_cat)
    rm(params_x_cnt, params_x_cat)
    # Compute per-sample log-likelihoods
    lik <- unique(params_x[, prod(lik) * cvg, by = .(obs, tree)])
    lik[is.na(V1), V1 := 0]
    if (isTRUE(log)) {
      out <- loglik[, log(mean(V1)), by = obs]
    } else {
      out <- loglik[, mean(V1), by = obs]
    }
    return(out)
  }
  if (k == 1L) {
    ll <- lik_fn(1)
  } else {
    if (isTRUE(parallel)) {
      ll <- foreach(fold = 1:k, .combine = rbind) %dopar% lik_fn(fold)
    } else {
      ll <- foreach(fold = 1:k, .combine = rbind) %do% lik_fn(fold)
    }
  }
  return(ll[order(obs), V1])
}


