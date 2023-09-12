#' Likelihood Estimation
#' 
#' Compute the likelihood of input data, optionally conditioned on some event(s).
#' 
#' @param params Circuit parameters learned via \code{\link{forde}}. 
#' @param query Data frame of samples, optionally comprising just a subset of 
#'   training features. Likelihoods will be computed for each sample. Missing
#'   features will be marginalized out. See Details.
#' @param evidence Optional set of conditioning events. This can take one of 
#'   three forms: (1) a partial sample, i.e. a single row of data with some but
#'   not all columns; (2) a data frame of conditioning events, which allows for 
#'   inequalities; or (3) a posterior distribution over leaves. See Details.
#' @param arf Pre-trained \code{\link{adversarial_rf}} or other object of class 
#'   \code{ranger}. This is not required but speeds up computation considerably
#'   for total evidence queries. (Ignored for partial evidence queries.)
#' @param oob Only use out-of-bag leaves for likelihood estimation? If 
#'   \code{TRUE}, \code{x} must be the same dataset used to train \code{arf}.
#'   Only applicable for total evidence queries.
#' @param log Return likelihoods on log scale? Recommended to prevent underflow.
#' @param batch Batch size. The default is to compute densities for all of 
#'   queries in one round, which is always the fastest option if memory allows. 
#'   However, with large samples or many trees, it can be more memory efficient 
#'   to split the data into batches. This has no impact on results.
#' @param parallel Compute in parallel? Must register backend beforehand, e.g. 
#'   via \code{doParallel}.
#'   
#'   
#' @details 
#' This function computes the likelihood of input data, optionally conditioned
#' on some event(s). Queries may be partial, i.e. covering some but not all
#' features, in which case excluded variables will be marginalized out. 
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
#' A vector of likelihoods, optionally on the log scale. 
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
#' # Estimate average log-likelihood
#' arf <- adversarial_rf(iris)
#' psi <- forde(arf, iris)
#' ll <- lik(psi, iris, arf = arf, log = TRUE)
#' mean(ll)
#' 
#' # Identical but slower
#' ll <- lik(psi, iris, log = TRUE)
#' mean(ll)
#' 
#' # Partial evidence query
#' lik(psi, query = iris[1, 1:3])
#' 
#' # Condition on Species = "setosa"
#' evi <- data.frame(Species = "setosa")
#' lik(psi, query = iris[1, 1:3], evidence = evi)
#' 
#' # Condition on Species = "setosa" and Petal.Width > 0.3
#' evi <- data.frame(variable = c("Species", "Petal.Width"),
#'                   relation = c("==", ">"), 
#'                   value = c("setosa", 0.3))
#' lik(psi, query = iris[1, 1:3], evidence = evi)
#' 
#' 
#' @seealso
#' \code{\link{adversarial_rf}}, \code{\link{forge}}
#' 
#'
#' @export
#' @import ranger 
#' @import data.table
#' @importFrom stats predict
#' @importFrom foreach foreach %do% %dopar%
#' @importFrom truncnorm dtruncnorm 
#' 

lik <- function(
    params, 
    query,
    evidence = NULL,
    arf = NULL,
    oob = FALSE,
    log = TRUE, 
    batch = NULL, 
    parallel = TRUE) {
  
  # To avoid data.table check issues
  tree <- cvg <- leaf <- variable <- mu <- sigma <- value <- obs <- prob <- 
    V1 <- relation <- f_idx <- wt <- val <- family <- fold <- . <- NULL
  
  # Check query
  x <- as.data.frame(query)
  colnames_x <- colnames(x)
  n <- nrow(x)
  d <- ncol(x)
  if (d == params$meta[, .N] & is.null(arf)) {
    warning('For total evidence queries, it is faster to include the ', 
            'pre-trained arf.')
  }
  if (any(!colnames(x) %in% params$meta$variable)) {
    err <- setdiff(colnames(x), params$meta$variable)
    stop('Unrecognized feature(s) among colnames: ', err)
  }
  x <- suppressWarnings(prep_x(x))
  factor_cols <- sapply(x, is.factor)
  
  # Prep evidence
  conj <- FALSE
  if (!is.null(evidence)) {
    evidence <- prep_evi(params, evidence)
    if (!all(c('f_idx', 'wt') %in% colnames(evidence))) {
      conj <- TRUE
    }
  } 
  
  # Check ARF
  if (d == params$meta[, .N] & !is.null(arf)) {
    num_trees <- arf$num.trees
    preds <- stats::predict(arf, x, type = 'terminalNodes')$predictions + 1L
    preds <- data.table('tree' = rep(seq_len(num_trees), each = n), 
                        'leaf' = as.vector(preds),
                        'obs' = rep(seq_len(n), times = num_trees))
    if (isTRUE(oob)) {
      preds <- na.omit(preds)
    }
    preds <- merge(preds, params$forest[, .(tree, leaf, f_idx)], 
                   by = c('tree', 'leaf'), sort = FALSE)
    setnames(x, params$meta$variable)
  } else {
    arf <- NULL
    setnames(x, colnames_x)
  }
  
  # PMF over leaves
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
  leaves <- omega[, f_idx]
  
  # Optional batching
  if (is.null(batch)) {
    batch <- n
  }
  k <- round(n/batch)
  if (k < 1) {
    k <- 1L
  }
  batch_idx <- suppressWarnings(split(seq_len(n), seq_len(k)))
  
  # Likelihood function 
  lik_fn <- function(fold, arf) {
    
    # Prep work
    psi_cnt <- psi_cat <- NULL
    pure <- all(factor_cols) | all(!factor_cols)
    if (is.null(arf) & !isTRUE(pure)) {
      omega_tmp <- rbindlist(lapply(batch_idx[[fold]], function(i) {
        omega$obs <- i
        omega$wt <- NULL
        return(omega)
      })) 
    }
    
    # Continuous data
    if (any(!factor_cols)) {
      fam <- params$meta[class == 'numeric', unique(family)]
      x_long <- melt(
        data.table(obs = batch_idx[[fold]], 
                   x[batch_idx[[fold]], !factor_cols, drop = FALSE]), 
        id.vars = 'obs', variable.factor = FALSE
      )
      if (is.null(arf)) {
        psi_cnt <- merge(params$cnt[f_idx %in% leaves], x_long, by = 'variable', 
                         sort = FALSE, allow.cartesian = TRUE)
        rm(x_long)
      } else {
        preds_cnt <- merge(preds[f_idx %in% leaves], x_long, by = 'obs', 
                           sort = FALSE, allow.cartesian = TRUE)
        rm(x_long)
        psi_cnt <- merge(params$cnt[f_idx %in% leaves], preds_cnt, 
                         by = c('f_idx', 'variable'), sort = FALSE)
        rm(preds_cnt)
      }
      if (fam == 'truncnorm') {
        psi_cnt[, lik := truncnorm::dtruncnorm(value, a = min, b = max, 
                                               mean = mu, sd = sigma)]
      } else if (fam == 'unif') {
        psi_cnt[, lik := stats::dunif(value, min = min, max = max)]
      }
      psi_cnt[value == min, lik := 0]
      psi_cnt[, lik := prod(lik), by = .(f_idx, obs)]
      psi_cnt <- unique(psi_cnt[lik > 0, .(f_idx, obs, lik)])
      if (is.null(arf) & !isTRUE(pure)) {
        omega_tmp <- merge(omega_tmp, psi_cnt[, .(f_idx, obs)], 
                           by = c('f_idx', 'obs'), sort = FALSE)
        leaves <- omega_tmp[, unique(f_idx)]
      }
    }
    
    # Categorical data
    if (any(factor_cols)) {
      x_tmp <- x[batch_idx[[fold]], factor_cols, drop = FALSE]
      n_tmp <- nrow(x_tmp)
      x_long <- melt(
        data.table(obs = batch_idx[[fold]], x_tmp), 
        id.vars = 'obs', value.name = 'val', variable.factor = FALSE
      )
      # Speedups are possible if there are many duplicates
      is_unique <- !duplicated(x_tmp)
      if (all(is_unique)) {
        x_unique <- x_long
        colnames(x_unique)[1] <- 's_idx'
      } else {
        x_unique <- unique(x_tmp)
        x_unique <- melt(
          data.table(s_idx = seq_len(nrow(x_unique)), x_unique),
          id.vars = 's_idx', value.name = 'val', variable.factor = FALSE
        )
        s_idx <- integer(length = n_tmp)
        s_idx[is_unique] <- seq_len(sum(is_unique))
        for (i in 2:n_tmp) {
          if (s_idx[i] == 0L) {
            s_idx[i] <- s_idx[i - 1L]
          }
        }
        idx_dt <- data.table(obs = batch_idx[[fold]], s_idx = s_idx)
      }
      if (is.null(arf)) {
        grd <- rbindlist(lapply(which(factor_cols), function(j) {
          expand.grid('f_idx' = leaves, 'variable' = colnames(x)[j],
                      'val' = x_long[variable == colnames(x)[j], unique(val)],
                      stringsAsFactors = FALSE)
        }))
        rm(x_long)
        psi_cat <- merge(params$cat[f_idx %in% leaves], grd, 
                         by = c('f_idx', 'variable', 'val'), 
                         sort = FALSE, all.y = TRUE)
        rm(grd)
        psi_cat[is.na(prob), prob := 0]
        psi_cat <- merge(psi_cat, x_unique, by = c('variable', 'val'), 
                         sort = FALSE, allow.cartesian = TRUE)
        psi_cat[, lik := prod(prob), by = .(f_idx, s_idx)]
        psi_cat <- unique(psi_cat[lik > 0, .(f_idx, s_idx, lik)])
        if (all(is_unique)) {
          setnames(psi_cat, 's_idx', 'obs')
        } else {
          if (!isTRUE(pure)) {
            omega_tmp <- merge(idx_dt, omega_tmp, by = 'obs', sort = FALSE)
            psi_cat <- merge(psi_cat, omega_tmp, by = c('f_idx', 's_idx'),
                             sort = FALSE)[, s_idx := NULL]
            rm(omega_tmp)
            setcolorder(psi_cat, c('f_idx', 'obs', 'lik'))
            psi_cnt <- merge(psi_cnt, psi_cat[, .(f_idx, obs)], 
                             by = c('f_idx', 'obs'), sort = FALSE)
          }
        }
      } else {
        preds_cat <- merge(preds[f_idx %in% leaves], x_long, by = 'obs', 
                           sort = FALSE, allow.cartesian = TRUE)
        rm(x_long)
        psi_cat <- merge(params$cat, preds_cat, by = c('f_idx', 'variable', 'val'),
                         sort = FALSE, allow.cartesian = TRUE, all.y = TRUE)
        rm(preds_cat)
        psi_cat[is.na(prob), prob := 0]
        psi_cat[, lik := prod(prob), by = .(f_idx, obs)]
        psi_cat <- unique(psi_cat[lik > 0, .(f_idx, obs, lik)])
      }
    }
    
    # Put it together
    psi_x <- rbind(psi_cnt, psi_cat)
    if (!isTRUE(pure)) {
      psi_x <- psi_x[, prod(lik), by = .(f_idx, obs)]
      setnames(psi_x, 'V1', 'lik')
    }
    return(psi_x)
  }
  if (isTRUE(parallel)) {
    out <- foreach(fold = seq_len(k), .combine = rbind) %dopar% lik_fn(fold, arf)
  } else {
    out <- foreach(fold = seq_len(k), .combine = rbind) %do% lik_fn(fold, arf)
  }
  
  # Compute per-sample likelihoods
  out <- merge(out, omega, by = 'f_idx', sort = FALSE)
  out <- out[, log(crossprod(wt, lik)), by = obs]
  setnames(out, 'V1', 'lik')
  
  # Anybody missing?
  zeros <- setdiff(seq_len(n), out[, obs])
  if (length(zeros) > 0L) {
    zero_dt <- data.table(obs = zeros, lik = -Inf)
    out <- rbind(out, zero_dt)
  }
  
  # Export
  if (!isTRUE(log)) {
    out[, lik := exp(lik)]
  }
  return(out[order(obs), lik])
}


