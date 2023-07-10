#' Likelihood Estimation
#' 
#' Compute the likelihood of input data, optionally conditioned on some event(s).
#' 
#' @param pc Probabilistic circuit learned via \code{\link{forde}}. 
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
#' evi <- data.frame(variable = c("Species", "Sepal.Length"),
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
#' @importFrom matrixStats logSumExp
#' 

lik <- function(
    pc, 
    query,
    evidence = NULL,
    arf = NULL,
    oob = FALSE,
    log = TRUE, 
    batch = NULL, 
    parallel = TRUE) {
  
  # To avoid data.table check issues
  tree <- cvg <- leaf <- variable <- mu <- sigma <- value <- obs <- prob <- 
    V1 <- family <- fold <- . <- NULL
  
  # Check query
  x <- as.data.frame(query)
  n <- nrow(x)
  d <- ncol(x)
  if (any(!colnames(x) %in% pc$meta$variable)) {
    err <- setdiff(colnames(x), pc$meta$variable)
    stop('Unrecognized feature(s) among colnames: ', err)
  }
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
      if (any(evidence$variable %in% colnames(x))) {
        warning('query and evidence features overlap. Likelihoods for these, ',
                'features will be set to one.')
      }
    }
  } else {
    conj <- FALSE
  }
  
  # Check ARF
  if (d == pc$meta[, .N] & !is.null(arf)) {
    if ('y' %in% colnames(x)) {
      colnames(x)[which(colnames(x) == 'y')] <- col_rename(x, 'y')
    }
    if ('obs' %in% colnames(x)) {
      colnames(x)[which(colnames(x) == 'obs')] <- col_rename(x, 'obs')
    }
    if ('tree' %in% colnames(x)) {
      colnames(x)[which(colnames(x) == 'tree')] <- col_rename(x, 'tree')
    } 
    if ('leaf' %in% colnames(x)) {
      colnames(x)[which(colnames(x) == 'leaf')] <- col_rename(x, 'leaf')
    } 
    preds <- stats::predict(arf, x, type = 'terminalNodes')$predictions + 1L
    preds <- rbindlist(lapply(seq_len(pc$forest[, max(tree)]), function(b) {
      data.table(tree = b, leaf = preds[, b], obs = seq_len(n))
    }))
    if (isTRUE(oob)) {
      preds <- preds[!is.na(leaf)]
    }
    preds <- merge(preds, pc$forest[, .(tree, leaf, f_idx)], 
                   by = c('tree', 'leaf'), sort = FALSE)
    colnames(x) <- pc$meta$variable
  } else {
    arf <- NULL
  }
  
  # PMF over leaves
  if (is.null(evidence)) {
    omega <- pc$forest
    omega[, wt := cvg / max(tree)]
    omega <- omega[, .(f_idx, wt)]
  } else if (conj) {
    omega <- leaf_posterior(pc, evidence, parallel)
  } else {
    omega <- evidence
  }
  omega <- omega[wt > 0]
  
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
    psi_cnt <- psi_cat <- NULL
    if (!is.null(arf)) {
      leaves <- omega[, f_idx]
    } else {
      omega_tmp <- rbindlist(lapply(batch_idx[[fold]], function(i) {
        omega$obs <- i
        return(omega)
      })) 
      leaves <- omega_tmp[, unique(f_idx)]
    }
    # Continuous data
    if (!is.null(pc$cnt)) {
      fam <- pc$meta[class == 'numeric', unique(family)]
      x_long <- melt(
        data.table(obs = batch_idx[[fold]], 
                   x[batch_idx[[fold]], !factor_cols, drop = FALSE]), 
        id.vars = 'obs', variable.factor = FALSE
      )
      if (is.null(arf)) {
        psi_cnt <- merge(pc$cnt[f_idx %in% leaves], x_long, by = 'variable', 
                         sort = FALSE, allow.cartesian = TRUE)
      } else {
        preds_cnt <- merge(preds[f_idx %in% leaves], x_long, by = 'obs', 
                           sort = FALSE, allow.cartesian = TRUE)
        psi_cnt <- merge(pc$cnt[f_idx %in% leaves], preds_cnt, 
                         by = c('f_idx', 'variable'), sort = FALSE)
        rm(preds_cnt)
      }
      rm(x_long)
      if (fam == 'truncnorm') {
        psi_cnt[, lik := log(truncnorm::dtruncnorm(value, a = min, b = max, 
                                                   mean = mu, sd = sigma))]
      } else if (fam == 'unif') {
        psi_cnt[, lik := stats::dunif(value, min = min, max = max, log = TRUE)]
      }
      psi_cnt[value == min, lik := -Inf]
      psi_cnt[, lik := sum(lik), by = .(obs, f_idx)]
      psi_cnt <- unique(psi_cnt[, .(f_idx, obs, lik)])
      if (is.null(arf)) {
        psi_cnt <- psi_cnt[is.finite(lik)]
        omega_tmp <- merge(omega_tmp, psi_cnt[, .(f_idx, obs)], 
                           by = c('f_idx', 'obs'), sort = FALSE)
        leaves <- omega_tmp[, unique(f_idx)]
      }
    }
    # Categorical data
    if (!is.null(pc$cat)) {
      if (is.null(arf)) {
        psi_cat <- rbindlist(lapply(which(factor_cols), function(j) {
          psi_j <- pc$cat[variable == colnames(x)[j] & f_idx %in% leaves]
          grd <- expand.grid(f_idx = leaves, val = psi_j[, unique(val)])
          psi_j <- merge(psi_j, grd, by = c('f_idx', 'val'), all.y = TRUE, sort = FALSE)
          rm(grd)
          psi_j[is.na(prob), prob := 0][, variable := NULL][, lik := log(prob)]
          x_long <- melt(
            data.table(obs = batch_idx[[fold]], 
                       x[batch_idx[[fold]], j, drop = FALSE]), 
            id.vars = 'obs', value.name = 'val', variable.factor = FALSE
          )
          tmp <- merge(x_long, omega_tmp[, .(f_idx, obs)], by = 'obs', sort = FALSE)
          psi_j <- merge(psi_j, tmp, by = c('f_idx', 'val'), sort = FALSE)
          return(psi_j[, .(f_idx, obs, lik)])
        }))
        psi_cat[, lik := sum(lik), by = .(obs, f_idx)]
        psi_cat <- unique(psi_cat[is.finite(lik), .(f_idx, obs, lik)])
        omega_tmp <- merge(omega_tmp, psi_cat[, .(f_idx, obs)], 
                           by = c('f_idx', 'obs'), sort = FALSE)
        if (!is.null(psi_cnt)) {
          psi_cnt <- merge(psi_cnt, omega_tmp[, .(f_idx, obs)], 
                           by = c('f_idx', 'obs'), sort = FALSE)
        }
      } else {
        x_long <- melt(
          data.table(obs = batch_idx[[fold]], 
                     x[batch_idx[[fold]], factor_cols, drop = FALSE]), 
          id.vars = 'obs', value.name = 'val', variable.factor = FALSE
        )
        preds_cat <- merge(preds[f_idx %in% leaves], x_long, by = 'obs', 
                           sort = FALSE, allow.cartesian = TRUE)
        rm(x_long)
        psi_cat <- merge(pc$cat, preds_cat, by = c('f_idx', 'variable', 'val'),
                         sort = FALSE, allow.cartesian = TRUE, all.y = TRUE)
        rm(preds_cat)
        psi_cat[is.na(prob), prob := 0][, lik := log(prob)]
        psi_cat[, lik := sum(lik), by = .(obs, f_idx)]
        psi_cat <- unique(psi_cat[, .(f_idx, obs, lik)])
      }
    }
    # Put it together
    psi_x <- rbind(psi_cnt, psi_cat)
    psi_x <- unique(psi_x[, lik := sum(lik), by = .(obs, f_idx)])
    return(psi_x)
  }
  if (isTRUE(parallel)) {
    out <- foreach(fold = seq_len(k), .combine = rbind) %dopar% lik_fn(fold, arf)
  } else {
    out <- foreach(fold = seq_len(k), .combine = rbind) %do% lik_fn(fold, arf)
  }
  
  # Compute per-sample likelihoods
  out <- merge(out, omega, by = 'f_idx', sort = FALSE)
  out[, ls := lik + log(wt)]
  out <- out[, matrixStats::logSumExp(ls), by = obs]
  
  # Anybody missing?
  zeros <- setdiff(seq_len(n), out[, obs])
  if (length(zeros) > 0L) {
    zero_dt <- data.table(obs = zeros, V1 = -Inf)
    out <- rbind(out, zero_dt)
  }
  
  # Export
  if (!isTRUE(log)) {
    out[, V1 := exp(V1)]
  }
  return(out[order(obs), V1])
}


