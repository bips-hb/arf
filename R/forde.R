#' Forests for Density Estimation
#' 
#' Uses a pre-trained ARF model to estimate leaf and distribution parameters.
#' 
#' @param arf Pre-trained adversarial random forest. Alternatively, any object
#'   of class \code{ranger}.
#' @param x Training data for estimating parameters.
#' @param oob Only use out-of-bag samples for parameter estimation? If 
#'   \code{TRUE}, \code{x} must be the same dataset used to train \code{arf}.
#' @param family Distribution to use for density estimation of continuous 
#'   features. Current options include truncated normal (the default
#'   \code{family = "truncnorm"}) and uniform (\code{family = "unif"}). See 
#'   Details.
#' @param epsilon Slack parameter on empirical bounds when \code{family = "unif"}.
#'   This avoids zero-density points when test data fall outside the support
#'   of training data. The gap between lower and upper bounds is expanded by 
#'   a factor of \code{epsilon}. Only used when a variable is never selected for
#'   splitting.
#' @param parallel Compute in parallel? Must register backend beforehand, e.g. 
#'   via \code{doParallel}.
#'   
#'   
#' @details 
#' \code{forde} extracts leaf parameters from a pretrained forest and learns
#' distribution parameters for data within each leaf. The former includes 
#' coverage (proportion of data falling into the leaf) and split criteria. The 
#' latter includes proportions for categorical features and mean/variance for
#' continuous features. These values are stored in a \code{data.table}, which 
#' can be used as input to various other functions.
#' 
#' Currently, \code{forde} only provides support for a limited number of 
#' distributional families: truncated normal or uniform for continuous data,
#' and multinomial for discrete data. Future releases will accommodate a larger 
#' class of options.
#' 
#' Though \code{forde} was designed to take an adversarial random forest as 
#' input, the function's first argument can in principle be any object of class 
#' \code{ranger}. This allows users to test performance with alternative 
#' pipelines (e.g., with supervised forest input). There is also no requirement 
#' that \code{x} be the data used to fit \code{arf}, unless \code{oob = TRUE}. 
#' In fact, using another dataset here may protect against overfitting. This 
#' connects with Wager & Athey's (2018) notion of "honest trees".
#' 
#' 
#' @return 
#' A \code{data.table} with leaf and distribution parameters for each feature. 
#' These are used for estimating likelihoods with \code{\link{lik}} and 
#' generating data with \code{\link{forge}}.
#' 
#' 
#' @references 
#' Watson, D., Blesch, K., Kapar, J., & Wright, M. (2022). Adversarial random 
#' forests for density estimation and generative modeling. \emph{arXiv} preprint,
#' 2205.09435.
#' 
#' Wager, S. & Athey, S. (2018). Estimation and inference of heterogeneous 
#' treatment effects using random forests. \emph{J. Am. Stat. Assoc.}, 
#' \emph{113}(523): 1228-1242.
#' 
#' 
#' @examples
#' arf <- adversarial_rf(iris)
#' psi <- forde(arf, iris)
#' head(psi)
#' 
#' 
#' @seealso
#' \code{\link{adversarial_rf}}, \code{\link{forge}}, \code{\link{lik}}
#' 
#'
#' @export
#' @import ranger 
#' @import data.table
#' @importFrom foreach foreach %do% %dopar%
#' 


forde <- function(
    arf, 
    x, 
    oob = FALSE,
    family = 'truncnorm', 
    epsilon = 0.1, 
    parallel = TRUE) {
  
  # To avoid data.table check issues
  tree <- n_oob <- cvg <- leaf <- variable <- count <- sd <- value <- . <- NULL
  
  # Prelimz
  if (isTRUE(oob) & !nrow(x) %in% c(arf$num.samples, arf$num.samples/2)) {
    stop('Forest must be trained on x when oob = TRUE.')
  }
  if (!family %in% c('truncnorm', 'unif')) {
    stop('family not recognized.')
  }
  
  # Prep data
  x <- as.data.frame(x)
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
  idx_intgr <- sapply(x, is.integer)
  if (any(idx_intgr)) {
    warning('Recoding integer data as ordered factors. To override this behavior, ',
            'explicitly code these variables as numeric.')
    for (j in which(idx_intgr)) {
      lvls <- sort(unique(x[, j]))
      x[, j] <- factor(x[, j], levels = lvls, ordered = TRUE)
    }
  }
  factor_cols <- sapply(x, is.factor)
  
  # Compute leaf bounds and coverage
  num_trees <- arf$num.trees
  pred <- predict(arf, x, type = 'terminalNodes')$predictions + 1L
  bnd_fn <- function(tree) {
    num_nodes <- length(arf$forest$split.varIDs[[tree]])
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
      left_child <- arf$forest$child.nodeIDs[[tree]][[1]][i] + 1L
      right_child <- arf$forest$child.nodeIDs[[tree]][[2]][i] + 1L
      splitvarID <- arf$forest$split.varIDs[[tree]][i] + 1L
      splitval <- arf$forest$split.value[[tree]][i]
      if (left_child > 1 & left_child != right_child) {
        ub[left_child, ] <- ub[right_child, ] <- ub[i, ]
        lb[left_child, ] <- lb[right_child, ] <- lb[i, ]
        ub[left_child, splitvarID] <- lb[right_child, splitvarID] <- splitval
      }
    }
    leaves <- which(arf$forest$child.nodeIDs[[tree]][[1]] == 0L) 
    colnames(lb) <- arf$forest$independent.variable.names
    colnames(ub) <- arf$forest$independent.variable.names
    merge(melt(data.table(tree = tree, leaf = leaves, lb[leaves, ]), 
               id.vars = c('tree', 'leaf'), value.name = 'min'), 
          melt(data.table(tree = tree, leaf = leaves, ub[leaves, ]), 
               id.vars = c('tree', 'leaf'), value.name = 'max'), 
          by = c('tree', 'leaf', 'variable'), sort = FALSE)
  }
  if (isTRUE(parallel)) {
    bnds <- foreach(tree = 1:num_trees, .combine = rbind) %dopar% bnd_fn(tree)
  } else {
    bnds <- foreach(tree = 1:num_trees, .combine = rbind) %do% bnd_fn(tree)
  }
  # Use only OOB data?
  if (isTRUE(oob)) {
    inbag <- (do.call(cbind, arf$inbag.counts) > 0L)[1:(arf$num.samples/2), ]
    pred[inbag] <- NA_integer_
    bnds[, n_oob := sum(!is.na(pred[, tree])), by = tree]
    bnds[, cvg := sum(pred[, tree] == leaf, na.rm = TRUE) / n_oob, by = .(tree, leaf)]
  } else {
    bnds[, cvg := sum(pred[, tree] == leaf) / n, by = .(tree, leaf)]
  }
  # Can't do anything with coverage 1/n
  if (any(!factor_cols)) {
    bnds[cvg == 1/n, cvg := 0]
  }
  
  # Calculate distribution parameters
  psi_cnt <- psi_cat <- NULL
  # Continuous case
  if (any(!factor_cols)) {
    psi_cnt_fn <- function(tree) {
      dt <- data.table(x[, !factor_cols, drop = FALSE], leaf = pred[, tree])
      if (isTRUE(oob)) {
        dt <- dt[!is.na(leaf)]
      }
      long <- melt(dt, id.vars = 'leaf', variable.factor = FALSE)
      if (family == 'truncnorm') {
        long[, list(cat = NA_character_, prob = NA_real_, 
                    mu = mean(value), sigma = sd(value), 
                    family = 'truncnorm', tree = tree), 
             by = .(leaf, variable)]
      } else if (family == 'unif') {
        long[, list(tree = tree, cat = NA_character_, prob = NA_real_, 
                    mu = NA_real_, sigma = NA_real_, family = 'unif'), 
             by = .(leaf, variable)]
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
      dt <- data.table(x[, factor_cols, drop = FALSE], leaf = pred[, tree])
      if (isTRUE(oob)) {
        dt <- dt[!is.na(leaf)]
      }
      long <- melt(dt, id.vars = 'leaf', variable.factor = FALSE,
                   value.factor = FALSE, value.name = 'cat')
      long[, count := .N, by = .(leaf, variable)]
      unique(long[, list(prob = .N/count, mu = NA_real_, sigma = NA_real_, 
                         family = 'multinom', tree = tree), 
                  by = .(leaf, variable, cat)])
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
  
  # Export
  return(psi)
}


