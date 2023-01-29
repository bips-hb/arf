#' Adversarial Random Forests
#' 
#' Implements an adversarial random forest to learn independence-inducing splits.
#' 
#' @param x Input data. Integer variables are recoded as ordered factors with
#'   a warning. See Details.
#' @param num_trees Number of trees to grow in each forest. The default works 
#'   well for most generative modeling tasks, but should be increased for 
#'   likelihood estimation. See Details.
#' @param min_node_size Minimal number of real data samples in leaf nodes.
#' @param delta Tolerance parameter. Algorithm converges when OOB accuracy is
#'   < 0.5 + \code{delta}. 
#' @param max_iters Maximum iterations for the adversarial loop.
#' @param early_stop Terminate loop if performance fails to improve from one 
#'   round to the next? 
#' @param verbose Print discriminator accuracy after each round?
#' @param parallel Compute in parallel? Must register backend beforehand, e.g. 
#'   via \code{doParallel}.
#' @param ... Extra parameters to be passed to \code{ranger}.
#' 
#' @details 
#' The adversarial random forest (ARF) algorithm partitions data into fully
#' factorized leaves where features are jointly independent. ARFs are trained
#' iteratively, with alternating rounds of generation and discrimination. In 
#' the first instance, synthetic data is generated via independent bootstraps of 
#' each feature, and a RF classifier is trained to distinguish between real and 
#' synthetic samples. In subsequent rounds, synthetic data is generated
#' separately in each leaf, using splits from the previous forest. This creates 
#' increasingly realistic data that satisfies local independence by construction. 
#' The algorithm converges when a RF cannot reliably distinguish between the two 
#' classes, i.e. when OOB accuracy falls below 0.5 + \code{delta}. 
#' 
#' ARFs are useful for several unsupservised learning tasks, such as density
#' estimation (see \code{\link{forde}}) and data synthesis (see 
#' \code{\link{forge}}). For the former, we recommend increasing the number of 
#' trees for improved performance (typically on the order of 100-1000 depending 
#' on sample size).
#' 
#' Integer variables are treated as ordered factors by default. If the ARF is
#' passed to \code{forde}, the estimated distribution for these variables will
#' only have support on observed factor levels (i.e., the output will be a pmf,
#' not a pdf). To override this behavior and assign nonzero density to 
#' intermediate values, explicitly recode the features as numeric. 
#' 
#' Note: convergence is not guaranteed in finite samples. The \code{max_iter} 
#' argument sets an upper bound on the number of training rounds. Similar 
#' results may be attained by increasing \code{delta}. Even a single round can 
#' often give good performance, but data with strong or complex dependencies may 
#' require more iterations. With the default \code{early_stop = TRUE}, the 
#' adversarial loop terminates if performance does not improve from one round 
#' to the next, in which case further training may be pointless. 
#' 
#' 
#' @return 
#' A random forest object of class \code{ranger}.
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
#' 
#' 
#' @seealso
#' \code{\link{forde}}, \code{\link{forge}}
#' 
#' 
#' @export
#' @import ranger 
#' @import data.table
#' @importFrom stats predict
#' @importFrom foreach foreach %do% %dopar%
#'

adversarial_rf <- function(
    x, 
    num_trees = 10L, 
    min_node_size = 2L, 
    delta = 0,
    max_iters = 10L,
    early_stop = TRUE,
    verbose = TRUE,
    parallel = TRUE,
    ...) {
  
  # To avoid data.table check issues
  i <- cnt <- obs <- tree <- leaf <- . <- NULL
  
  # Prep data
  x_real <- as.data.frame(x)
  n <- nrow(x_real)
  if ('y' %in% colnames(x_real)) {
    colnames(x_real)[which(colnames(x_real) == 'y')] <- col_rename(x_real, 'y')
  }
  if ('obs' %in% colnames(x_real)) {
    colnames(x_real)[which(colnames(x_real) == 'obs')] <- col_rename(x_real, 'obs')
  }
  if ('tree' %in% colnames(x_real)) {
    colnames(x_real)[which(colnames(x_real) == 'tree')] <- col_rename(x_real, 'tree')
  } 
  if ('leaf' %in% colnames(x_real)) {
    colnames(x_real)[which(colnames(x_real) == 'leaf')] <- col_rename(x_real, 'leaf')
  } 
  idx_char <- sapply(x_real, is.character)
  if (any(idx_char)) {
    x_real[, idx_char] <- as.data.frame(
      lapply(x_real[, idx_char, drop = FALSE], as.factor)
    )
  }
  idx_logical <- sapply(x_real, is.logical)
  if (any(idx_logical)) {
    x_real[, idx_logical] <- as.data.frame(
      lapply(x_real[, idx_logical, drop = FALSE], as.factor)
    )
  }
  idx_intgr <- sapply(x_real, is.integer)
  if (any(idx_intgr)) {
    warning('Recoding integer data as ordered factors. To override this behavior, ',
            'explicitly code these variables as numeric.')
    for (j in which(idx_intgr)) {
      lvls <- sort(unique(x_real[, j]))
      x_real[, j] <- factor(x_real[, j], levels = lvls, ordered = TRUE)
    }
  }
  factor_cols <- sapply(x_real, is.factor)
  if (any(!factor_cols) & min_node_size == 1L) {
    warning('Variance is undefined when a leaf contains just a single observation. ', 
            'Consider increasing min_node_size.')
  }
  # Sample from marginals to get naive synthetic data
  x_synth <- as.data.frame(lapply(x_real, sample, n, replace = TRUE))
  # Merge real and synthetic data
  dat <- rbind(data.frame(y = 1L, x_real),
               data.frame(y = 0L, x_synth))
  # Train unsupervised random forest
  if (isTRUE(parallel)) {
    rf0 <- ranger(y ~ ., dat, keep.inbag = TRUE, classification = TRUE, 
                  num.trees = num_trees, min.node.size = 2L * min_node_size, 
                  respect.unordered.factors = TRUE, ...)
  } else {
    rf0 <- ranger(y ~ ., dat, keep.inbag = TRUE, classification = TRUE, 
                  num.trees = num_trees, min.node.size = 2L * min_node_size, 
                  respect.unordered.factors = TRUE, num.threads = 1L, ...)
  }
  
  # Recurse
  iters <- 0L
  acc <- acc0 <- 1 - rf0$prediction.error
  if (isTRUE(verbose)) {
    cat(paste0('Iteration: ', iters, 
               ', Accuracy: ', round(acc0 * 100, 2), '%\n'))
  }
  if (acc0 > 0.5 + delta & iters < max_iters) {
    sample_by_class <- function(x, n) {
      if (is.numeric(x)) {
        as.numeric(sample(x, n, replace = TRUE))
      } else {
        sample(x, n, replace = TRUE)
      }
    }
    converged <- FALSE
    while (!isTRUE(converged)) {
      # Create synthetic data
      nodeIDs <- stats::predict(rf0, x_real, type = 'terminalNodes')$predictions
      tmp <- melt(as.data.table(nodeIDs), measure.vars = 1:num_trees,
                  variable.name = 'tree', value.name = 'leaf')
      tmp[, tree := as.numeric(gsub('V', '', tree))][, obs := rep(1:n, num_trees)]
      x_real_dt <- as.data.table(x_real)[, obs := 1:n] 
      x_real_dt <- merge(x_real_dt, tmp, by = 'obs', sort = FALSE)
      tmp[, obs := NULL]
      tmp <- tmp[sample(.N, n, replace = TRUE)]
      tmp <- unique(tmp[, cnt := .N, by = .(tree, leaf)])
      draw_from <- merge(tmp, x_real_dt, by = c('tree', 'leaf'), sort = FALSE)
      x_synth <- draw_from[, lapply(.SD[, -c('cnt', 'obs')], sample_by_class, .SD[, max(cnt)]), 
                           by = .(tree, leaf)][, c('tree', 'leaf') := NULL]
      rm(nodeIDs, tmp, x_real_dt, draw_from)
      # Merge real and synthetic data
      dat <- rbind(data.frame(y = 1L, x_real),
                   data.frame(y = 0L, x_synth))
      # Train unsupervised random forest
      if (isTRUE(parallel)) {
        rf1 <- ranger(y ~ ., dat, keep.inbag = TRUE, classification = TRUE, 
                      num.trees = num_trees, min.node.size = 2 * min_node_size, 
                      respect.unordered.factors = TRUE, ...)
      } else {
        rf1 <- ranger(y ~ ., dat, keep.inbag = TRUE, classification = TRUE, 
                      num.trees = num_trees, min.node.size = 2 * min_node_size, 
                      respect.unordered.factors = TRUE, num.threads = 1, ...)
      }
      # Evaluate
      acc0 <- 1 - rf1$prediction.error
      acc <- c(acc, acc0)
      iters <- iters + 1L
      plateau <- ifelse(isTRUE(early_stop), 
                        acc[iters] <= acc[iters + 1L], FALSE)
      if (acc0 <= 0.5 + delta | iters >= max_iters | plateau) {
        converged <- TRUE
      } else {
        rf0 <- rf1
      }
      if (isTRUE(verbose)) {
        cat(paste0('Iteration: ', iters, 
                   ', Accuracy: ', round(acc0 * 100, 2), '%\n'))
      }
    }
  }
  
  # Prune leaves to ensure min_node_size w.r.t. real data
  pred <- stats::predict(rf0, x_real, type = 'terminalNodes')$predictions + 1L
  for (tree in 1:num_trees) {
    leaves <- which(rf0$forest$child.nodeIDs[[tree]][[1]] == 0L)
    to_prune <- leaves[!(leaves %in% which(tabulate(pred[, tree]) >= min_node_size))]
    while(length(to_prune) > 0) {
      for (tp in to_prune) {
        # Find parents
        parent <- which((rf0$forest$child.nodeIDs[[tree]][[1]] + 1L) == tp)
        if (length(parent) > 0) {
          # Left child
          rf0$forest$child.nodeIDs[[tree]][[1]][parent] <- rf0$forest$child.nodeIDs[[tree]][[2]][parent]
        } else {
          # Right child
          parent <- which((rf0$forest$child.nodeIDs[[tree]][[2]] + 1L) == tp)
          rf0$forest$child.nodeIDs[[tree]][[2]][parent] <- rf0$forest$child.nodeIDs[[tree]][[1]][parent]
        }
      }
      to_prune <- which((rf0$forest$child.nodeIDs[[tree]][[1]] + 1L) %in% to_prune)
    }
  }
  
  # Export
  rf0$acc <- acc
  return(rf0)
}


