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
#' require more iterations.
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
#' @importFrom foreach foreach %do% %dopar%
#'

adversarial_rf <- function(
    x, 
    num_trees = 10, 
    min_node_size = 2, 
    delta = 0,
    max_iters = 10,
    verbose = TRUE,
    parallel = TRUE,
    ...) {
  
  # Prelimz
  x_real <- as.data.frame(x)
  n <- nrow(x_real)
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
  # Sample from marginals to get naive synthetic data
  x_synth <- as.data.frame(lapply(x_real, function(x) {
    sample(x, n, replace = TRUE)
  }))
  # Merge real and synthetic data
  dat <- rbind(data.frame(y = 1L, x_real),
               data.frame(y = 0L, x_synth))
  # Train unsupervised random forest
  if (isTRUE(parallel)) {
    rf0 <- ranger(y ~ ., dat, keep.inbag = TRUE, classification = TRUE, 
                  num.trees = num_trees, min.node.size = 2 * min_node_size, 
                  respect.unordered.factors = TRUE, ...)
  } else {
    rf0 <- ranger(y ~ ., dat, keep.inbag = TRUE, classification = TRUE, 
                  num.trees = num_trees, min.node.size = 2 * min_node_size, 
                  respect.unordered.factors = TRUE, num.threads = 1, ...)
  }
  
  # Recurse
  iters <- 0L
  acc <- 1 - rf0$prediction.error
  if (isTRUE(verbose)) {
    cat(paste0('Iteration: ', iters, 
                 ', Accuracy: ', round(acc * 100, 2), '%\n'))
  }
  if (acc > 0.5 + delta & iters < max_iters) {
    converged <- FALSE
    while (!isTRUE(converged)) {
      nodeIDs <- predict(rf0, x_real, type = 'terminalNodes')$predictions
      # Create synthetic data
      trees <- sample(1:num_trees, n, replace = TRUE)
      leaves <- sapply(1:n, function(i) sample(nodeIDs[, trees[i]], 1))
      synth <- function(i) {
        idx <- nodeIDs[, trees[i]] == leaves[i]
        as.data.frame(lapply(x_real[idx, ], function(x) sample(x, 1)))
      }
      if (isTRUE(parallel)) {
        x_synth <- foreach(i = 1:n, .combine = rbind) %dopar% synth(i)
      } else {
        x_synth <- foreach(i = 1:n, .combine = rbind) %do% synth(i)
      }
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
      acc <- 1 - rf1$prediction.error
      if (acc <= 0.5 + delta | iters >= max_iters) {
        converged <- TRUE
      } else {
        rf0 <- rf1
      }
      iters <- iters + 1L
      if (isTRUE(verbose)) {
        cat(paste0('Iteration: ', iters, 
                     ', Accuracy: ', round(acc * 100, 2), '%\n'))
      }
    }
  }
  
  # Prune leaves to ensure min_node_size w.r.t. real data
  pred <- predict(rf0, x_real, type = 'terminalNodes')$predictions + 1L
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
  return(rf0)
}


