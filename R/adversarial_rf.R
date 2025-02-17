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
#' @param prune Impose \code{min_node_size} by pruning? 
#' @param verbose Print discriminator accuracy after each round? Will also show additional warnings.
#' @param parallel Compute in parallel? Must register backend beforehand, e.g. 
#'   via \code{doParallel} or \code{doFuture}; see examples.
#' @param ... Extra parameters to be passed to \code{ranger}.
#' 
#' @details 
#' The adversarial random forest (ARF) algorithm partitions data into fully
#' factorized leaves where features are jointly independent. ARFs are trained
#' iteratively, with alternating rounds of generation and discrimination. In 
#' the first instance, synthetic data is generated via independent bootstraps of 
#' each feature, and a RF classifier is trained to distinguish between real and 
#' fake samples. In subsequent rounds, synthetic data is generated separately in 
#' each leaf, using splits from the previous forest. This creates increasingly 
#' realistic data that satisfies local independence by construction. The 
#' algorithm converges when a RF cannot reliably distinguish between the two 
#' classes, i.e. when OOB accuracy falls below 0.5 + \code{delta}. 
#' 
#' ARFs are useful for several unsupervised learning tasks, such as density
#' estimation (see \code{\link{forde}}) and data synthesis (see 
#' \code{\link{forge}}). For the former, we recommend increasing the number of 
#' trees for improved performance (typically on the order of 100-1000 depending 
#' on sample size).
#' 
#' Integer variables are recoded with a warning. Default behavior is to convert
#' those with six or more unique values to numeric, while those with up to five
#' unique values are treated as ordered factors. To override this behavior, 
#' explicitly recode integer variables to the target type prior to training.
#' 
#' Note: convergence is not guaranteed in finite samples. The \code{max_iters} 
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
#' Watson, D., Blesch, K., Kapar, J., & Wright, M. (2023). Adversarial random 
#' forests for density estimation and generative modeling. In \emph{Proceedings 
#' of the 26th International Conference on Artificial Intelligence and 
#' Statistics}, pp. 5357-5375.
#' 
#' 
#' @examples
#' # Train ARF and estimate leaf parameters
#' arf <- adversarial_rf(iris)
#' psi <- forde(arf, iris)
#' 
#' # Generate 100 synthetic samples from the iris dataset
#' x_synth <- forge(psi, n_synth = 100)
#'
#' # Condition on Species = "setosa" and Sepal.Length > 6
#' evi <- data.frame(Species = "setosa",
#'                   Sepal.Length = "(6, Inf)")
#' x_synth <- forge(psi, n_synth = 100, evidence = evi)
#' 
#' # Estimate average log-likelihood
#' ll <- lik(psi, iris, arf = arf, log = TRUE)
#' mean(ll)
#' 
#' # Expectation of Sepal.Length for class setosa
#' evi <- data.frame(Species = "setosa")
#' expct(psi, query = "Sepal.Length", evidence = evi)
#' 
#' \dontrun{
#' # Parallelization with doParallel
#' doParallel::registerDoParallel(cores = 4)
#'
#' # ... or with doFuture
#' doFuture::registerDoFuture()
#' future::plan("multisession", workers = 4)
#' }
#' 
#' @seealso
#' \code{\link{arf}}, \code{\link{forde}}, \code{\link{forge}}, \code{\link{expct}}, \code{\link{lik}}
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
    prune = TRUE,
    verbose = TRUE,
    parallel = TRUE,
    ...) {
  
  # To avoid data.table check issues
  i <- b <- cnt <- obs <- tree <- leaf <- N <- . <- NULL
  
  # Prep data
  x_real <- prep_x(x, verbose)
  n <- nrow(x_real)
  d <- ncol(x_real)
  factor_cols <- sapply(x_real, is.factor)
  lvls <- lapply(x_real[factor_cols], levels)
  
  # Fit initial model: sample from marginals, concatenate data, train RF
  x_synth <- setDF(lapply(x_real, sample, n, replace = TRUE))
  dat <- rbind(data.frame(y = 1L, x_real),
               data.frame(y = 0L, x_synth))
  if (isTRUE(parallel)) {
    num.threads <- NULL
  } else {
    num.threads <- 1L
  }
  if (utils::packageVersion("ranger") >= "0.16.1") {
    min.bucket <- c(min_node_size, 0)
  } else {
    min.bucket <- min_node_size
  }
  rf0 <- ranger(y ~ ., dat, keep.inbag = TRUE, classification = TRUE, 
                num.trees = num_trees, min.bucket = min.bucket, 
                respect.unordered.factors = TRUE, num.threads = num.threads, ...)
  
  # Recurse
  iters <- 0L
  acc <- acc0 <- 1 - rf0$prediction.error
  if (isTRUE(verbose)) {
    cat(paste0('Iteration: ', iters, 
               ', Accuracy: ', round(acc0 * 100, 2), '%\n'))
  }
  if (acc0 > 0.5 + delta & iters < max_iters) {
    converged <- FALSE
    while (!isTRUE(converged)) { # Adversarial loop begins...
      # Create synthetic data by sampling from intra-leaf marginals
      nodeIDs <- stats::predict(rf0, x_real, type = 'terminalNodes')$predictions
      tmp <- data.table('tree' = rep(seq_len(num_trees), each = n), 
                        'leaf' = as.integer(nodeIDs))
      tmp2 <- tmp[sample(.N, n, replace = TRUE)]
      tmp2 <- unique(tmp2[, cnt := .N, by = .(tree, leaf)])
      draw_from <- rbindlist(lapply(seq_len(num_trees), function(b) {
        x_real_b <- cbind(x_real, tmp[tree == b])
        x_real_b[, factor_cols] <- lapply(x_real_b[, factor_cols, drop = FALSE], as.numeric)
        merge(tmp2, x_real_b, by = c('tree', 'leaf'), 
              sort = FALSE)[, N := .N, by = .(tree, leaf)]
      }))
      rm(nodeIDs, tmp, tmp2)
      draw_params_within <- unique(draw_from, by = c('tree','leaf'))[, .(cnt, N)]
      adj_absolut_col <- rep(c(0, draw_params_within[-.N, cumsum(N)]), 
                             times = draw_params_within$cnt)
      adj_absolut <- rep(adj_absolut_col, d) + rep(seq(0, d - 1) * nrow(draw_from), each = n)
      idx_drawn_within <- ceiling(runif(n * d, 0, rep(draw_params_within$N, draw_params_within$cnt)))
      idx_drawn <- idx_drawn_within + adj_absolut
      draw_from_stacked <- unlist(draw_from[, -c('tree', 'leaf', 'cnt', 'N')], 
                                  use.names = FALSE)
      values_drawn_stacked <- data.table('col_id' = rep(seq_len(d), each = n), 
                                         'values' = draw_from_stacked[idx_drawn])
      x_synth <- as.data.table(split(values_drawn_stacked, by = 'col_id', keep.by = FALSE))
      setnames(x_synth, names(x_real))
      if (any(factor_cols)) {
        x_synth[, names(which(factor_cols))] <- lapply(names(which(factor_cols)), function(j) {
          lvls[[j]][x_synth[[j]]]
        })
      }
      rm(draw_from, draw_params_within, adj_absolut_col, 
         adj_absolut, idx_drawn_within, idx_drawn, draw_from_stacked)
      # Concatenate real and synthetic data
      dat <- rbind(data.frame(y = 1L, x_real),
                   data.frame(y = 0L, x_synth))
      # Train discriminator
      rf1 <- ranger(y ~ ., dat, keep.inbag = TRUE, classification = TRUE, 
                    num.trees = num_trees, min.bucket = min.bucket, 
                    respect.unordered.factors = TRUE, num.threads = num.threads, ...)
      # Evaluate
      acc0 <- 1 - rf1$prediction.error
      acc <- c(acc, acc0)
      iters <- iters + 1L
      plateau <- fifelse(isTRUE(early_stop), 
                         acc[iters] <= acc[iters + 1L], FALSE)
      if (acc0 <= 0.5 + delta | iters >= max_iters | plateau) {
        converged <- TRUE
      } else {
        # Discriminator becomes the new generator
        rf0 <- rf1
      }
      if (isTRUE(verbose)) {
        cat(paste0('Iteration: ', iters, 
                   ', Accuracy: ', round(acc0 * 100, 2), '%\n'))
      }
    }
  }
  
  # Prune leaves to ensure min_node_size w.r.t. real data
  if (isTRUE(prune)) {
    pred <- stats::predict(rf0, x_real, type = 'terminalNodes')$predictions + 1L
    prune <- function(tree) {
      # Nodes to prune are leaves which contain fewer than min_node_size real samples
      out <- rf0$forest$child.nodeIDs[[tree]]
      leaves <- which(out[[1]] == 0L)
      to_prune <- leaves[!(leaves %in% which(tabulate(pred[, tree]) >= min_node_size))]
      while(length(to_prune) > 0) {
        if (1 %in% to_prune) {
          # Never prune the root
          break
        }
        for (tp in to_prune) {
          # Find parent
          parent <- which((out[[1]] + 1L) == tp)
          if (length(parent) > 0) {
            # If node to prune (tp) is the left child of parent, replace left child with right child
            out[[1]][parent] <- out[[2]][parent]
          } else {
            # If node to prune (tp) is the right child of parent, replace right child with left child
            parent <- which((out[[2]] + 1L) == tp)
            out[[2]][parent] <- out[[1]][parent]
          }
        }
        # If both children of a parent are to be pruned, prune the parent in the next round
        # This happens if both children have been pruned
        to_prune <- which((out[[1]] + 1L) %in% to_prune)
      }
      return(out)
    }
    if (isTRUE(parallel)) {
      rf0$forest$child.nodeIDs <- foreach(b = seq_len(num_trees)) %dopar% prune(b)
    } else {
      rf0$forest$child.nodeIDs <- foreach(b = seq_len(num_trees)) %do% prune(b)
    }
  }
  
  # Export
  rf0$acc <- acc
  return(rf0)
}


