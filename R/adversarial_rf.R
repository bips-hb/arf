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
#' @param generator Synthetic data generator. The default is to use intra-leaf 
#'   marginals (\code{generator = "bootstrap"}). Alternatively, compile the 
#'   forest into a probabilistic circuit and synthesize data from the model 
#'   (\code{generator = "forge"}). See Details.
#' @param prune Impose \code{min_node_size} by pruning? 
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
#' fake samples. In subsequent rounds, synthetic data is generated separately in 
#' each leaf, using splits from the previous forest. This can be done either by 
#' bootstrapping (\code{generator = "bootstrap"}) or fitting a circuit model 
#' (\code{generator = "forge"}). In either case, synthetic data satisfies local 
#' independence by construction. The algorithm converges when a RF cannot 
#' reliably distinguish between the two classes, i.e. when OOB accuracy falls 
#' below 0.5 + \code{delta}. 
#' 
#' ARFs are useful for several unsupervised learning tasks, such as density
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
    generator = 'bootstrap',
    prune = TRUE,
    verbose = TRUE,
    parallel = TRUE,
    ...) {
  
  # To avoid data.table check issues
  i <- b <- cnt <- obs <- tree <- leaf <- . <- NULL
  
  # Prep data
  x_real <- prep_x(x)
  n <- nrow(x_real)
  d <- ncol(x_real)
  factor_cols <- sapply(x_real, is.factor)
  lvls <- lapply(x_real[factor_cols], levels)
  if (any(!factor_cols) & min_node_size == 1L) {
    warning('Variance is undefined when a leaf contains just a single observation. ', 
            'Consider increasing min_node_size.')
  }
  if (!generator %in% c('bootstrap', 'forge')) {
    stop('generator not recognized.')
  }
  
  # Fit initial model: sample from marginals, concatenate data, train RF
  x_synth <- setDF(lapply(x_real, sample, n, replace = TRUE))
  dat <- rbind(data.frame(y = 1L, x_real),
               data.frame(y = 0L, x_synth))
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
    converged <- FALSE
    while (!isTRUE(converged)) { # Adversarial loop begins...
      if (generator == 'bootstrap') {
        # Create synthetic data by sampling from intra-leaf marginals
        nodeIDs <- stats::predict(rf0, x_real, type = 'terminalNodes')$predictions
        tmp <- data.table('tree' = rep(seq_len(num_trees), each = n), 
                          'leaf' = as.vector(nodeIDs), 
                          'obs' = rep(seq_len(n), num_trees))
        x_real_dt <- as.data.table(x_real)[, obs := seq_len(n)]
        x_real_dt <- merge(x_real_dt, tmp, by = 'obs', sort = FALSE)
        tmp[, obs := NULL]
        tmp <- tmp[sample(.N, n, replace = TRUE)]
        tmp <- unique(tmp[, cnt := .N, by = .(tree, leaf)])
        draw_from <- merge(tmp, x_real_dt, by = c('tree', 'leaf'), 
                           sort = FALSE)[, N := .N, by = .(tree, leaf)]
        rm(nodeIDs, tmp, x_real_dt)
        draw_params_within <- unique(draw_from, by = c('tree','leaf'))[, .(cnt, N)]
        adj_absolut_col <- rep(c(0, cumsum(draw_params_within$N[-nrow(draw_params_within)])), 
                               draw_params_within$cnt)
        adj_absolut <- rep(adj_absolut_col, d) + rep(seq(0, d - 1) * nrow(draw_from), each = n)
        idx_drawn_within <- ceiling(runif(n * d, 0, rep(draw_params_within$N, draw_params_within$cnt)))
        idx_drawn <- idx_drawn_within + adj_absolut
        draw_from_stacked <- unlist(draw_from[, -c('obs', 'tree', 'leaf', 'cnt', 'N')], 
                                    use.names = FALSE)
        values_drawn_stacked <- data.table('col_id' = rep(seq_len(d), each = n), 
                                           'values' = draw_from_stacked[idx_drawn])
        x_synth <- as.data.table(split(values_drawn_stacked, by = 'col_id', keep.by = FALSE))
        setnames(x_synth, names(x_real))
        if (any(factor_cols)) {
          x_synth[, (names(which(factor_cols)))] <- lapply(names(which(factor_cols)), function(x) {
            lvls[[x]][unlist(x_synth[, ..x])]
          })
        }
        rm(draw_from, draw_params_within, adj_absolut_col, 
           adj_absolut, idx_drawn_within, idx_drawn, draw_from_stacked)
      } else if (generator == 'forge') {
        # Create synthetic data using forge
        psi <- forde(rf0, x_real)
        x_synth <- forge(psi, n)
      } 
      # Concatenate real and synthetic data
      dat <- rbind(data.frame(y = 1L, x_real),
                   data.frame(y = 0L, x_synth))
      # Train discriminator
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
      out <- rf0$forest$child.nodeIDs[[tree]]
      leaves <- which(out[[1]] == 0L)
      to_prune <- leaves[!(leaves %in% which(tabulate(pred[, tree]) >= min_node_size))]
      while(length(to_prune) > 0) {
        for (tp in to_prune) {
          # Find parents
          parent <- which((out[[1]] + 1L) == tp)
          if (length(parent) > 0) {
            # Left child
            out[[1]][parent] <- out[[2]][parent]
          } else {
            # Right child
            parent <- which((out[[2]] + 1L) == tp)
            out[[2]][parent] <- out[[1]][parent]
          }
        }
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


