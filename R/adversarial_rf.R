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
#' @param verbose Print discriminator accuracy after each round? Will also show 
#'   additional warnings.
#' @param parallel Compute in parallel? Requires a registered \code{foreach}
#'   backend (\code{doParallel}, \code{doFuture}) or active \code{mirai}
#'   daemons. See \code{\link{arf-options}}.
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
#' Integer variables are recoded with a warning (set \code{verbose = FALSE} to 
#' silence these). Default behavior is to convert integer variables with six or
#' more unique values to numeric, while those with up to five unique values are 
#' treated as ordered factors. To override this behavior, explicitly recode 
#' integer variables to the target type prior to training.
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
#'
#' # ... or with mirai (shares the learned circuit across workers via mori)
#' mirai::daemons(4)
#' }
#' 
#' @seealso
#' \code{\link{arf}}, \code{\link{forde}}, \code{\link{forge}}, 
#' \code{\link{expct}}, \code{\link{lik}}
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
      x_synth <- sample_from_leaves(rf0, x_real, factor_cols = factor_cols, lvls = lvls, prep = FALSE)
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
  
  # Prune leaves to ensure min_node_size w.r.t. real data. Per-tree work lives in
  # arf_prune_tree() (prune_workers.R). This is a meaningful share of runtime once
  # ranger training is threaded (the prune loop is pure R, unaffected by
  # ranger's num.threads), so it gets the same backend treatment as the rest.
  if (isTRUE(prune)) {
    pred <- stats::predict(rf0, x_real, type = 'terminalNodes')$predictions + 1L
    prune_one <- function(b) {
      arf_prune_tree(b, rf0$forest$child.nodeIDs, pred, min_node_size)
    }
    use_mirai <- FALSE
    if (num_trees > 1) {
      backend <- arf_select_backend(parallel)
      use_mirai <- identical(backend, 'mirai')
    }
    if (use_mirai) {
      # arf_prune_tree's body is base-R, so pass it as an object (daemons need no
      # arf loaded). Share the two big read-only objects (pred is n x num_trees;
      # child.nodeIDs is the forest) once via mori. Chunk contiguously and c() so
      # the flat result keeps tree order 1..num_trees (unname: mirai_map names
      # chunks, but child.nodeIDs must stay an unnamed list).
      pred_shared <- mori::share(pred)
      child_shared <- mori::share(rf0$forest$child.nodeIDs)
      n_chunks <- max(1L, min(as.integer(mirai::status()$connections), num_trees))
      chunks <- split(seq_len(num_trees),
                      sort(rep(seq_len(n_chunks), length.out = num_trees)))
      chunk_fn <- function(trees, worker, child_nodeIDs, pred, min_node_size) {
        lapply(trees, worker, child_nodeIDs = child_nodeIDs,
               pred = pred, min_node_size = min_node_size)
      }
      # Both functions are base-R with explicit args: strip their environments
      # before shipping. chunk_fn's would otherwise be THIS frame (rf0, dat,
      # x_real: hundreds of MB serialized into every task); arf_prune_tree's
      # namespace env would force daemons to load arf.
      # See the closure note above arf_mirai_tree_map() in mirai_helpers.R.
      environment(chunk_fn) <- globalenv()
      prune_worker <- arf_prune_tree
      environment(prune_worker) <- globalenv()
      res <- mirai::mirai_map(chunks, chunk_fn,
                              .args = list(worker = prune_worker,
                                           child_nodeIDs = child_shared,
                                           pred = pred_shared,
                                           min_node_size = min_node_size))[]
      rf0$forest$child.nodeIDs <- unname(do.call(c, res))
    } else if (isTRUE(parallel) && num_trees > 1) {
      rf0$forest$child.nodeIDs <- foreach(b = seq_len(num_trees)) %dopar% prune_one(b)
    } else {
      rf0$forest$child.nodeIDs <- foreach(b = seq_len(num_trees)) %do% prune_one(b)
    }
  }
  
  # Export
  rf0$acc <- acc
  return(rf0)
}


