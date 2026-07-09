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
#' @param parallel Compute in parallel? Requires a registered \code{foreach}
#'   backend (\code{doParallel}, \code{doFuture}) or active \code{mirai}
#'   daemons. See \code{\link{arf-options}}.
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
#' # Train ARF and estimate leaf parameters
#' arf <- adversarial_rf(iris)
#' psi <- forde(arf, iris)
#' 
#' # Estimate average log-likelihood
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
#' evi <- data.frame(Species = "setosa", 
#'                   Petal.Width = ">0.3")
#' lik(psi, query = iris[1, 1:3], evidence = evi)
#' 
#' \dontrun{
#' # Parallelization with doParallel
#' doParallel::registerDoParallel(cores = 4)
#'
#' # ... or with doFuture
#' doFuture::registerDoFuture()
#' future::plan("multisession", workers = 4)
#'
#' # ... or with mirai (shares large read-only inputs across workers via mori)
#' mirai::daemons(4)
#' }
#' 
#' @seealso
#' \code{\link{arf}}, \code{\link{adversarial_rf}}, \code{\link{forde}}, \code{\link{forge}}, \code{\link{expct}}
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
    V1 <- relation <- f_idx <- wt <- val <- family <- fold <- f_idx_uncond <- 
    . <- NULL
  
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
  pure <- all(factor_cols) | all(!factor_cols)  # worker arg (was computed inside lik_fn)
  
  # Prep evidence
  conj <- !is.null(evidence) && !(ncol(evidence) == 2 && all(c("f_idx", "wt") %in% colnames(evidence)))
  
  # Check ARF
  preds <- NULL  # set below iff arf-based leaf assignment is used
  if (d == params$meta[, .N] & !is.null(arf)) {
    num_trees <- arf$num.trees
    preds <- stats::predict(arf, x, type = 'terminalNodes')$predictions + 1L
    preds <- data.table('tree' = rep(seq_len(num_trees), each = n), 
                        'leaf' = as.vector(preds),
                        'obs' = rep(seq_len(n), times = num_trees))
    if (isTRUE(oob)) {
      preds <- stats::na.omit(preds)
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
    omega <- cforde(params, evidence, "or")$forest[, .(f_idx = f_idx_uncond, wt = cvg)]
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
  
  # Per-fold work lives once in arf_lik_fold() (lik_workers.R). This closure
  # adapts it to foreach; `arf` is used only via is.null() (preds carries the
  # arf-derived leaf assignments), passed to the worker as a boolean.
  lik_fn <- function(fold, arf) {
    arf_lik_fold(fold, params, x, factor_cols, leaves, omega, preds,
                 batch_idx, pure, !is.null(arf))
  }
  # Parallelism is across folds: a single fold is inherently serial. mirai only
  # for k > 1.
  use_mirai <- FALSE
  if (k > 1) {
    backend <- arf_select_backend(parallel)
    use_mirai <- identical(backend, "mirai")
  }
  if (use_mirai) {
    arf_load_on_daemons()  # daemons need arf (worker uses bare data.table verbs)
    out <- arf_mirai_tree_map(k, arf_lik_fold, list(
      params = mori::share(params), x = mori::share(x),
      factor_cols = factor_cols, leaves = leaves, omega = mori::share(omega),
      preds = if (!is.null(preds)) mori::share(preds) else NULL,
      batch_idx = batch_idx, pure = pure, has_arf = !is.null(arf)))
  } else if (isTRUE(parallel) && k > 1) {
    out <- foreach(fold = seq_len(k), .combine = rbind) %dopar% lik_fn(fold, arf)
  } else {
    out <- foreach(fold = seq_len(k), .combine = rbind) %do% lik_fn(fold, arf)
  }
  
  # Folds return per-obs log-likelihoods, already reduced against omega inside
  # the worker (folds cover disjoint obs); nothing left to aggregate here.

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


