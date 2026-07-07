#' Forests for Density Estimation
#' 
#' Uses a pre-trained ARF model to estimate leaf and distribution parameters.
#' 
#' @param arf Pre-trained \code{\link{adversarial_rf}}. Alternatively, any 
#'   object of class \code{ranger}.
#' @param x Training data for estimating parameters.
#' @param oob Only use out-of-bag samples for parameter estimation? If 
#'   \code{TRUE}, \code{x} must be the same dataset used to train \code{arf}. 
#'   Set to \code{"inbag"} to only use in-bag samples. Default is \code{FALSE}, 
#'   i.e. use all observations.
#' @param family Distribution to use for density estimation of continuous 
#'   features. Current options include truncated normal (the default
#'   \code{family = "truncnorm"}) and uniform (\code{family = "unif"}). See 
#'   Details.
#' @param finite_bounds Impose finite bounds on all continuous variables? If
#'   \code{"local"}, infinite bounds are set to empirical extrema within leaves.
#'   If \code{"global"}, infinite bounds are set to global empirical extrema. 
#'   if \code{"no"} (the default), infinite bounds are left unchanged.
#' @param alpha Optional pseudocount for Laplace smoothing of categorical 
#'   features. This avoids zero-mass points when test data fall outside the 
#'   support of training data. Effectively parameterizes a flat Dirichlet prior
#'   on multinomial likelihoods.
#' @param epsilon Optional slack parameter on empirical bounds when 
#'   \code{finite_bounds != "no"}. This avoids zero-density points when test 
#'   data fall outside the support of training data. The gap between lower and 
#'   upper bounds is expanded by a factor of \code{1 + epsilon}. 
#' @param parallel Compute in parallel? Requires a registered \code{foreach}
#'   backend (\code{doParallel}, \code{doFuture}) or active \code{mirai}
#'   daemons. See \code{\link{arf-options}}.
#'   
#'   
#' @details 
#' \code{forde} extracts leaf parameters from a pretrained forest and learns
#' distribution parameters for data within each leaf. The former includes 
#' coverage (proportion of data falling into the leaf) and split criteria. The 
#' latter includes proportions for categorical features and mean/variance for
#' continuous features. The result is a probabilistic circuit, stored as a 
#' \code{data.table}, which can be used for various downstream inference tasks.
#' 
#' Currently, \code{forde} only provides support for a limited number of 
#' distributional families: truncated normal or uniform for continuous data,
#' and multinomial for discrete data. 
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
#' A \code{list} with 5 elements: (1) parameters for continuous data; (2) 
#' parameters for discrete data; (3) leaf indices and coverage; (4) metadata on
#' variables; and (5) the data input class. This list is used for estimating 
#' likelihoods with \code{\link{lik}} and generating data with \code{\link{forge}}.
#' 
#' 
#' @references 
#' Watson, D., Blesch, K., Kapar, J., & Wright, M. (2023). Adversarial random 
#' forests for density estimation and generative modeling. In \emph{Proceedings 
#' of the 26th International Conference on Artificial Intelligence and 
#' Statistics}, pp. 5357-5375.
#' 
#' Wager, S. & Athey, S. (2018). Estimation and inference of heterogeneous 
#' treatment effects using random forests. \emph{J. Am. Stat. Assoc.}, 
#' \emph{113}(523): 1228-1242.
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
#' 
#' @seealso
#' \code{\link{arf}}, \code{\link{adversarial_rf}}, \code{\link{forge}}, 
#' \code{\link{expct}}, \code{\link{lik}}
#' 
#'
#' @export
#' @import ranger 
#' @import data.table
#' @importFrom stats predict runif
#' @importFrom foreach foreach %do% %dopar%
#' 

forde <- function(
    arf, 
    x, 
    oob = FALSE,
    family = 'truncnorm', 
    finite_bounds = c('no', 'local', 'global'),
    alpha = 0,
    epsilon = 0,
    parallel = TRUE) {
  
  # To avoid data.table check issues
  tree <- n_oob <- cvg <- leaf <- variable <- count <- sd <- value <- psi_cnt <- 
    psi_cat <- f_idx <- sigma <- new_min <- new_max <- mid <- sigma0 <- prob <- 
    val <- val_count <- level <- all_na <- i <- k <- cnt <- . <- NA_share <-
    mu <- length_emp <- max_emp <- min_emp <- inbag <- n_inbag <- NULL
  
  # Prelimz
  if (isTRUE(oob) & !nrow(x) %in% c(arf$num.samples, arf$num.samples/2)) {
    stop('Forest must be trained on x when oob = TRUE.')
  }
  if (!family %in% c('truncnorm', 'unif')) {
    stop('family not recognized.')
  }
  
  finite_bounds <- match.arg(finite_bounds)
  
  # Uniform distribution requires finite bounds
  if (family == 'unif' & finite_bounds == 'no') {
    finite_bounds <- 'local'
    warning('Density estimation with uniform distribution requires finite bounds. ',
            'Resetting finite_bounds to "local".')
  }
  
  if (alpha < 0) {
    stop('alpha must be nonnegative.')
  }
  if (epsilon < 0) {
    stop('epsilon must be nonnegative.')
  }
  
  # Prep data
  input_class <- class(x)
  x <- as.data.frame(x)
  inf_flag <- sapply(seq_along(x), function(j) any(is.infinite(x[[j]])))
  if (any(inf_flag)) {
    stop('x contains infinite values.')
  }
  n <- nrow(x)
  d <- ncol(x)
  colnames_x <- colnames(x)
  classes <- sapply(x, class)
  x <- suppressWarnings(prep_x(x))
  factor_cols <- sapply(x, is.factor)
  if (any(factor_cols)) {
    # Store levels used in rf (used for internal calculations with all-NA leaves)
    lvls_rf <- arf$forest$covariate.levels[factor_cols]
    lvl_df_rf <- data.table(variable = colnames_x[factor_cols], val = lvls_rf)[
      , .(val = unlist(val), level = seq_len(length(unlist(val)))), by = variable]
    # Store levels used in data (used for forde output to post-process synthetic data)
    lvl_df_data <- data.table(x)[, .(variable = colnames_x[factor_cols], val = lapply(.SD, levels)) ,.SDcols = factor_cols][
      , .(val = unlist(val)), by = variable]
  } else {
    lvl_df_rf <- lvl_df_data <- data.table()
  }
  names(factor_cols) <- colnames_x
  deci <- rep(NA_integer_, d) 
  if (any(!factor_cols)) {
    deci[!factor_cols] <- sapply(which(!factor_cols), function(j) {
      if (any(grepl('\\.', x[[j]]))) {
        tmp <- x[grepl('\\.', x[[j]]), j]
        out <- max(nchar(sub('.*[.]', '', tmp)))
      } else {
        out <- 0L
      }
      return(out)
    })
  }
  
  # Compute leaf bounds and coverage
  num_trees <- arf$num.trees
  # Pick the parallel backend (mirai daemons vs a registered foreach backend)
  # and report the choice. Only relevant when parallel = TRUE; see
  # arf_select_backend() in mirai_helpers.R.
  backend <- arf_select_backend(parallel)
  use_mirai <- identical(backend, 'mirai')
  if (use_mirai) {
    # Load (not attach) data.table on every daemon so its S3 methods
    # (e.g. `[.data.table`) are registered for the worker bodies.
    mirai::everywhere(requireNamespace('data.table', quietly = TRUE))
  }
  # Leaf assignments first: the per-tree bounds worker computes its tree's
  # coverage from its pred column, replacing the old n x num_trees `keep`
  # table and global merge on this process (a large-grid memory hog).
  # predict() needs prep_x's column names (the forest was trained on them).
  pred <- stats::predict(arf, x, type = 'terminalNodes')$predictions + 1L
  # Restore original column names before x is shared with workers: the psi
  # workers match x's melted column names against bnds$variable, which uses
  # colnames_x (the bounds worker accesses x by index, never by name). The old
  # code renamed only after the bounds pass, so mirai's x_shared kept prep
  # names -- a latent variable-name mismatch for pathologically named columns.
  setnames(x, colnames_x)
  # Per-tree workers live once in forde_workers.R, shared by all backends; the
  # closures below adapt them to foreach's one-argument iteration.
  bnd_fn <- function(tree) {
    arf_bnd_fn(tree, arf$forest, d, finite_bounds, factor_cols, x, epsilon,
               colnames_x, pred, arf$inbag.counts, n, oob)
  }
  if (use_mirai) {
    forest_slice <- mori::share(
      arf$forest[c('split.varIDs', 'child.nodeIDs', 'split.values')])
    x_shared <- mori::share(x)
    pred_shared <- mori::share(pred)
    inbag_shared <- if (!is.null(arf$inbag.counts)) {
      mori::share(arf$inbag.counts)
    } else NULL
    bnds <- arf_mirai_tree_map(num_trees, arf_bnd_fn, list(
      forest = forest_slice, d = d, finite_bounds = finite_bounds,
      factor_cols = factor_cols, x = x_shared, epsilon = epsilon,
      colnames_x = colnames_x, pred = pred_shared,
      inbag.counts = inbag_shared, n = n, oob = oob))
  } else if (isTRUE(parallel)) {
    bnds <- foreach(tree = seq_len(num_trees), .combine = rbind) %dopar% bnd_fn(tree)
  } else {
    bnds <- foreach(tree = seq_len(num_trees), .combine = rbind) %do% bnd_fn(tree)
  }
  # Create forest index
  setkey(bnds, tree, leaf)
  bnds[, f_idx := .GRP, by = key(bnds)]
  
  # Calculate distribution parameters for each variable: one fused dispatch
  # per tree computes continuous and categorical params together (they share
  # every input), halving the round-trips of the old separate cnt/cat passes.
  psi_fn <- function(tree) {
    arf_psi_fn(tree, arf_psi_cnt_fn, arf_psi_cat_fn, x, factor_cols, pred,
               arf$inbag.counts, n, oob, bnds, finite_bounds, epsilon, family,
               lvl_df_rf, alpha)
  }
  if (use_mirai) {
    bnds_shared <- mori::share(bnds)
    psi_pair <- arf_mirai_tree_map(num_trees, arf_psi_fn, list(
      cnt_fn = arf_psi_cnt_fn, cat_fn = arf_psi_cat_fn,
      x = x_shared, factor_cols = factor_cols, pred = pred_shared,
      inbag.counts = inbag_shared, n = n, oob = oob,
      bnds = bnds_shared, finite_bounds = finite_bounds,
      epsilon = epsilon, family = family,
      lvl_df_rf = mori::share(lvl_df_rf), alpha = alpha),
      combine = arf_combine_psi)
  } else if (isTRUE(parallel)) {
    psi_pair <- arf_combine_psi(
      foreach(tree = seq_len(num_trees)) %dopar% psi_fn(tree))
  } else {
    psi_pair <- arf_combine_psi(
      foreach(tree = seq_len(num_trees)) %do% psi_fn(tree))
  }
  # Continuous case
  if (any(!factor_cols)) {
    psi_cnt <- psi_pair$cnt
    setkey(psi_cnt, f_idx, variable)
    setcolorder(psi_cnt, c('f_idx', 'variable'))
  } else {
    psi_cnt <- data.table(f_idx = integer(), variable = character(), min = numeric(), max = numeric(),
                          mu = numeric(), sigma = numeric(), NA_share = numeric())
  }

  # Categorical case
  if (any(factor_cols)) {
    psi_cat <- psi_pair$cat
    lvl_df_rf[, level := NULL]
    setkey(psi_cat, f_idx, variable)
    setcolorder(psi_cat, c('f_idx', 'variable'))
  } else {
    psi_cat <- data.table(f_idx = integer(), variable = character(), val = character(), prob = numeric(),
                          NA_share = numeric())
  }
  
  # Add metadata, export
  psi <- list(
    'cnt' = psi_cnt, 
    'cat' = psi_cat, 
    'forest' = unique(bnds[, .(f_idx, tree, leaf, cvg)]),
    'meta' = data.table('variable' = colnames_x, 'class' = classes, 
                        'family' = fifelse(factor_cols, 'multinom', family),
                        'decimals' = deci), 
    'levels' = lvl_df_data, 
    'input_class' = input_class
  )
  return(psi)
}
