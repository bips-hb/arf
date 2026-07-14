#' Expected Value
#' 
#' Compute the expectation of some query variable(s), optionally conditioned
#' on some event(s).
#' 
#' @param params Circuit parameters learned via \code{\link{forde}}. 
#' @param query Optional character vector of variable names. Estimates will be
#'   computed for each. If \code{NULL}, all variables other than those in 
#'   \code{evidence} will be estimated. If \code{evidence} contains \code{NA}s, 
#'   those values will be imputed and a full dataset is returned.
#' @param evidence Optional set of conditioning events. This can take one of 
#'   three forms: (1) a partial sample, i.e. a single row of data with
#'   some but not all columns; (2) a data frame of conditioning events, 
#'   which allows for inequalities and intervals; or (3) a posterior 
#'   distribution over leaves. See Details and Examples.
#' @param evidence_row_mode Interpretation of rows in multi-row evidence. If 
#'   \code{"separate"}, each row in \code{evidence} is a unique conditioning 
#'   event for which \code{n_synth} synthetic samples are generated. If 
#'   \code{"or"}, the rows are combined with a logical OR. See Examples.
#' @param round Round continuous variables to their respective maximum precision 
#'   in the real data set?
#' @param nomatch What to do if no leaf matches a condition in \code{evidence}?
#'   Options are to force sampling from a random leaf (\code{"force"}) or return 
#'   \code{NA} (\code{"na"}). The default is \code{"force"}.
#' @param verbose Show warnings, e.g. when no leaf matches a condition?   
#' @param stepsize How many rows of evidence should be handled at each step? 
#'   Defaults to \code{nrow(evidence)} divided by the number of registered
#'   workers or daemons for 
#'   \code{parallel == TRUE}.
#' @param parallel Compute in parallel? Requires a registered \code{foreach}
#'   backend (\code{doParallel}, \code{doFuture}) or active \code{mirai}
#'   daemons. See \code{\link{arf-options}}. With
#'   \code{evidence_row_mode = "or"}, parallelization happens inside the
#'   conditional circuit computation; in benchmarks this gave little speedup
#'   while raising peak memory, so consider \code{parallel = FALSE} for large
#'   \code{"or"} queries.
#'
#' @details 
#' This function computes expected values for any subset of features, optionally 
#' conditioned on some event(s). 
#' 
#' There are three methods for (optionally) encoding conditioning events via the 
#' \code{evidence} argument. The first is to provide a partial sample, where
#' some columns from the training data are missing or set to \code{NA}. The 
#' second is to provide a data frame with condition events. This supports 
#' inequalities and intervals. Alternatively, users may directly input a 
#' pre-calculated posterior distribution over leaves, with columns \code{f_idx} 
#' and \code{wt}. This may be preferable for complex constraints. See Examples.
#' 
#' Please note that results for continuous features which are both included in 
#' \code{query} and in \code{evidence} with an interval condition are currently 
#' inconsistent.
#' 
#' @return 
#' A one row data frame with values for all query variables.
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
#' # What is the expected value of Sepal.Length?
#' expct(psi, query = "Sepal.Length")
#' 
#' # What if we condition on Species = "setosa"?
#' evi <- data.frame(Species = "setosa")
#' expct(psi, query = "Sepal.Length", evidence = evi)
#' 
#' # Compute expectations for all features other than Species
#' expct(psi, evidence = evi)
#' 
#' # Condition on Species = "setosa" and Petal.Width > 0.3
#' evi <- data.frame(Species = "setosa", 
#'                   Petal.Width = ">0.3")
#' expct(psi, evidence = evi)
#' 
#' # Condition on first two rows with some missing values
#' evi <- iris[1:2,]
#' evi[1, 1] <- NA_real_
#' evi[1, 5] <- NA_character_
#' evi[2, 2] <- NA_real_
#' x_synth <- expct(psi, evidence = evi)
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
#' \code{\link{arf}}, \code{\link{adversarial_rf}}, \code{\link{forde}}, 
#' \code{\link{forge}}, \code{\link{lik}}
#' 
#'
#' @export
#' @import data.table
#' @importFrom truncnorm etruncnorm
#' 

expct <- function(
    params, 
    query = NULL, 
    evidence = NULL,
    evidence_row_mode = c("separate", "or"),
    round = FALSE,
    nomatch = c("force", "na"),
    verbose = TRUE,
    stepsize = 0,
    parallel = TRUE) {
  
  evidence_row_mode <- match.arg(evidence_row_mode)
  nomatch <- match.arg(nomatch)
  
  # To avoid data.table check issues
  variable <- tree <- f_idx <- cvg <- wt <- V1 <- value <- val <- family <-
    mu <- sigma <- obs <- prob <- f_idx_uncond <- step <- c_idx <- idx <- 
    NA_share <- . <- NULL

  # Defaults so the extracted per-step worker always receives these (set below
  # for the conditional cases).
  stepsize_cforde <- 0L
  parallel_cforde <- FALSE
  
  # Prepare evidence and stepsize
  if (is.null(evidence)) {
    step_no <- 1
  } else {
    evidence <- as.data.table(evidence)
    if (stepsize == 0) {
      if (parallel) {
        stepsize <- ceiling(nrow(evidence)/arf_n_workers())
      } else {
        stepsize <- nrow(evidence)
      }
    } else if (stepsize > nrow(evidence)) {
      stepsize <- nrow(evidence)
    }
    if (ncol(evidence) == 2 && all(colnames(evidence) == c("f_idx", "wt"))) {
      stepsize <- nrow(evidence)
    } else if (evidence_row_mode == "separate") {
      # For "separate", parallelize in expct (not in cforde)
      stepsize_cforde <- 0
      parallel_cforde = FALSE
    } else {
      # For "or", parallelize in cforde (not in expct)
      parallel_cforde <- parallel
      stepsize_cforde <- stepsize
      parallel <- FALSE
      stepsize <- nrow(evidence)
    }
    step_no <- ceiling(nrow(evidence)/stepsize)
  } 
  
  # Check query
  if (is.null(query)) {
    if (any(is.na(evidence))) {
      query <- params$meta$variable
    } else {
      query <- setdiff(params$meta$variable, colnames(evidence))
    }
  } else if (any(!query %in% params$meta$variable)) {
    err <- setdiff(query, params$meta$variable)
    stop('Unrecognized feature(s) in query: ', err)
  }
  factor_cols <- params$meta[variable %in% query, family == 'multinom']
  
  # Per-step work lives once in arf_expct_step() (expct_workers.R). This closure
  # adapts it to foreach's one-argument iteration.
  par_fun <- function(step_) {
    arf_expct_step(step_, params, evidence, query, factor_cols,
                   evidence_row_mode, nomatch, verbose, round, stepsize,
                   stepsize_cforde, parallel_cforde)
  } 
  # Parallelism is across steps: 1 step is inherently serial for any backend.
  # mirai only for step_no > 1; "or" mode already set parallel <- FALSE.
  use_mirai <- FALSE
  if (step_no > 1) {
    backend <- arf_select_backend(parallel)
    use_mirai <- identical(backend, "mirai")
  } 
  if (use_mirai) {
    arf_load_on_daemons()  # daemons need arf: worker calls cforde/post_x/which.max.random
    params_shared <- mori::share(params)
    evidence_shared <- if (!is.null(evidence)) mori::share(evidence) else NULL
    x_synth_ <- arf_mirai_tree_map(step_no, arf_expct_step, list(
      params = params_shared, evidence = evidence_shared, query = query,
      factor_cols = factor_cols, evidence_row_mode = evidence_row_mode,
      nomatch = nomatch, verbose = verbose, round = round, stepsize = stepsize,
      stepsize_cforde = stepsize_cforde, parallel_cforde = parallel_cforde),
      # package-level combine: an inline closure here would serialize expct's
      # whole frame (params included) into every task
      # See the closure note in mirai_helpers.R.
      combine = arf_rbind_steps)
  } else if (isTRUE(parallel) && step_no > 1) {
    x_synth_ <- foreach(step = 1:step_no, .combine = "rbind") %dopar% par_fun(step)
  } else {
    x_synth_ <- foreach(step = 1:step_no, .combine = "rbind") %do% par_fun(step)
  }
  
  return(x_synth_)
}
