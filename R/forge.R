#' Forests for Generative Modeling
#' 
#' Uses pre-trained FORDE model to simulate synthetic data.
#' 
#' @param params Circuit parameters learned via \code{\link{forde}}. 
#' @param evidence Optional set of conditioning events. This can take one of 
#'   three forms: (1) a partial sample, i.e. a single row of data with some but
#'   not all columns; (2) a data frame of conditioning events, which allows for 
#'   inequalities; or (3) a posterior distribution over leaves. See Details.
#' @param evidence_row_mode Interpretation of rows in multi-row evidence. If 
#'   \code{"separate"}, each row in \code{evidence} is a unique conditioning 
#'   event for which \code{n_synth} synthetic samples are generated. If 
#'   \code{"or"}, the rows are combined with a logical OR. See Examples.
#' @param round Round continuous variables to their respective maximum precision 
#'   in the real data set?
#' @param sample_NAs Sample \code{NA}s respecting the probability for missing 
#'   values in the original data?
#' @param nomatch What to do if no leaf matches a condition in \code{evidence}?
#'   Options are to force sampling from a random leaf (\code{"force"}) or return 
#'   \code{NA} (\code{"na"}). The default is \code{"force"}.
#' @param verbose Show warnings, e.g. when no leaf matches a condition?   
#' @param stepsize How many rows of evidence should be handled at each step? 
#'   Defaults to \code{nrow(evidence) / num_registered_workers} for 
#'   \code{parallel == TRUE}.
#' @param parallel Compute in parallel? Requires a registered \code{foreach}
#'   backend (\code{doParallel}, \code{doFuture}) or active \code{mirai}
#'   daemons. See \code{\link{arf-options}}.
#' @param n_synth Number of synthetic samples to generate.
#'
#' @details  
#' \code{forge} simulates a synthetic dataset of \code{n_synth} samples. First,
#' leaves are sampled in proportion to either their coverage (if 
#' \code{evidence = NULL}) or their posterior probability. Then, each feature is 
#' sampled independently within each leaf according to the probability mass or 
#' density function learned by \code{\link{forde}}. This will create realistic 
#' data so long as the adversarial RF used in the previous step satisfies the 
#' local independence criterion. See Watson et al. (2023).
#' 
#' There are three methods for (optionally) encoding conditioning events via the 
#' \code{evidence} argument. The first is to provide a partial sample, where
#' some columns from the training data are missing or set to \code{NA}. The 
#' second is to provide a data frame with condition events. This supports 
#' inequalities and intervals. Alternatively, users may directly input a 
#' pre-calculated posterior distribution over leaves, with columns \code{f_idx} 
#' and \code{wt}. This may be preferable for complex constraints. See Examples.
#' 
#' @return  
#' A dataset of \code{n_synth} synthetic samples. 
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
#' # Condition on Species = "setosa"
#' evi <- data.frame(Species = "setosa")
#' x_synth <- forge(psi, n_synth = 100, evidence = evi)
#' 
#' # Condition on Species = "setosa" and Sepal.Length > 6
#' evi <- data.frame(Species = "setosa",
#'                   Sepal.Length = "(6, Inf)")
#' x_synth <- forge(psi, n_synth = 100, evidence = evi)
#' 
#' # Alternative syntax for </> conditions
#' evi <- data.frame(Sepal.Length = ">6")
#' x_synth <- forge(psi, n_synth = 100, evidence = evi)
#' 
#' # Negation condition, i.e. all classes except "setosa"
#' evi <- data.frame(Species = "!setosa")
#' x_synth <- forge(psi, n_synth = 100, evidence = evi)
#' 
#' # Condition on first two data rows with some missing values
#' evi <- iris[1:2,]
#' evi[1, 1] <- NA_real_
#' evi[1, 5] <- NA_character_
#' evi[2, 2] <- NA_real_
#' x_synth <- forge(psi, n_synth = 1, evidence = evi)
#' 
#' # Or just input some distribution on leaves
#' # (Weights that do not sum to unity are automatically scaled)
#' n_leaves <- nrow(psi$forest)
#' evi <- data.frame(f_idx = psi$forest$f_idx, wt = rexp(n_leaves))
#' x_synth <- forge(psi, n_synth = 100, evidence = evi)
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
#' \code{\link{arf}}, \code{\link{adversarial_rf}}, \code{\link{forde}}, 
#' \code{\link{expct}}, \code{\link{lik}}
#' 
#' @export
#' @import data.table
#' @importFrom foreach foreach %dopar% getDoParWorkers
#' @importFrom truncnorm rtruncnorm 
#' @importFrom stats rbinom
#'

forge <- function(
    params, 
    n_synth,
    evidence = NULL,
    evidence_row_mode = c("separate", "or"),
    round = TRUE,
    sample_NAs = FALSE,
    nomatch = c("force", "na"),
    verbose = TRUE,
    stepsize = 0,
    parallel = TRUE) {
  
  evidence_row_mode <- match.arg(evidence_row_mode)
  nomatch <- match.arg(nomatch)
  
  # To avoid data.table check issues
  tree <- cvg <- leaf <- idx <- family <- mu <- sigma <- prob <- dat <- 
    variable <- relation <- wt <- j <- f_idx <- val <- . <- step_ <- c_idx <-
    f_idx_uncond <- N <- step <- V1 <- NULL

  # Defaults so the extracted per-step worker always receives these (the
  # evidence branches below set them for the conditional cases).
  stepsize_cforde <- 0L
  parallel_cforde <- FALSE

  factor_cols <- params$meta[, family == 'multinom']
  
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
    }
    if (ncol(evidence) == 2 && all(colnames(evidence) == c("f_idx", "wt"))) {
      stepsize <- nrow(evidence)
    } else if (evidence_row_mode == "separate") {
      # For "separate", parallelize in forge (not in cforde)
      stepsize_cforde <- 0
      parallel_cforde = FALSE
    } else {
      # For "or", parallelize in cforde (not in forge)
      parallel_cforde <- parallel
      stepsize_cforde <- stepsize
      parallel <- FALSE
      stepsize <- nrow(evidence)
    }
    step_no <- ceiling(nrow(evidence)/stepsize)
  } 
  
  # Per-step work lives once in arf_forge_step() (forge_workers.R), shared by all
  # backends; this closure adapts it to foreach's one-argument iteration.
  par_fun <- function(step_) {
    arf_forge_step(step_, params, evidence, n_synth, factor_cols,
                   evidence_row_mode, nomatch, verbose, round, sample_NAs,
                   stepsize, stepsize_cforde, parallel_cforde)
  }
  # Parallelism is across steps, so a single step is inherently serial regardless
  # of backend (1-task %dopar% is pure overhead). Only pick a backend when
  # step_no > 1; the "or" branch already set parallel <- FALSE (cforde
  # parallelizes there instead).
  use_mirai <- FALSE
  if (step_no > 1) {
    backend <- arf_select_backend(parallel)
    use_mirai <- identical(backend, 'mirai')
  }
  if (use_mirai) {
    arf_load_on_daemons()  # daemons need arf: worker calls cforde/post_x/resample
    params_shared <- mori::share(params)
    evidence_shared <- if (!is.null(evidence)) mori::share(evidence) else NULL
    x_synth_ <- arf_mirai_tree_map(step_no, arf_forge_step, list(
      params = params_shared, evidence = evidence_shared, n_synth = n_synth,
      factor_cols = factor_cols, evidence_row_mode = evidence_row_mode,
      nomatch = nomatch, verbose = verbose, round = round,
      sample_NAs = sample_NAs, stepsize = stepsize,
      stepsize_cforde = stepsize_cforde, parallel_cforde = parallel_cforde),
      # package-level combine: an inline closure here would serialize forge's
      # whole frame (params included) into every task
      combine = arf_rbind_steps)
  } else if (isTRUE(parallel) && step_no > 1) {
    x_synth_ <- foreach(step = 1:step_no, .combine = "rbind") %dopar% par_fun(step)
  } else {
    x_synth_ <- foreach(step = 1:step_no, .combine = "rbind") %do% par_fun(step)
  }

  return(x_synth_)
}


