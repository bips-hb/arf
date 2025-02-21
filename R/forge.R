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
#'   Options are to force sampling from a random leaf, either with a warning 
#'   (\code{"force_warning"}) or without (\code{"force"}); or to return 
#'   \code{NA}, also with or without a warning (\code{"na_warning"} and 
#'   (\code{"na"}, respectively). The default is \code{"force_warning"}.
#' @param stepsize How many rows of evidence should be handled at each step? 
#'   Defaults to \code{nrow(evidence) / num_registered_workers} for 
#'   \code{parallel == TRUE}.
#' @param parallel Compute in parallel? Must register backend beforehand, e.g. 
#'   via \code{doParallel} or \code{doFuture}; see examples.
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
    nomatch = c("force_warning", "force", "na_warning", "na"),
    stepsize = 0,
    parallel = TRUE) {
  
  evidence_row_mode <- match.arg(evidence_row_mode)
  nomatch <- match.arg(nomatch)
  
  # To avoid data.table check issues
  tree <- cvg <- leaf <- idx <- family <- mu <- sigma <- prob <- dat <- 
    variable <- relation <- wt <- j <- f_idx <- val <- . <- step_ <- c_idx <-
    f_idx_uncond <- N <- step <- V1 <- NULL
  
  factor_cols <- params$meta[, family == 'multinom']
  
  # Prepare evidence and stepsize
  if (is.null(evidence)) {
    step_no <- 1
  } else {
    evidence <- as.data.table(evidence)
    if (stepsize == 0) {
      if (parallel) {
        stepsize <- ceiling(nrow(evidence)/foreach::getDoParWorkers())
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
  
  # Run in parallel for each step
  par_fun <- function(step_) {
    
    # Prepare the event space
    if (is.null(evidence) || ( ncol(evidence) == 2 && all(colnames(evidence) == c("f_idx", "wt")))) {
      cparams <- NULL
    } else {
      # Call cforde with part of the evidence for this step
      index_start <- (step_-1)*stepsize + 1
      index_end <- min(step_*stepsize, nrow(evidence))
      evidence_part <- evidence[index_start:index_end,]
      cparams <- cforde(params, evidence_part, evidence_row_mode, nomatch, stepsize_cforde, parallel_cforde)
      if (is.null(cparams)) {
        n_synth <- n_synth * nrow(evidence_part)
      }
    } 

    # omega contains the weight (wt) for each leaf (f_idx) for each condition (c_idx)
    if (is.null(cparams)) {
      if (is.null(evidence)) {
        num_trees <- params$forest[, max(tree)]
        omega <- params$forest[, .(f_idx, f_idx_uncond = f_idx, cvg)]
        omega[, `:=` (c_idx = 1, wt = cvg / num_trees)]
        omega[, cvg := NULL]
      } else {
        omega <- copy(evidence)
        omega[, f_idx_uncond := f_idx]
        omega[, c_idx := 1]
      }
    } else {
      omega <- cparams$forest[, .(c_idx, f_idx, f_idx_uncond, wt = cvg)]
    } 
    omega <- omega[wt > 0, ]
    
    # For each synthetic sample and condition, draw a leaf according to the leaf weights
    if (nrow(omega) == 1) {
      omega <- omega[rep(1, n_synth),][, idx := .I]
    } else {
      if (evidence_row_mode == "or") {
        draws <- omega[, .(f_idx = resample(f_idx, size = n_synth, replace = TRUE, prob = wt))]
        omega <- merge(draws, omega, by = "f_idx", sort = FALSE)[, idx := .I]
      } else {
        draws <- omega[, .(f_idx = resample(f_idx, size = n_synth, replace = TRUE, prob = wt)), by = c_idx]
        omega <- merge(draws, omega, by = c("c_idx", "f_idx"), sort = FALSE)[, idx := .I]
      }
      setcolorder(omega, "idx")
    }
    
    # Simulate continuous data
    synth_cnt <- synth_cat <- NULL
    if (any(!factor_cols)) {
      fam <- params$meta[family != 'multinom', unique(family)]
      if (is.null(cparams)) {
        psi_cond <- data.table()
      } else {
        psi_cond <- merge(omega, cparams$cnt[,-c("cvg_factor", "f_idx_uncond")], by = c('c_idx', 'f_idx'), 
                          sort = FALSE, allow.cartesian = TRUE)[prob > 0,]
        # draw sub-leaf areas (resulting from within-row or-conditions)
        if(any(psi_cond[,prob != 1])) {
          psi_cond[, I := .I]
          psi_cond <- psi_cond[sort(c(psi_cond[prob == 1, I],
                          psi_cond[prob > 0 & prob < 1, fifelse(.N > 1, resample(I, 1, prob = prob), 0), by = .(variable, idx)][,V1])), -"I"]
        }
        psi_cond[, prob := NULL]
      } 
      psi <- unique(rbind(psi_cond,
                          merge(omega, params$cnt, by.x = 'f_idx_uncond', by.y = 'f_idx',
                                sort = FALSE, allow.cartesian = TRUE)[,val := NA_real_]), 
                    by = c("idx", "variable"))
      if (fam == 'truncnorm') {
        psi[is.na(val), val := truncnorm::rtruncnorm(.N, a = min, b = max, mean = mu, sd = sigma)]
        psi[is.na(val), val := mu]
      } else if (fam == 'unif') {
        psi[is.na(val), val := stats::runif(.N, min = min, max = max)]
      }
      NA_share_cnt <- psi[,.(idx, variable, NA_share)]
      synth_cnt <- dcast(psi, idx ~ variable, value.var = 'val')[, idx := NULL]
    }
    
    # Simulate categorical data
    if (any(factor_cols)) {
      if (is.null(cparams)) {
        psi <- merge(omega, params$cat, by.x = 'f_idx_uncond', by.y = 'f_idx', sort = FALSE, allow.cartesian = TRUE)
      } else {
        psi_cond <- merge(omega, cparams$cat[,-c("cvg_factor", "f_idx_uncond")], by = c('c_idx', 'f_idx'), 
                          sort = FALSE, allow.cartesian = TRUE)
        psi_uncond <- merge(omega, params$cat, by.x = 'f_idx_uncond', by.y = 'f_idx',
                            sort = FALSE, allow.cartesian = TRUE)
        psi_uncond_relevant <- psi_uncond[!psi_cond[,.(idx, variable)], on = .(idx, variable), all = FALSE]
        psi <- rbind(psi_cond, psi_uncond_relevant)
      }
      psi[prob < 1, val := sample(val, 1, prob = prob), by = .(variable, idx)]
      
      psi <- unique(psi[, .(idx, variable, val, NA_share)])
      NA_share_cat <- psi[,.(idx, variable, NA_share)]
      synth_cat <- dcast(psi, idx ~ variable, value.var = 'val')[, idx := NULL]
    }
    
    # Combine, optionally impose constraint(s)
    x_synth <- cbind(synth_cnt, synth_cat)
    if (length(x_synth) == 0) {
      x_synth <- evidence_part[FALSE,]
    }
    
    # Clean up, export
    x_synth <- post_x(x_synth, params, round)
    
    if (sample_NAs) {
      setDT(x_synth)
      NA_share <- rbind(NA_share_cnt, NA_share_cat)
      setorder(NA_share[,variable := factor(variable, levels = params$meta[,variable])], variable, idx)
      NA_share[,dat := rbinom(.N, 1, prob = NA_share)]
      x_synth[dcast(NA_share,formula =  idx ~ variable, value.var = "dat")[,-"idx"] == 1] <- NA
      x_synth <- post_x(x_synth, params, round)
    }
    
    if (evidence_row_mode == "separate" & any(omega[, is.na(f_idx)])) {
      setDT(x_synth)
      indices_na <- cparams$forest[is.na(f_idx), c_idx]
      indices_sampled <- cparams$forest[!is.na(f_idx), unique(c_idx)]
      evidence_part_long <- dcast(rbind(data.table(c_idx = 0, variable = params$meta[,variable]),
                                        cparams$evidence_prepped,
                                        fill = T),
                                  c_idx ~ variable, value.var = "val")[c_idx != 0,-"c_idx"]
      rows_na <- evidence_part_long[indices_na, ]
      rows_na[, idx := indices_na]
      rows_na <- rbindlist(replicate(n_synth, rows_na, simplify = FALSE))
      x_synth[, idx := rep(indices_sampled, each = n_synth)]
      x_synth <- rbind(x_synth, rows_na, fill = T)
      setorder(x_synth, idx)[, idx :=  NULL]
      x_synth <- post_x(x_synth, params, round)
    }
    x_synth
  }
  if (isTRUE(parallel)) {
    x_synth_ <- foreach(step = 1:step_no, .combine = "rbind") %dopar% par_fun(step)
  } else {
    x_synth_ <- foreach(step = 1:step_no, .combine = "rbind") %do% par_fun(step)
  }
  
  return(x_synth_)
}


