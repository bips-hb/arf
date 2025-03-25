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
#'   Defaults to \code{nrow(evidence) / num_registered_workers} for 
#'   \code{parallel == TRUE}.
#' @param parallel Compute in parallel? Must register backend beforehand, e.g. 
#'   via \code{doParallel} or \code{doFuture}; see Examples.
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
      cparams <- cforde(params, evidence_part, evidence_row_mode, nomatch, verbose, 
                        stepsize_cforde, parallel_cforde)
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
    omega[, idx := .I]
    
    synth_cnt <- synth_cat <- NULL
    # Continuous data
    if (any(!factor_cols)) {
      if (is.null(cparams) || nrow(cparams$cnt) == 0){
        psi_cond <- data.table()
      } else {
        psi_cond <- merge(omega, cparams$cnt[variable %in% query, -c("cvg_factor", "f_idx_uncond")], by = c('c_idx', 'f_idx'), 
                          sort = FALSE, allow.cartesian = TRUE)[prob > 0,]
        # calculate absolute weights for sub-leaf areas (resulting from within-row or-conditions)
        if(any(psi_cond[,prob != 1])) {
          psi_cond[, wt := wt*prob]
          psi_cond[, I := seq_len(.N), by = .(variable, idx)]
        } else {
          psi_cond[, I := 1]
        }
        psi_cond[, prob := NULL]
      } 
      psi <- unique(rbind(psi_cond,
                          merge(omega, params$cnt[variable %in% query, ], by.x = 'f_idx_uncond', by.y = 'f_idx',
                                sort = FALSE, allow.cartesian = TRUE)[,`:=` (val = NA_real_, I = 1)]), by = c("c_idx", "f_idx", "variable", "I"))[, I := NULL]
      psi[NA_share == 1, wt := 0]
      cnt <- psi[is.na(val), val := sum(wt * mu)/sum(wt), by = .(c_idx, variable)]
      cnt <- unique(cnt[, .(c_idx, variable, val)])
      synth_cnt <- dcast(cnt, c_idx ~ variable, value.var = 'val')[, c_idx := NULL]
    }
    
    
    # Categorical data
    if (any(factor_cols)) {
      if (is.null(cparams) || nrow(cparams$cat) == 0) {
        psi <- merge(omega, params$cat[variable %in% query, ], by.x = 'f_idx_uncond', by.y = 'f_idx', sort = FALSE, allow.cartesian = TRUE)
      } else {
        psi_cond <- merge(omega, cparams$cat[variable %in% query, -c("cvg_factor", "f_idx_uncond")], by = c('c_idx', 'f_idx'), 
                          sort = FALSE, allow.cartesian = TRUE)
        psi_uncond <- merge(omega, params$cat[variable %in% query, ], by.x = 'f_idx_uncond', by.y = 'f_idx',
                            sort = FALSE, allow.cartesian = TRUE)
        psi_uncond_relevant <- psi_uncond[!psi_cond[,.(idx, variable)], on = .(idx, variable), all = FALSE]
        psi <- rbind(psi_cond, psi_uncond_relevant)
      }
      psi[NA_share == 1, wt := 0]
      cat <- psi[, sum(wt * prob), by = .(c_idx, variable, val)]
      cat <- setDT(cat)[, .SD[which.max.random(V1)], by = .(c_idx, variable)]
      synth_cat <- dcast(cat, c_idx ~ variable, value.var = 'val')[, c_idx := NULL]
    }
    
    # Create dataset with expectations
    x_synth <- cbind(synth_cnt, synth_cat)
    x_synth <- post_x(x_synth, params, round)
    
    if (evidence_row_mode == "separate" & any(omega[, is.na(f_idx)])) {
      setDT(x_synth)
      indices_na <- cparams$forest[is.na(f_idx), c_idx]
      indices_sampled <- cparams$forest[!is.na(f_idx), unique(c_idx)]
      rows_na <- dcast(rbind(data.table(c_idx = 0, variable = params$meta[,variable]),
                             cparams$evidence_prepped[c_idx == indices_na,],
                             fill = T),
                       c_idx ~ variable, value.var = "val")[c_idx != 0,]
      if (nomatch == "force") {
        rows_na_sampled <- expct(params, parallel = parallel, stepsize = stepsize)
        rows_na[is.na(rows_na)] <- rows_na_sampled[is.na(rows_na[,-1])]
      }
      x_synth[, c_idx := indices_sampled]
      x_synth <- rbind(x_synth, rows_na, fill = T)
      setorder(x_synth, c_idx)[, c_idx :=  NULL]
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