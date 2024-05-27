#' Expected Value
#' 
#' Compute the expectation of some query variable(s), optionally conditioned
#' on some event(s).
#' 
#' @param params Circuit parameters learned via \code{\link{forde}}. 
#' @param query Optional character vector of variable names. Estimates will be
#'   computed for each. If \code{NULL}, all variables other than those in 
#'   \code{evidence} will be estimated. If evidence contains \code{NA}s, those
#'   variables will be estimated and a full dataset is returned.
#' @param evidence Optional set of conditioning events. This can take one of 
#'   three forms: (1) a partial sample, i.e. a single row of data with
#'   some but not all columns; (2) a data frame of conditioning events, 
#'   which allows for inequalities and intervals; or (3) a posterior distribution over leaves;
#'   see Details and Examples.
#' @param evidence_row_mode Interpretation of rows in multi-row evidence. If \code{'separate'},
#'   each row in \code{evidence} is a separate conditioning event for which \code{n_synth} synthetic samples
#'   are generated. If \code{'or'}, the rows are combined with a logical or; see Examples.
#' @param stepsize Stepsize defining number of evidence rows handled in one for each step.
#'   Defaults to nrow(evidence)/num_registered_workers for \code{parallel == TRUE}.
#' @param parallel Compute in parallel? Must register backend beforehand, e.g. 
#'   via \code{doParallel}.
#'   
#' @details 
#' This function computes expected values for any subset of features, optionally 
#' conditioned on some event(s). 
#' 
#' There are three methods for (optionally) encoding conditioning events via the 
#' \code{evidence} argument. The first is to provide a partial sample, where
#' some columns from the training data are missing or set to \code{NA}. The second is to 
#' provide a data frame with condition events. This supports inequalities and intervals. 
#' Alternatively, users may directly input a pre-calculated posterior 
#' distribution over leaves, with columns \code{f_idx} and \code{wt}. This may 
#' be preferable for complex constraints. See Examples.
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
#' # Condition on first two data rows with some missing values
#' evi <- iris[1:2,]
#' evi[1, 1] <- NA_real_
#' evi[1, 5] <- NA_character_
#' evi[2, 2] <- NA_real_
#' x_synth <- expct(psi, evidence = evi)
#' 
#' @seealso
#' \code{\link{arf}}, \code{\link{adversarial_rf}}, \code{\link{forde}}, \code{\link{forge}}, \code{\link{lik}}
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
    stepsize = 0,
    parallel = TRUE) {
  
  evidence_row_mode <- match.arg(evidence_row_mode)
  
  # To avoid data.table check issues
  variable <- tree <- f_idx <- cvg <- wt <- V1 <- value <- val <- family <-
    mu <- sigma <- obs <- prob <- f_idx_uncond <- step <- c_idx <- idx <- 
    NA_share <- . <- NULL
  
  # Prepare evidence and stepsize
  if (is.null(evidence)) {
    step_no <- 1
  } else {
    evidence <- as.data.table(evidence)
    if (ncol(evidence) == 2 && all(colnames(evidence) == c("f_idx", "wt"))) {
      stepsize <- nrow(evidence)
      step_no <- 1
    } else if (parallel & evidence_row_mode == "separate") {
      # For "separate", parallelize in forge (not in cforde)
      if (stepsize == 0) {
        stepsize <- ceiling(nrow(evidence)/foreach::getDoParWorkers())
      }
      stepsize_cforde <- 0
      parallel_cforde = FALSE
      step_no <- ceiling(nrow(evidence)/stepsize)
    } else {
      # For "or", parallelize in cforde (not in forge)
      if (stepsize == 0) {
        stepsize <- nrow(evidence)
      }
      stepsize_cforde <- stepsize
      parallel_cforde <- parallel
      stepsize <- nrow(evidence)
      step_no <- 1
    }
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
      cparams <- cforde(params, evidence_part, evidence_row_mode, stepsize_cforde, parallel_cforde)
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
        # draw sub-leaf areas (resulting from within-row or-conditions)
        if(any(psi_cond[,prob != 1])) {
          psi_cond[, I := .I]
          psi_cond <- psi_cond[sort(c(psi_cond[prob == 1, I],
                                      psi_cond[prob > 0 & prob < 1, fifelse(.N > 1, resample(I, 1, prob = prob), 0), by = .(variable, idx)][,V1])), -"I"]
        }
        psi_cond[, prob := NULL]
      } 
      psi <- unique(rbind(psi_cond,
                          merge(omega, params$cnt[variable %in% query, ], by.x = 'f_idx_uncond', by.y = 'f_idx',
                                sort = FALSE, allow.cartesian = TRUE)[,val := NA_real_]), by = c("c_idx", "f_idx", "variable"))
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
      psi[prob < 1, prob := sum(wt * prob)/sum(wt), by = .(c_idx, variable, val)]
      cat <- setDT(psi)[, .SD[which.max.random(prob)], by = .(c_idx, variable)]
      cat <- unique(cat[, .(c_idx, variable, val)]) 
      synth_cat <- dcast(cat, c_idx ~ variable, value.var = 'val')[, c_idx := NULL]
    }
    
    # Create dataset with expectations
    x_synth <- cbind(synth_cnt, synth_cat)
    x_synth <- post_x(x_synth, params)
    
    x_synth
  }
  if (isTRUE(parallel)) {
    x_synth_ <- foreach(step = 1:step_no, .combine = "rbind") %dopar% par_fun(step)
  } else {
    x_synth_ <- foreach(step = 1:step_no, .combine = "rbind") %do% par_fun(step)
  }
  
  return(x_synth_)
}


