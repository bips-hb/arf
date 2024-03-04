#' Forests for Generative Modeling
#' 
#' Uses pre-trained FORDE model to simulate synthetic data.
#' 
#' @param params Circuit parameters learned via \code{\link{forde}}. 
#' @param condition Optional set of conditioning events.
#' @param condition_row_mode Interpretation of rows in multi-row conditions.
#' @param sample_NAs Sample NAs respecting the probability for missing values in the original data.
#' @param stepsize Stepsize defining number of condition rows handled in one for each step.
#'   Defaults to nrow(condition)/num_registered_workers for parallel == TRUE.
#' @param parallel Compute in parallel? Must register backend beforehand, e.g. 
#'   via \code{doParallel}.
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
#' some but not all columns from the training data are present. The second is to 
#' provide a data frame with three columns: \code{variable}, \code{relation}, 
#' and \code{value}. This supports inequalities via \code{relation}. 
#' Alternatively, users may directly input a pre-calculated posterior 
#' distribution over leaves, with columns \code{f_idx} and \code{wt}. This may 
#' be preferable for complex constraints. See Examples.
#' 
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
#' arf <- adversarial_rf(iris)
#' psi <- forde(arf, iris)
#' x_synth <- forge(psi, n_synth = 100)
#'
#' # Condition on Species = "setosa"
#' evi <- data.frame(Species = "setosa")
#' x_synth <- forge(psi, n_synth = 100, evidence = evi)
#' 
#' # Condition on Species = "setosa" and Sepal.Length > 6
#' evi <- data.frame(variable = c("Species", "Sepal.Length"),
#'                   relation = c("==", ">"), 
#'                   value = c("setosa", 6))
#' x_synth <- forge(psi, n_synth = 100, evidence = evi)
#' 
#' # Condition on first two data rows with some missing values
#' condition <- iris[1:2,]
#' condition[1,1] <- NA_real_
#' condition[1,5] <- NA_character_
#' condition[2,2] <- NA_real_
#' x_synth <- forge(psi, n_synth = 1, condition = condition)
#' 
#' # Or just input some distribution on leaves
#' # (Weights that do not sum to unity are automatically scaled)
#' n_leaves <- nrow(psi$forest)
#' evi <- data.frame(f_idx = psi$forest$f_idx, wt = rexp(n_leaves))
#' x_synth <- forge(psi, n_synth = 100, evidence = evi)
#' 
#'
#' @seealso
#' \code{\link{adversarial_rf}}, \code{\link{forde}}
#' 
#' 
#' @export
#' @import data.table
#' @importFrom foreach foreach %dopar% getDoPar getDoParWorkers
#' @importFrom truncnorm rtruncnorm 
#' @importFrom stats rbinom
#'

forge <- function(
    params, 
    n_synth,
    condition = NULL,
    condition_row_mode = c("separate", "or"),
    sample_NAs = F,
    stepsize = 0,
    parallel = TRUE) {
  
  condition_row_mode <- match.arg(condition_row_mode)
  
  doParRegistered <- getDoParRegistered()
  num_workers <- getDoParWorkers()
  if(!parallel & doParRegistered & (num_workers > 1)) {
    registerDoSEQ()
  }
  
  # To avoid data.table check issues
  tree <- cvg <- leaf <- idx <- family <- mu <- sigma <- prob <- dat <- 
    variable <- relation <- wt <- j <- f_idx <- val <- . <- step_ <- c_idx <-
    f_idx_uncond <- N <- NULL
  
  factor_cols <- params$meta[, family == 'multinom']
  
  if(!is.null(condition)) {
    condition <- as.data.table(condition)
    if(parallel & condition_row_mode == "separate") {
      if(stepsize == 0) {
        stepsize <- ceiling(nrow(condition)/getDoParWorkers())
      }
      stepsize_cforde <- 0
      parallel_cforde = F
      step_no <- ceiling(nrow(condition)/stepsize)
    } else {
      if(stepsize == 0) {
        stepsize = nrow(condition)
      }
      stepsize_cforde <- stepsize
      parallel_cforde <- parallel
      stepsize <- nrow(condition)
      step_no <- 1
    }
  } else {
    step_no <- 1
  }
  
  x_synth_ <- foreach(step_ = 1:step_no, .combine = "rbind") %dopar% {
    
    # Prepare the event space
    if (!is.null(condition)) {
      index_start <- (step_-1)*stepsize + 1
      index_end <- min(step_*stepsize, nrow(condition))
      condition_part <- condition[index_start:index_end,]
      cparams <- cforde(params, condition_part, condition_row_mode, stepsize_cforde, parallel_cforde)
      if(is.null(cparams)) {
        n_synth <- n_synth * nrow(condition_part)
      }
    } else {
      cparams <- NULL
    }
    if(!is.null(cparams)) {
      omega <- cparams$forest[, .(c_idx, f_idx, f_idx_uncond, wt = cvg)]
    } else {
      num_trees <- params$forest[, max(tree)]
      omega <- params$forest[, .(f_idx, f_idx_uncond = f_idx, cvg)]
      omega[, `:=` (c_idx = 1, wt = cvg / num_trees)]
      omega[, cvg := NULL]
    } 
    omega <- omega[wt > 0]
    
    if (nrow(omega) == 1) {
      omega <- omega[rep(1, n_synth),][, idx := .I]
    } else {
      if(condition_row_mode == "or") {
        omega <- merge(draws, omega, by = "f_idx", sort = FALSE)[, idx := .I]
      } else {
        omega[, N := .N, by = c_idx]
        draws <- setorder(rbind(
          omega[N != 1, .(f_idx = sample(f_idx, size = n_synth, replace = TRUE, prob = wt)), by = c_idx],
          omega[N == 1, .(f_idx, c_idx)]
        ))
        omega <- merge(draws, omega[,-"N"], by = c("c_idx", "f_idx"), sort = FALSE)[, idx := .I]
      }
      setcolorder(omega, "idx")
    }
    
    # Simulate continuous data
    synth_cnt <- synth_cat <- NULL
    if (any(!factor_cols)) {
      fam <- params$meta[family != 'multinom', unique(family)]
      if(!is.null(cparams)) {
        psi_cond <- merge(omega, cparams$cnt[,-c("cvg_factor", "f_idx_uncond")], by = c('c_idx', 'f_idx'), 
                                                sort = FALSE, allow.cartesian = TRUE)
      } else {
        psi_cond <- data.table()
      }
      psi <- unique(rbind(psi_cond,
                   merge(omega, params$cnt, by.x = 'f_idx_uncond', by.y = 'f_idx',
                         sort = FALSE, allow.cartesian = TRUE)[,val := NA_real_]
                   ), by = c("idx", "variable")
      )
      if (fam == 'truncnorm') {
        psi[is.na(val), val := truncnorm::rtruncnorm(.N, a = min, b = max, mean = mu, sd = sigma)]
      } else if (fam == 'unif') {
        psi[is.na(val), val := stats::runif(.N, min = min, max = max)]
      }
      NA_share_cnt <- psi[,.(idx, variable, NA_share)]
      synth_cnt <- dcast(psi, idx ~ variable, value.var = 'val')[, idx := NULL]
    }
    
    # Simulate categorical data
    if (any(factor_cols)) {
      if(!is.null(cparams)) {
        psi_cond <- merge(omega, cparams$cat[,-c("cvg_factor", "f_idx_uncond")], by = c('c_idx', 'f_idx'), 
                          sort = FALSE, allow.cartesian = TRUE)
      } else {
        psi_cond <- data.table()
      }
      psi <- unique(rbind(psi_cond,
                          merge(omega, params$cat, by.x = 'f_idx_uncond', by.y = 'f_idx',
                                sort = FALSE, allow.cartesian = TRUE)
      ), by = c("idx", "variable")
      )
      psi[prob < 1, val := sample(val, 1, prob = prob), by = .(variable, idx)]
      
      psi <- unique(psi[, .(idx, variable, val, NA_share)])
      NA_share_cat <- psi[,.(idx, variable, NA_share)]
      synth_cat <- dcast(psi, idx ~ variable, value.var = 'val')[, idx := NULL]
    }
    
    # Combine, optionally impose constraint(s)
    x_synth <- cbind(synth_cnt, synth_cat)
    if(length(x_synth) == 0) {
      x_synth <- condition_part[F,]
    }
    
    # Clean up, export
    x_synth <- post_x(x_synth, params)
    
    if(sample_NAs) {
      setDT(x_synth)
      NA_share <- rbind(NA_share_cnt, NA_share_cat)
      setorder(NA_share[,variable := factor(variable, levels = params$meta[,variable])], variable, idx)
      NA_share[,dat := rbinom(.N, 1, prob = NA_share)]
      x_synth[dcast(NA_share,formula =  idx ~ variable, value.var = "dat")[,-"idx"] == 1] <- NA
      
      x_synth <- post_x(x_synth, params)
    }
    
    if(condition_row_mode == "separate" & any(omega[,is.na(f_idx)])){
      setDT(x_synth)
      indices_na <- cparams$forest[is.na(f_idx), c_idx]
      indices_sampled <- cparams$forest[!is.na(f_idx), unique(c_idx)]
      rows_na <- condition_part[indices_na,]
      if(!all(factor_cols)){
        rows_na <- rows_na[,(names(x_synth)[!factor_cols]) := lapply(.SD,as.numeric),.SDcols=!factor_cols]
      }
      rows_na[, idx:= indices_na]
      rows_na <- rbindlist(replicate(n_synth, rows_na, simplify = F))
      x_synth[, idx:= rep(indices_sampled, each = n_synth)]
      x_synth <- rbind(x_synth, rows_na)
      setorder(x_synth, idx)[, idx:= NULL]
      x_synth <- post_x(x_synth, params)
    }
    x_synth
  }
  
  if(!parallel & doParRegistered & (num_workers > 1)) {
    registerDoParallel(num_workers)
  }
  
  return(x_synth_)
}


