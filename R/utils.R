#' Adaptive column renaming
#' 
#' This function renames columns in case the input data.frame includes any
#' colnames required by internal functions (e.g., \code{"y"}).
#' 
#' @param df Input data.frame.
#' @param old_name Name of column to be renamed.
#' 

col_rename <- function(df, old_name) {
  k <- 1L
  converged <- FALSE
  while (!isTRUE(converged)) {
    new_name <- paste0(old_name, k)
    if (!new_name %in% colnames(df)) {
      converged <- TRUE
    } else {
      k <- k + 1L
    }
  }
  return(new_name)
}

#' Compute leaf posterior
#' 
#' This function returns a posterior distribution on leaves, conditional on some 
#' evidence. 
#' 
#' @param params Parameters learned via \code{\link{forde}}. 
#' @param evidence Data frame of conditioning event(s).
#' 
#' @import data.table
#' @importFrom truncnorm dtruncnorm ptruncnorm 
#' 

leaf_posterior <- function(params, evidence) {
  
  # To avoid data.table check issues
  variable <- operator <- value <- prob <- f_idx <- cvg <- wt <- . <- NULL
  
  # Likelihood per leaf-event
  leaf_list <- lapply(seq_len(nrow(evidence)), function(k) {
    j <- evidence$variable[k]
    op <- evidence$operator[k]
    value <- evidence$value[k]
    if (params$meta[variable == j, class == 'numeric']) {
      value <- as.numeric(value)
      psi <- params$cnt[variable == j]
      if (op == '==') {
        psi[, prob := truncnorm::dtruncnorm(value, a = min, b = max, mean = mu, sd = sigma)]
      } else if (op %in% c('<', '<=')) {
        psi[, prob := truncnorm::ptruncnorm(value, a = min, b = max, mean = mu, sd = sigma)]
      } else if (op %in% c('>', '>=')) {
        psi[, prob := 1 - truncnorm::ptruncnorm(value, a = min, b = max, mean = mu, sd = sigma)]
      }
    } else {
      psi <- params$cat[variable == j]
      grd <- expand.grid(f_idx = params$forest$f_idx, val = psi[, unique(val)])
      psi <- merge(psi, grd, by = c('f_idx', 'val'), all.y = TRUE, sort = FALSE)
      psi[is.na(prob), prob := 0][is.na(variable), variable := j]
      if (op == '==') {
        psi <- psi[val == value]
      } else if (op == '!=') {
        psi <- psi[val != value]
        psi[, prob := sum(prob), by = f_idx]
        psi <- unique(psi[, .(f_idx, variable, prob)])
      }
    }
    out <- psi[, .(f_idx, variable, prob)]
    return(out)
  })
  
  # Weight is proportional to coverage times product of likelihoods
  psi <- rbindlist(leaf_list)
  psi <- merge(psi, params$forest[, .(f_idx, cvg)], by = 'f_idx', sort = FALSE)
  psi[, wt := cvg * prod(prob), by = f_idx] # Worth doing in log space?
  
  # Normalize, export
  out <- unique(psi[, .(f_idx, wt)])
  out[, wt := wt / sum(wt)]
  return(out)
}


