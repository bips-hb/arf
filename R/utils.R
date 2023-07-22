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

#' Preprocess evidence
#' 
#' This function prepares the evidence for computing leaf posteriors.
#' 
#' @param evidence 
#' 
#' @import data.table
#' 

prep_evi <- function(pc, evidence) {
  
  # To avoid data.table check issues
  relation <- N <- NULL
  
  # Prep
  evidence <- as.data.table(evidence)
  part <- all(colnames(evidence) %in% pc$meta$variable)
  conj <- all(c('variable', 'relation', 'value') %in% colnames(evidence))
  post <- all(c('f_idx', 'wt') %in% colnames(evidence))
  if (part + conj + post != 1L) {
    stop('evidence must either be a partial sample, a data frame of conjuncts, ', 
         'or a posterior distribution over leaves.')
  }
  if (part) {
    if (!all(colnames(evidence) %in% pc$meta$variable)) {
      err <- setdiff(colnames(evidence), pc$meta$variable)
      stop('Unrecognized feature(s) among colnames: ', err)
    }
    evidence <- suppressWarnings(melt(evidence, measure.vars = colnames(evidence)))
    evidence[, relation := '==']
    conj <- TRUE
  }
  if (conj) {
    if (max(evidence[, .N, by = variable]$N > 1L)) {
      stop('Only one constraint per variable allowed when using conjuncts.')
    }
    evi <- merge(pc$meta, evidence, by = 'variable', sort = FALSE)
    if (evi[class == 'numeric' & relation == '!=', .N] > 0) {
      evidence <- evidence[!(class == 'numeric' & relation == '!=')]
      warning('With continuous features, "!=" is not a valid relation. ', 
              'This constraint has been removed.')
    }
    if (evi[class != 'numeric' & !relation %in% c('==', '!='), .N] > 0) {
      stop('With categorical features, the only valid relations are ',
           '"==" or "!=".')
    }
  }
  return(evidence)
}


#' Compute leaf posterior
#' 
#' This function returns a posterior distribution on leaves, conditional on some 
#' evidence. 
#' 
#' @param params Parameters learned via \code{\link{forde}}. 
#' @param evidence Data frame of conditioning event(s).
#' @param parallel Compute in parallel?
#' 
#' @import data.table
#' @importFrom truncnorm dtruncnorm ptruncnorm 
#' 

leaf_posterior <- function(params, evidence, parallel) {
  
  # To avoid data.table check issues
  variable <- relation <- value <- prob <- f_idx <- cvg <- wt <- 
    mu <- sigma <- val <- k <- . <- NULL
  
  # Likelihood per leaf-event combo
  psi_cnt <- psi_cat <- NULL
  evidence <- merge(evidence, params$meta, by = 'variable', sort = FALSE)
  if (any(evidence$class == 'numeric')) { # Continuous features
    evi <- evidence[class == 'numeric']
    psi <- merge(evi, params$cnt, by = 'variable')
    if (any(evi$relation == '==')) {
      psi[relation == '==', prob := 
            truncnorm::dtruncnorm(value, a = min, b = max, mean = mu, sd = sigma)]
    }
    if (any(evi$relation %in% c('<', '<='))) {
      psi[relation %in% c('<', '<='), prob := 
            truncnorm::ptruncnorm(value, a = min, b = max, mean = mu, sd = sigma)]
    }
    if (any(evi$relation %in% c('>', '>='))) {
      psi[relation %in% c('>', '>='), prob := 
            1 - truncnorm::ptruncnorm(value, a = min, b = max, mean = mu, sd = sigma)]
    }
    psi[value == min, prob := 0]
    psi_cnt <- psi[, .(f_idx, variable, prob)]
  }
  if (any(evidence$class != 'numeric')) { # Categorical features
    evi <- evidence[class != 'numeric']
    cat_upd8 <- function(k) {
      j <- evi$variable[k]
      op <- evi$relation[k]
      value <- evi$value[k]
      psi <- params$cat[variable == j]
      grd <- expand.grid(f_idx = params$forest$f_idx, val = psi[, unique(val)])
      psi <- merge(psi, grd, by = c('f_idx', 'val'), all.y = TRUE, sort = FALSE)
      psi[is.na(prob), prob := 0][is.na(variable), variable := j]
      if (op == '==') {
        out <- psi[val == value, .(f_idx, variable, prob)]
      } else if (op == '!=') {
        psi <- psi[val != value]
        psi[, prob := sum(prob), by = f_idx]
        out <- unique(psi[, .(f_idx, variable, prob)])
      }
      return(out)
    }
    if (isTRUE(parallel)) {
      psi_cat <- foreach(k = seq_len(evi[, .N]), .combine = rbind) %dopar% cat_upd8(k)
    } else {
      psi_cat <- foreach(k = seq_len(evi[, .N]), .combine = rbind) %do% cat_upd8(k)
    }
  }
  psi <- rbind(psi_cnt, psi_cat)
  
  # Weight is proportional to coverage times product of likelihoods
  psi <- merge(psi, params$forest[, .(f_idx, cvg)], by = 'f_idx', sort = FALSE)
  psi[, wt := cvg * prod(prob), by = f_idx] # Worth doing in log space?
  
  # Normalize, export
  out <- unique(psi[, .(f_idx, wt)])
  out[, wt := wt / sum(wt)]
  return(out)
}


