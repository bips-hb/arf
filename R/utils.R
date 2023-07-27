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

#' Preprocess input data
#' 
#' This function prepares input data for ARFs.
#' 
#' @param x Input data.frame.
#' 

prep_x <- function(x) {
  # Reclass all non-numeric features as factors
  x <- as.data.frame(x)
  idx_char <- sapply(x, is.character)
  if (any(idx_char)) {
    x[, idx_char] <- setDF(
      lapply(x[, idx_char, drop = FALSE], as.factor)
    )
  }
  idx_logical <- sapply(x, is.logical)
  if (any(idx_logical)) {
    x[, idx_logical] <- setDF(
      lapply(x[, idx_logical, drop = FALSE], as.factor)
    )
  }
  idx_integer <- sapply(x, is.integer)
  if (any(idx_integer)) {
    warning('Recoding integer data as ordered factors. To override this behavior, ',
            'explicitly code these variables as numeric.')
    for (j in which(idx_integer)) {
      lvls <- sort(unique(x[, j]))
      x[, j] <- factor(x[, j], levels = lvls, ordered = TRUE)
    }
  }
  # Rename annoying columns
  if ('y' %in% colnames(x)) {
    colnames(x)[which(colnames(x) == 'y')] <- col_rename(x, 'y')
  }
  if ('obs' %in% colnames(x)) {
    colnames(x)[which(colnames(x) == 'obs')] <- col_rename(x, 'obs')
  }
  if ('tree' %in% colnames(x)) {
    colnames(x)[which(colnames(x) == 'tree')] <- col_rename(x, 'tree')
  } 
  if ('leaf' %in% colnames(x)) {
    colnames(x)[which(colnames(x) == 'leaf')] <- col_rename(x, 'leaf')
  } 
  return(x)
}

#' Preprocess evidence
#' 
#' This function prepares the evidence for computing leaf posteriors.
#' 
#' @param pc Probabilistic circuit learned via \code{\link{forde}}. 
#' @param evidence Optional set of conditioning events.
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
  if (isTRUE(part)) {
    if (!all(colnames(evidence) %in% pc$meta$variable)) {
      err <- setdiff(colnames(evidence), pc$meta$variable)
      stop('Unrecognized feature(s) among colnames: ', err)
    }
    evidence <- suppressWarnings(
      melt(evidence, measure.vars = colnames(evidence))
    )
    evidence[, relation := '==']
    conj <- TRUE
  }
  if (isTRUE(conj)) {
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
#' @param pc Probabilistic circuit learned via \code{\link{forde}}.
#' @param evidence Data frame of conditioning event(s).
#' @param parallel Compute in parallel?
#' 
#' @import data.table
#' @importFrom truncnorm dtruncnorm ptruncnorm 
#' 

leaf_posterior <- function(pc, evidence, parallel) {
  
  # To avoid data.table check issues
  variable <- relation <- value <- prob <- f_idx <- cvg <- wt <- 
    mu <- sigma <- val <- k <- . <- NULL
  
  # Likelihood per leaf-event combo
  psi_cnt <- psi_cat <- NULL
  evidence <- merge(evidence, pc$meta, by = 'variable', sort = FALSE)
  if (any(evidence$class == 'numeric')) { # Continuous features
    evi <- evidence[class == 'numeric']
    psi <- merge(evi, pc$cnt, by = 'variable')
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
    cat_upd8 <- function(k) {             # Redo as rbindlist?
      j <- evi$variable[k]
      op <- evi$relation[k]
      value <- evi$value[k]
      psi <- pc$cat[variable == j]
      grd <- expand.grid(f_idx = pc$forest$f_idx, val = psi[, unique(val)])
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
  psi <- merge(psi, pc$forest[, .(f_idx, cvg)], by = 'f_idx', sort = FALSE)
  psi[, wt := cvg * prod(prob), by = f_idx] # Worth doing in log space?
  
  # Normalize, export
  out <- unique(psi[, .(f_idx, wt)])
  out[, wt := wt / sum(wt)]
  return(out)
}


#' Post-process data
#' 
#' This function prepares output data for map and forge.
#' 
#' @param x Input data.frame.
#' @param pc Probabilistic circuit learned via \code{\link{forde}}.
#' 

post_x <- function(x, pc) {
  
  # Order, classify features
  meta_tmp <- pc$meta[variable %in% colnames(x)]
  setcolorder(x, match(colnames(x), meta_tmp$variable))
  setDF(x)
  idx_factor <- meta_tmp[, which(class == 'factor')]
  idx_logical <- meta_tmp[, which(class == 'logical')]
  idx_integer <- meta_tmp[, which(class == 'integer')]
  
  # Recode
  if (sum(idx_character) > 0L) {
    x[, idx_character] <- setDF(
      lapply(x[, idx_character, drop = FALSE], function(j) as.factor(j))
    )
  }
  if (sum(idx_logical) > 0L) {
    x[, idx_logical] <- setDF(
      lapply(x[, idx_logical, drop = FALSE], function(j) as.logical(j))
    )
  }
  if (sum(idx_integer) > 0L) {
    x[, idx_integer] <- setDF(
      lapply(x[, idx_integer, drop = FALSE], function(j) as.integer(j)) 
    )
  }
  
  # Export
  if ('data.table' %in% pc$input_class) {
    setDT(x)
  } else if ('tbl_df' %in% pc$input_class) {
    x <- tibble::as_tibble(x)
  } else if ('matrix' %in% pc$input_class) {
    x <- as.matrix(x)
  }
  return(x)
}


