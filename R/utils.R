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
    x[, idx_char] <- lapply(x[, idx_char, drop = FALSE], as.factor)
  }
  idx_logical <- sapply(x, is.logical)
  if (any(idx_logical)) {
    x[, idx_logical] <- lapply(x[, idx_logical, drop = FALSE], as.factor)
  }
  idx_integer <- sapply(x, is.integer)
  if (any(idx_integer)) {
    warning('Recoding integer data as ordered factors. To override this behavior, ',
            'explicitly code these variables as numeric.')
    x[, idx_integer] <- lapply(which(idx_integer), function(j) {
        lvls <- sort(unique(x[[j]]))
        factor(x[[j]], levels = lvls, ordered = TRUE)
    })
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
  if ('cnt' %in% colnames(x)) {
    colnames(x)[which(colnames(x) == 'cnt')] <- col_rename(x, 'cnt')
  }
  if ('N' %in% colnames(x)) {
    colnames(x)[which(colnames(x) == 'N')] <- col_rename(x, 'N')
  }
  return(x)
}

#' Preprocess evidence
#' 
#' This function prepares the evidence for computing leaf posteriors.
#' 
#' @param params Circuit parameters learned via \code{\link{forde}}. 
#' @param evidence Optional set of conditioning events.
#' 
#' @import data.table
#' 

prep_evi <- function(params, evidence) {
  
  # To avoid data.table check issues
  variable <- relation <- N <- n <- family <- NULL
  
  # Prep
  setDT(evidence)
  part <- all(colnames(evidence) %in% params$meta$variable)
  conj <- all(c('variable', 'relation', 'value') %in% colnames(evidence))
  post <- all(c('f_idx', 'wt') %in% colnames(evidence))
  if (part + conj + post != 1L) {
    stop('evidence must either be a partial sample, a data frame of conjuncts, ', 
         'or a posterior distribution over leaves.')
  }
  if (isTRUE(part)) {
    if (!all(colnames(evidence) %in% params$meta$variable)) {
      err <- setdiff(colnames(evidence), params$meta$variable)
      stop('Unrecognized feature(s) among colnames: ', err)
    }
    evidence <- suppressWarnings(
      melt(evidence, measure.vars = colnames(evidence), variable.factor = FALSE)
    )
    evidence[, relation := '==']
    conj <- TRUE
  }
  if (isTRUE(conj)) {
    evi <- merge(params$meta, evidence, by = 'variable', sort = FALSE)
    evi[, n := .N, by = variable]
    if (any(evi[n > 1L, relation == '=='])) {
      culprit <- evi[n > 1L & relation == '==', variable]
      stop(paste('Inconsistent conditioning events for the following variable(s):', 
                 culprit))
    }
    if (any(evi[, family == 'multinom'])) {
      evi_tmp <- evi[family == 'multinom']
      if (any(evi_tmp[, !relation %in% c('==', '!=')])) {
        stop('With categorical features, the only valid relations are ',
             '"==" or "!=".')
      }
    }
    if (any(evi[, family != 'multinom'])) {
      evi_tmp <- evi[family != 'multinom']
      if (any(evi_tmp[, relation == '!='])) {
        evidence <- evidence[!(family != 'multinom' & relation == '!=')]
        warning('With continuous features, "!=" is not a valid relation. ', 
                'This constraint has been removed.')
      }
      #if (any(evi_tmp[, n > 2L])) {
      #  inf <- blah
      #  sup <- blah
      #}
    }
  }
  
  ### ALSO: Reduce redundant events to most informative condition
  ###       and check for inconsistencies, e.g. >3 & <2
  
  
  return(evidence)
}


#' Compute leaf posterior
#' 
#' This function returns a posterior distribution on leaves, conditional on some 
#' evidence. 
#' 
#' @param params Circuit parameters learned via \code{\link{forde}}.
#' @param evidence Data frame of conditioning event(s).
#' 
#' @import data.table
#' @importFrom truncnorm dtruncnorm ptruncnorm 
#' 

leaf_posterior <- function(params, evidence) {
  
  # To avoid data.table check issues
  variable <- relation <- value <- prob <- f_idx <- cvg <- wt <- 
    mu <- sigma <- val <- k <- family <- n <- compl <- lik2 <- . <- NULL
  
  # Likelihood per leaf-event combo
  psi_cnt <- psi_cat <- NULL
  evidence <- merge(evidence, params$meta, by = 'variable', sort = FALSE)
  
  # Continuous features
  if (any(evidence$family != 'multinom')) { 
    evi <- evidence[family != 'multinom']
    psi <- merge(evi, params$cnt, by = 'variable')
    evi[, n := .N, by = variable]
    if (any(evi$relation == '==')) {
      psi[relation == '==', lik := 
            truncnorm::dtruncnorm(value, a = min, b = max, mean = mu, sd = sigma)]
    }
    if (any(evi$relation != '==')) {
      psi[relation != '==', lik := 
            truncnorm::ptruncnorm(value, a = min, b = max, mean = mu, sd = sigma)]
    }
    if (any(evi[, n > 1L])) {
      interval_vars <- evi[n > 1L, variable]
      psi1 <- psi[!variable %in% interval_vars]
      psi2 <- psi[variable %in% interval_vars]
      psi2[, lik := rep(psi2[relation %in% c('<', '<='), lik] - 
                          psi2[relation %in% c('>', '>='), lik], 2)]
      psi <- rbind(psi1, psi2)
    }
    if (any(evi[, n == 1 & relation %in% c('>', '>=')])) {
      psi[n == 1 & relation %in% c('>', '>='), lik := 1 - lik]
    }
    psi[value == min, lik := 0]
    psi_cnt <- unique(psi[, .(f_idx, variable, lik)])
  }
  
  # Categorical features
  psi_eq <- psi_ineq <- NULL
  if (any(evidence$family == 'multinom')) { 
    evi <- evidence[family == 'multinom']
    grd <- rbindlist(lapply(evi[, variable], function(j) {
      expand.grid('f_idx' = params$forest$f_idx, 'variable' = j,
                  'val' = params$cat[variable == j, unique(val)],
                  stringsAsFactors = FALSE)
    }))
    psi <- merge(params$cat, grd, by = c('f_idx', 'variable', 'val'),
                 sort = FALSE, all.y = TRUE)
    psi[is.na(prob), prob := 0]
    setnames(psi, 'prob', 'lik')
    if (any(evi[, relation == '=='])) {
      evi_tmp <- evi[relation == '==', .(variable, value)]
      setnames(evi_tmp, 'value', 'val')
      psi_eq <- merge(psi, evi_tmp, by = c('variable', 'val'), sort = FALSE)
      psi_eq <- psi_eq[, .(f_idx, variable, lik)]
    }
    if (any(evi[, relation == '!='])) {
      evi_tmp <- evi[relation == '!=', .(variable, value)]
      psi_ineq <- rbindlist(lapply(evi_tmp[, .I], function(k) {
        psi[variable == evi_tmp$variable[k] & val != evi_tmp$value[k]]
      }))
      psi_ineq[, lik := sum(lik), by = .(f_idx, variable)]
      psi_ineq <- unique(psi_ineq[, .(f_idx, variable, lik)])
    }
    psi_cat <- rbind(psi_eq, psi_ineq)
  }
  psi <- rbind(psi_cnt, psi_cat)
  
  # Weight is proportional to coverage times product of likelihoods
  psi <- merge(psi, params$forest[, .(f_idx, cvg)], by = 'f_idx', sort = FALSE)
  psi[, wt := cvg * prod(lik), by = f_idx] 
  
  # Normalize, export
  out <- unique(psi[wt > 0, .(f_idx, wt)])
  out[, wt := wt / sum(wt)]
  return(out)
}


#' Post-process data
#' 
#' This function prepares output data for forge.
#' 
#' @param x Input data.frame.
#' @param params Circuit parameters learned via \code{\link{forde}}.
#' 
#' @import data.table
#' @importFrom tibble as_tibble
#'

post_x <- function(x, params) {
  
  # To avoid data.table check issues
  variable <- NULL
  
  # Order, classify features
  meta_tmp <- params$meta[variable %in% colnames(x)]
  setcolorder(x, match(meta_tmp$variable, colnames(x)))
  setDF(x)
  idx_numeric <- meta_tmp[, which(class == 'numeric')]
  idx_factor <- meta_tmp[, which(class == 'factor')]
  idx_ordered <- meta_tmp[, which(grepl('ordered', class))]
  idx_logical <- meta_tmp[, which(class == 'logical')]
  idx_integer <- meta_tmp[, which(class == 'integer')]
  
  # Recode
  if (sum(idx_numeric) > 0L) {
    x[, idx_numeric] <- lapply(idx_numeric, function(j) {
        round(as.numeric(x[[j]]), meta_tmp$precision[j])
    })
  }
  if (sum(idx_factor) > 0L) {
    x[, idx_factor] <- lapply(idx_factor, function(j) {
        factor(x[[j]], levels = params$levels[[colnames(x)[j]]])
    })
  }
  if (sum(idx_ordered) > 0L) {
    x[, idx_ordered] <- lapply(idx_ordered, function(j) {
        factor(x[[j]], levels = params$levels[[colnames(x)[j]]], ordered = TRUE)
    })
  }
  if (sum(idx_logical) > 0L) {
    x[, idx_logical] <- lapply(x[, idx_logical, drop = FALSE], as.logical)
  }
  if (sum(idx_integer) > 0L) {
    x[, idx_integer] <- lapply(idx_integer, function(j) {
      as.integer(as.character(x[[j]]))
    }) 
  }
  
  # Export
  if ('data.table' %in% params$input_class) {
    setDT(x)
  } else if ('tbl_df' %in% params$input_class) {
    x <- tibble::as_tibble(x)
  } else if ('matrix' %in% params$input_class) {
    x <- as.matrix(x)
  }
  return(x)
}


