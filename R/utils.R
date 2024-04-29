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
    # Recoding integers with > 5 levels as numeric
    to_numeric <- sapply(seq_len(ncol(x)), function(j) {
      idx_integer[j] & length(unique(x[[j]])) > 5
    })
    if (any(to_numeric)) {
      warning('Recoding integers with more than 5 unique values as numeric. ', 
              'To override this behavior, explicitly code these variables as factors.')
      x[, to_numeric] <- lapply(x[, to_numeric, drop = FALSE], as.numeric)
    }
    to_factor <- sapply(seq_len(ncol(x)), function(j) {
      idx_integer[j] & length(unique(x[[j]])) < 6
    })
    if (any(to_factor)) {
      warning('Recoding integers with fewer than 6 unique values as ordered factors. ', 
              'To override this behavior, explicitly code these variables as numeric.')
      x[, to_factor] <- lapply(which(to_factor), function(j) {
        lvls <- sort(unique(x[[j]]))
        factor(x[[j]], levels = lvls, ordered = TRUE)
      })
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
  variable <- relation <- N <- n <- family <- wt <- NULL
  
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
  if (isTRUE(post)) {
    if (evidence[, sum(wt)] != 1) {
      evidence[, wt := wt / sum(wt)]
      warning('Posterior weights have been normalized to sum to unity.')
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
    evi[, n := .N, by = variable]
    psi <- merge(evi, params$cnt, by = 'variable')
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
    evi[, value := as.character(value)]
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
  psi[, wt:= {
    if (any(lik == 0)) {
      0
    } else {
      exp(mean(log(c(cvg[1], lik))))
    }
  }
  , by = f_idx]
  
  # Normalize, export
  out <- unique(psi[wt > 0, .(f_idx, wt)])
  if (nrow(out) == 0) {
    # If all leaves have zero weight, choose one randomly
    warning("All leaves have zero likelihood. This is probably because evidence contains an (almost) impossible combination. For categorical data, consider setting alpha>0 in forde().")
    out <- unique(psi[, .(f_idx)])
    out[, wt := 1]
  }
  out[, wt := (wt / max(wt, na.rm = T))^(nrow(evidence) + 1)][wt > 0, wt := wt / sum(wt)]
  return(out[])
}


#' Post-process data
#' 
#' This function prepares output data for forge.
#' 
#' @param x Input data.frame.
#' @param params Circuit parameters learned via \code{\link{forde}}.
#' 
#' @import data.table
#'

post_x <- function(x, params) {
  
  # To avoid data.table check issues
  variable <- val <- NULL
  
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
        round(as.numeric(x[[j]]), meta_tmp$decimals[j])
    })
  }
  if (sum(idx_factor) > 0L) {
    x[, idx_factor] <- lapply(idx_factor, function(j) {
      factor(x[[j]], levels = params$levels[variable == colnames(x)[j], val])
    })
  }
  if (sum(idx_ordered) > 0L) {
    x[, idx_ordered] <- lapply(idx_ordered, function(j) {
        factor(x[[j]], levels = params$levels[variable == colnames(x)[j], val], ordered = TRUE)
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
    setDT(x)[]
  } else if ('tbl_df' %in% params$input_class & requireNamespace("tibble", quietly = TRUE)) {
    x <- tibble::as_tibble(x)
  } else if ('matrix' %in% params$input_class) {
    x <- as.matrix(x)
  }
  return(x)
}


#' Compute conditional circuit parameters
#' 
#' @param params Circuit parameters learned via \code{\link{forde}}.
#' @param condition Data frame of conditioning event(s).
#' @param row_mode Interpretation of rows in multi-row conditions.
#' @param stepsize Stepsize defining number of condition rows handled in one for each step.
#' 
#' @import data.table
#' @importFrom foreach foreach %dopar% getDoParRegistered getDoParWorkers
#' @importFrom truncnorm dtruncnorm ptruncnorm 
#' @importFrom stats dunif punif
#' 

cforde <- function(params, condition, row_mode = c("separate", "or"), stepsize = 0, parallel = TRUE) {
  
  row_mode <- match.arg(row_mode)
  
  doParRegistered <- getDoParRegistered()
  num_workers <- getDoParWorkers()
  if(!parallel & doParRegistered & (num_workers > 1)) {
    registerDoSEQ()
  }
  
  # To avoid data.table check issues
  . <- c_idx <- cvg <- cvg_arf <- cvg_factor <- f_idx <- f_idx_uncond <- i.max <-
    i.min <- leaf <- max.x <- max.y <- min.x <- min.y <- mu <- prob <- sigma <-
    step_ <-	tree <-	V1 <- val <- variable <- NULL
  
  meta <- params$meta
  family <- meta[family != "multinom", unique(family)]
  forest <- params$forest
  cat <- params$cat
  cnt <- params$cnt
  cnt_cols <-meta[family != "multinom", variable]
  cat_cols <-meta[family == "multinom", variable]
  
  # format c, calculate DNF and output disjoint hyperrectangles
  condition_long <- prep_cond(condition,params, row_mode)
  setkey(condition_long, c_idx)
  
  if(nrow(condition_long) == 0){
    return(NULL)
  }
  
  # store resulting number of disjoint hyperrectangles
  
  if(row_mode == "or") {
    nconds <- nconds_conditioned <- condition_long[,max(c_idx)]
  } else {
    nconds <- nrow(condition)
    nconds_conditioned <- condition_long[,uniqueN(c_idx)]
  }
  
  conds_conditioned <- condition_long[,unique(c_idx)]
  
  if(stepsize == 0) {
    if(parallel) {
      stepsize <- ceiling(nconds_conditioned/getDoParWorkers())
    } else {
      stepsize <- nconds_conditioned
    }
  }
  step_no <- ceiling(nconds_conditioned/stepsize)
  updates_relevant_leaves <- foreach(step_ = 1:step_no, .combine = "rbind") %dopar% {
    
    index_start <- conds_conditioned[(step_ - 1)*stepsize + 1]
    index_end <- conds_conditioned[min(step_ * stepsize, nconds_conditioned)]
    condition_long_step <- condition_long[.(index_start:index_end),nomatch = NULL]
    
    # store conditions for cat and cnt separately
    cat_conds <- condition_long_step[variable %in% cat_cols,c("c_idx","variable","val")][,variable := factor(variable)]
    cnt_conds <- condition_long_step[variable %in% cnt_cols,c("c_idx","variable","min", "max","val")][,`:=` (variable = factor(variable),
                                                                                                         val = as.numeric(val))]
    
    if (nrow(cat_conds) != 0) {
      cat_relevant <- cat[,.(.(f_idx)), by=.(variable,val)]
      setkey(cat_relevant, variable, val)
      setkey(cat_conds, variable, val)
      cat_relevant <- cat_conds[cat_relevant, on = .(variable, val), nomatch = NULL]
      setkey(cat_relevant, c_idx)
      relevant_leaves_changed_cat <- cat_relevant[, Reduce(intersect,V1),by = c_idx][,.(c_idx, f_idx = V1)]
      setorder(relevant_leaves_changed_cat)
      
      conditions_unchanged_cat <- setdiff(condition_long_step[, c_idx], cat_conds[, c_idx])
      relevant_leaves_unchanged_cat <- data.table(c_idx = rep(conditions_unchanged_cat, each = nrow(forest) ), f_idx = rep(forest[,f_idx],length(conditions_unchanged_cat)))
      relevant_leaves_cat <- rbind(relevant_leaves_changed_cat, relevant_leaves_unchanged_cat, fill = T)
      relevant_leaves_cat_list <- relevant_leaves_cat[,.(f_idx = .(f_idx)),by=c_idx]
    } else {
      relevant_leaves_cat <- data.table(c_idx = integer(), f_idx = integer())
    }
    
    
    if (nrow(cnt_conds) != 0) {
      cnt_conds_compact <- copy(cnt_conds)
      cnt_conds_compact[!is.na(val), `:=`(min = val, max = val)][,val := NULL]
      
      cnt_relevant <- cnt[,.(min = .(min), max = .(max)),by = variable]
      cnt_relevant <- cnt_conds_compact[cnt_relevant, on = .(variable), nomatch = NULL]
      setkey(cnt_relevant,c_idx)

      if (nrow(cat_conds) != 0) {
        cnt_relevant <- cnt_relevant[relevant_leaves_cat_list, on = .(c_idx)]
      } else {
        cnt_relevant[, f_idx := NA]
      }

      relevant_leaves_changed_cnt <- cnt_relevant[, .(
        c_idx,
        variable,
        f_idx = Map(\(f_idx, min, max, i.min, i.max) {
          if (!inherits(f_idx, "logical")) {
            rel_cat_min <- i.min[f_idx]
            rel_cat_max <- i.max[f_idx]
            rel_min <- f_idx[which(max > rel_cat_min)]
            rel_max <- f_idx[which(min <= rel_cat_max)]
          } else {
            rel_cat_min <- i.min
            rel_cat_max <- i.max
            rel_min <- which(max > rel_cat_min)
            rel_max <- which(min <= rel_cat_max)
          }
          intersect(rel_min,rel_max) 
        }, f_idx = f_idx, min = min, max = max, i.min = i.min, i.max = i.max))
        ][, Reduce(intersect,f_idx),by = c_idx][,.(c_idx, f_idx = V1)]
      
      conditions_unchanged_cnt <- setdiff(condition_long_step[, c_idx], cnt_conds[, c_idx])
      relevant_leaves_unchanged_cnt <- data.table(c_idx = rep(conditions_unchanged_cnt, each = nrow(forest)), f_idx = rep(forest[,f_idx],length(conditions_unchanged_cnt)))
      relevant_leaves_cnt <- rbind(relevant_leaves_changed_cnt, relevant_leaves_unchanged_cnt)
      
      cnt_new <- merge(merge(relevant_leaves_cnt, cnt_conds, by = "c_idx", allow.cartesian = T, sort = F), cnt, by = c("f_idx", "variable"), all.x = T, allow.cartesian = T, sort = F)
      cnt_new[!is.na(val),`:=` (min = min.y,
                                max = max.y)]
      cnt_new[is.na(val),`:=` (min = pmax(min.x, min.y, na.rm = T),
                               max = pmin(max.x, max.y, na.rm = T))]
      cnt_new[,cvg_factor := NA_real_]
      if (family == "truncnorm") {
        cnt_new[!is.na(val), cvg_factor := dtruncnorm(val, a=min.y, b=max.y, mean=mu, sd=sigma)*(val != min.y)]
        cnt_new[is.na(val) & (min == min.y) & (max == max.y), cvg_factor := 1]
        cnt_new[is.na(val) & is.na(cvg_factor), cvg_factor := ptruncnorm(max, a=min.y, b=max.y, mean=mu, sd=sigma) - ptruncnorm(min, a=min.y, max.y, mean=mu,sd=sigma)] 
      } else if (family == "unif") {
        cnt_new[!is.na(val), cvg_factor := dunif(val, min=min.y, max=max.y)*(val != min.y)]
        cnt_new[is.na(val) & (min == min.y) & (max == max.y), cvg_factor := 1]
        cnt_new[is.na(val) & is.na(cvg_factor), cvg_factor := punif(max, min=min.y, max=max.y) - punif(min, min=min.y, max.y)]  
      }
      cnt_new[,c("min.x","max.x","min.y","max.y") := NULL]
      if (nrow(cat_conds) > 0) {
          relevant_leaves <- merge(relevant_leaves_cnt, relevant_leaves_cat, by = c("c_idx", "f_idx"))[,.(c_idx, f_idx)] #necessary?
      } else {
        relevant_leaves <- relevant_leaves_cnt[,.(c_idx, f_idx)]
      }
    } else {
      relevant_leaves <- relevant_leaves_cat[,.(c_idx, f_idx)]
      cnt_new <- cbind(cnt[F,], data.table(cvg_factor = numeric(), c_idx = integer(), val = numeric()))
    }
    cat_new <- merge(merge(relevant_leaves, cat_conds, by = "c_idx", allow.cartesian = T), cat, by = c("f_idx","variable", "val")) 
    cat_new[,`:=` (cvg_factor = prob, prob = 1)]
    
    list(cnt_new = cnt_new, cat_new = cat_new, relevant_leaves = relevant_leaves)
  }
  
  if(is.matrix(updates_relevant_leaves)) {
    updates_relevant_leaves <- lapply(as.data.table(updates_relevant_leaves), rbindlist) 
  }

  relevant_leaves <- updates_relevant_leaves$relevant_leaves[,`:=` (f_idx = .I, f_idx_uncond = f_idx)][]
  
  cnt_new <- setcolorder(merge(relevant_leaves, updates_relevant_leaves$cnt_new, by.x = c("c_idx", "f_idx_uncond"), by.y = c("c_idx", "f_idx"), sort = F), c("f_idx","c_idx","variable","min","max","val","cvg_factor"))[]
  cat_new <- setcolorder(merge(relevant_leaves, updates_relevant_leaves$cat_new, by.x = c("c_idx", "f_idx_uncond"), by.y = c("c_idx", "f_idx"), sort = F), c("f_idx","c_idx","variable","val","prob","cvg_factor"))[]
    
  if(relevant_leaves[,uniqueN(c_idx)] < nconds_conditioned) {
    if(relevant_leaves[,uniqueN(c_idx)] == 0 & row_mode == "or") {
      stop("For all entered evidence rows, no matching leaves could be found. This is probably because evidence lies outside of the distribution calculated by FORDE. For continuous data, consider setting epsilon>0 or finite_bounds=FALSE in forde(). For categorical data, consider setting alpha>0 in forde()")
    } else {
      warning("For some entered evidence rows, no matching leaves could be found. This is probably because evidence lies outside of the distribution calculated by FORDE. For continuous data, consider setting epsilon>0 or finite_bounds=FALSE in forde(). For categorical data, consider setting alpha>0 in forde()")
      conds_impossible <- conds_conditioned[!(conds_conditioned %in% relevant_leaves[,unique(c_idx)])]
      relevant_leaves <- setorder(rbind(relevant_leaves, data.table(c_idx = conds_impossible, f_idx = NA_integer_, f_idx_uncond = NA_integer_)))
    }
  }
  
  forest_new <- merge(relevant_leaves, forest, by.x = "f_idx_uncond", by.y = "f_idx", all.x = T, sort = F)
  setnames(forest_new,"cvg","cvg_arf")
  
  cvg_new <- rbind(cat_new[,.(f_idx, c_idx, cvg_factor)], cnt_new[, .(f_idx, c_idx, cvg_factor)])

  if(nrow(cvg_new) > 0) {
    cvg_new[,cvg_factor := log(cvg_factor)]
    cvg_new <- cvg_new[, .(cvg_factor = sum(cvg_factor)), keyby = f_idx]
    cvg_new <- cbind(cvg_new, forest_new[, .(c_idx, cvg_arf = log(cvg_arf))])
    cvg_new[,`:=` (cvg = cvg_factor + cvg_arf, cvg_factor = NULL, cvg_arf = NULL)]
    
    if(row_mode == "or") {
      if(cvg_new[,all(cvg == - Inf)]) {
        warning("All leaves have zero likelihood. This is probably because evidence contains an (almost) impossible combination.")
        cvg_new[, cvg := 1/.N]
      } else {
        cvg_new[, cvg := exp(cvg - max(cvg))]
        cvg_new <- cvg_new[cvg > 0,][,cvg := cvg / sum(cvg)]
      }
    } else {
      cvg_new[, leaf_zero_lik := all(cvg == -Inf), by = c_idx]
      if(any(cvg_new[, leaf_zero_lik])) {
        warning("All leaves have zero likelihood for some entered evidence rows. This is probably because evidence contains an (almost) impossible combination.")
        cvg_new[leaf_zero_lik == T, cvg := 1/.N, by = c_idx]
      }
      cvg_new[leaf_zero_lik == F, scale := max(cvg), by = c_idx]
      cvg_new[leaf_zero_lik == F, cvg := exp(cvg - scale)]
      cvg_new[leaf_zero_lik == F, scale := sum(cvg), by = c_idx]
      cvg_new[leaf_zero_lik == F, cvg := cvg / scale]
      cvg_new[, `:=` (leaf_zero_lik = NULL, scale = NULL)]
    }
  }
  
  forest_new_noleaf <- data.table(c_idx = setdiff(unique(forest_new[,c_idx]), unique(cvg_new[,c_idx])))[,f_idx := NA_integer_]
  forest_new <- cbind(forest_new, cvg_new[,.(cvg)])
  forest_new <- forest_new[cvg > 0,]
  forest_new <- rbind(forest_new, forest_new_noleaf,fill = T)
  if (row_mode == "or") {
    if(forest_new[,all(is.na(f_idx))]) {
      forest_new[is.na(f_idx), cvg := 1/.N]
    } else {
      forest_new[is.na(f_idx), cvg := 0]
    }
  } else {
    forest_new[is.na(f_idx), cvg := 1]
  }
  if(row_mode == "separate" & (nconds != nconds_conditioned)) {
    conds_unconditioned <- (1:nconds)[!(1:nconds) %in% conds_conditioned]
    forest_new_unconditioned <- copy(forest)
    forest_new_unconditioned <- rbindlist(replicate(length(conds_unconditioned), forest, simplify = F))
    forest_new_unconditioned[, `:=` (c_idx = rep(conds_unconditioned,each = nrow(forest)), f_idx_uncond = f_idx, cvg_arf = cvg)]
    forest_new <- rbind(forest_new, forest_new_unconditioned)
  }

  setorder(setcolorder(forest_new,c("f_idx","c_idx","f_idx_uncond","tree","leaf","cvg_arf","cvg")), c_idx, f_idx, f_idx_uncond, tree, leaf)
  
  if(!parallel & doParRegistered & (num_workers > 1)) {
    registerDoParallel(num_workers)
  }
  
  list(condition_input = condition, condition_prepped = condition_long, cnt = cnt_new, cat = cat_new, forest = forest_new)
}


#' Preprocess conditions
#' 
#' This function prepares conditions for computing conditional circuit paramaters via cforde
#' 
#' @param params Circuit parameters learned via \code{\link{forde}}. 
#' @param condition Optional set of conditioning events.
#' @param row_mode Interpretation of rows in multi-row conditions.
#' 
#' @import data.table
#' @import stringr
#' 

prep_cond <- function(condition, params, row_mode) {
  
  # To avoid data.table check issues
  c_idx <- family <- val <- variable <- NULL
  
  n_row_cond <- nrow(condition)
  meta <- params$meta
  cat <- params$cat
  cnt_cols <- intersect(meta[family != "multinom", variable], colnames(condition))
  cat_cols <- intersect(meta[family == "multinom", variable], colnames(condition))
  
  cond <- copy(condition)
  cond <- setDT(cond)
  if (length(cat_cols > 0)) {
    cond[,(cat_cols) := lapply(.SD,as.character),.SDcols = cat_cols]  
  }
  
  cols_check_range <- cond[,sapply(.SD, \(x) sum((str_sub(x,,1) == "(") | is.na(x))), .SDcols = cnt_cols]
  cols_check_or <- cond[,sapply(.SD, \(x) sum(str_detect(x, "\\|")))]
  if (row_mode == "or") {
    if (any(cols_check_range > 0 & cols_check_range < nrow(cond))){
      stop("Condition vector contains columns with both range and scalar entries. No valid conditional density can be calculated.")
    }
    if (length(cnt_cols) > 0 & n_row_cond > 1) {
      cond[,(cnt_cols) := lapply(.SD,\(col) replace(col, which(is.na(col)), "(-Inf,Inf)")),.SDcols = cnt_cols]
    }
    if (length(cat_cols) > 0 & n_row_cond > 1) {
      cond[,(cat_cols) := mapply(\(colname,col) {
        lvls_str <- paste(levels(as.factor(cat[cat$variable == colname]$val)),collapse="|")
        list(replace(col,which(is.na(col)),lvls_str))
      },colname = colnames(.SD),
      col = .SD),
      .SDcols = cat_cols]
    }
    cond <- apply(cond,1,str_split,"\\|")
    cond <- rbindlist(lapply(cond,expand.grid))
    cond <- unique(cond)
    names(cond) <- meta[, variable]
  } else if (row_mode == "separate") {
    if (any(cols_check_or > 0,na.rm = T)) {
      stop("Please use the option row_mode = 'or' when including logical disjunctions.")
    }
  }
  
  cond[,c_idx := .I]
  
  suppressWarnings(
    condition_long <- melt(cond,id.vars = "c_idx", value.name = "val")[!is.na(val),]
  )
  
  if(row_mode == "or") {
    condition_long[,c_idx:= .GRP, by = c_idx]
  }
  condition_long[(variable %in% cnt_cols) & str_detect(val,"\\("),c("val", "min", "max") := cbind(c(NA_real_,transpose(strsplit(substr(val, 2, nchar(val) - 1), split = ","))))]
  condition_long[, c("min", "max") := lapply(.SD,as.numeric), .SDcols = c("min","max")]
  setcolorder(condition_long, c("c_idx", "variable", "min", "max"))
  
  if (row_mode == 'or' & nrow(cond) > 1 & !all(cols_check_range == 0)) {
    condition_long <- unoverlap_hyperrectangles(condition_long, cat_cols)
    condition_long <- condition_long[!(min == -Inf & max == Inf & is.na(val))]
  }
  setorder(condition_long, c_idx)
  condition_long[]
}



#' Unoverlap hyperrectangles for or-ed conditions
#' 
#' Algorithm is not optimal ( = number of resulting hyperrectangles is not minimal).
#' 
#' @param hyperrectangles Input data.table of hyperrectangles in long format
#' @param cat_cols Vector of names of categorical features
#' 
#' @import data.table
#'

unoverlap_hyperrectangles <- function(hyperrectangles, cat_cols) {
  
  # To avoid data.table check issues
  . <- val <- variable <- val_fac <- bound <- c_idx <- n_partials_var <- N <-
    n_subhyperrectangles <- n_subhyperrectangles_unoverlapped <- subvolume_id <-
    NULL
  
  scalar_cols <- hyperrectangles[!is.na(val) & !(variable %in% cat_cols),unique(as.character(variable))]
  hyperrectangles[variable %in% c(scalar_cols, cat_cols), val_fac := as.numeric(as.factor(val))]
  hyperrectangles[variable %in% c(scalar_cols, cat_cols), c("min", "max") := .(val_fac - 0.1, val_fac + 0.1)]
  
  d <- uniqueN(hyperrectangles[,variable]) 
  
  partitioned_space <- setorder(unique(data.table(variable = hyperrectangles[,variable], bound = hyperrectangles[,c(min,max)])))
  
  hyperrectangles_split <- hyperrectangles[,{
    var <- variable
    indices <- partitioned_space[(variable == var) & (bound %in% c(min,max)),,which=T]
    indices_min <- indices[1]
    indices_max <- indices[2]
    list(min = partitioned_space[indices_min:(indices_max-1),bound],max = partitioned_space[(indices_min+1):indices_max,bound],n_partials_var = indices_max - indices_min)
  },by = .(c_idx,variable)]
  
  hyperrectangles_split_cartesian <- hyperrectangles_split[,{
    lengths <- unique(.SD[,.(variable,n_partials_var)])[,n_partials_var]
    cartesian <- expand.grid(lapply(lengths,seq_len))
    subsetted <- .SD[as.vector(unlist(transpose(cartesian))) + rep(c(0,cumsum(lengths)[-d]),nrow(cartesian))]
    subsetted[, `:=` (n_partials_var = NULL, subvolume_id = rep(seq_len(.N/d),each=d), n_subhyperrectangles = .N/d)]
  }, by = c_idx]
  
  
  hyperrectangles_split_cartesian_wide <- dcast(hyperrectangles_split_cartesian, c_idx + subvolume_id + n_subhyperrectangles ~ variable, value.var =  c("min","max"))
  hyperrectangles_split_cartesian_wide[,N:=.N,by = eval(names(hyperrectangles_split_cartesian_wide)[4:(3+2*d)])][,N:=sum(N > 1), by = c_idx][,N:=N/n_subhyperrectangles]
  setorder(hyperrectangles_split_cartesian_wide, N, n_subhyperrectangles)
  
  hyperrectangles_unoverlapped_fragmented_wide <- unique(hyperrectangles_split_cartesian_wide,by = 4:(3+2*d))[,n_subhyperrectangles_unoverlapped := .N, by = c_idx]
  
  hyperrectangles_unchanged <- hyperrectangles[c_idx %in% hyperrectangles_unoverlapped_fragmented_wide[n_subhyperrectangles == n_subhyperrectangles_unoverlapped,c_idx]]
  hyperrectangles_changed_wide <- hyperrectangles_unoverlapped_fragmented_wide[n_subhyperrectangles != n_subhyperrectangles_unoverlapped]
  hyperrectangles_changed <-  hyperrectangles_split_cartesian[c_idx %in% hyperrectangles_changed_wide[,c_idx] & subvolume_id %in% hyperrectangles_changed_wide[,subvolume_id],]
  hyperrectangles_changed[,c_idx := hyperrectangles_unchanged[,max(c_idx)] + .GRP, by = .(c_idx, subvolume_id)]
  hyperrectangles_unoverlapped <- rbind(hyperrectangles_unchanged[,1:4],hyperrectangles_changed[,1:4])
  hyperrectangles_unoverlapped[variable %in% c(cat_cols, scalar_cols), c("min", "max", "val_fac") := .(NA_real_, NA_real_, round(min))]
  hyperrectangles_unoverlapped[variable %in% c(scalar_cols,cat_cols), val := merge(hyperrectangles_unoverlapped, unique(hyperrectangles[variable %in% c(scalar_cols, cat_cols),.(variable, val_fac, val)]), by = c("variable", "val_fac"))[,val]]
  setcolorder(hyperrectangles_unoverlapped, names(hyperrectangles))  
  hyperrectangles_unoverlapped[,-c("val_fac")]
}

