#' Adaptive column renaming
#' 
#' This function renames columns in case the input colnames includes any
#' colnames required by internal functions (e.g., \code{"y"}).
#' 
#' @param cn Column names.
#' @param old_name Name of column to be renamed.
#' 

col_rename <- function(cn, old_name) {
  k <- 1L
  converged <- FALSE
  while (!isTRUE(converged)) {
    new_name <- paste0(old_name, k)
    if (!new_name %in% cn) {
      converged <- TRUE
    } else {
      k <- k + 1L
    }
  }
  return(new_name)
}

#' Rename all problematic columns with col_rename().
#'
#' @param Old column names.
#'
#' @return New columns names.
#'
#' @examples
col_rename_all <- function(cn) {
  if ('y' %in% cn) {
    cn[which(cn == 'y')] <- col_rename(cn, 'y')
  }
  if ('obs' %in% cn) {
    cn[which(cn == 'obs')] <- col_rename(cn, 'obs')
  }
  if ('tree' %in% cn) {
    cn[which(cn == 'tree')] <- col_rename(cn, 'tree')
  }
  if ('leaf' %in% cn) {
    cn[which(cn == 'leaf')] <- col_rename(cn, 'leaf')
  }
  if ('cnt' %in% cn) {
    cn[which(cn == 'cnt')] <- col_rename(cn, 'cnt')
  }
  if ('N' %in% cn) {
    cn[which(cn == 'N')] <- col_rename(cn, 'N')
  }
  cn
}

#' Safer version of sample()
#'
#' @param x A vector of one or more elements from which to choose.
#' @param ... Further arguments for sample().
#'
#' @return A vector of length size with elements drawn from x.

resample <- function(x, ...) {
  x[sample.int(length(x), ...)]
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
  colnames(x) <- col_rename_all(colnames(x))
  return(x)
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
      if (is.numeric(x[[j]])) {
        as.integer(round(x[[j]]))
      } else {
        as.integer(as.character(x[[j]]))
      }
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
#' @param evidence Data frame of conditioning event(s).
#' @param row_mode Interpretation of rows in multi-row conditions.
#' @param stepsize Stepsize defining number of condition rows handled in one for each step.
#' @param parallel Compute in parallel? Must register backend beforehand, e.g. 
#'   via \code{doParallel}.
#'   
#' @return List with conditions (\code{evidence_input}), prepared conditions (\code{evidence_prepped})
#'   and leaves that match the conditions in evidence with continuous data (\code{cnt}) 
#'   and categorical data (\code{cat}) as well as leaf info (\code{forest}).
#' 
#' @import data.table
#' @importFrom foreach foreach %dopar%
#' @importFrom truncnorm dtruncnorm ptruncnorm 
#' @importFrom stats dunif punif
#' 

cforde <- function(params, evidence, row_mode = c("separate", "or"), stepsize = 0, parallel = TRUE) {
  
  row_mode <- match.arg(row_mode)
  
  # To avoid data.table check issues
  . <- c_idx <- cvg <- cvg_arf <- cvg_factor <- f_idx <- f_idx_uncond <- i.max <-
    i.min <- leaf <- max.x <- max.y <- min.x <- min.y <- mu <- prob <- sigma <-
    step_ <-	tree <-	V1 <- val <- variable <- leaf_zero_lik <- step <- NULL
  
  # Store informations of params as variables
  meta <- params$meta
  family <- meta[family != "multinom", unique(family)]
  forest <- params$forest
  cat <- params$cat
  cnt <- params$cnt
  cnt_cols <- meta[family != "multinom", variable]
  cat_cols <- meta[family == "multinom", variable]
  
  # Calculate long format of evidence depending on row_mode
  condition_long <- prep_cond(evidence, params, row_mode)
  setkey(condition_long, c_idx)
  
  # If evidence does not any conditions (i.e. all entries equal NA), return NULL
  if (nrow(condition_long) == 0){
    return(NULL)
  }
  
  # Store number of evidence rows and number of evidence rows that do not consist of NA only
  if (row_mode == "or") {
    nconds <- nconds_conditioned <- condition_long[, max(c_idx)]
  } else {
    nconds <- condition_long[, max(c_idx)]
    nconds_conditioned <- condition_long[,uniqueN(c_idx)]
  }
  
  # Store set of condition (from evidence rows that do not consist of NA only)
  conds_conditioned <- condition_long[, unique(c_idx)]
  
  # Calculate stepsize for parallelization depending on number of conditions and registered workers
  if (stepsize == 0) {
    if (parallel) {
      stepsize <- ceiling(nconds_conditioned/getDoParWorkers())
    } else {
      stepsize <- nconds_conditioned
    }
  }
  step_no <- ceiling(nconds_conditioned/stepsize)
  
  # Loop through conditions with defined stepsize to determined matching leaves and updates for cat and cnt params
  update_fun <-function(step_) {
    
    # Define subset of conditions for step_
    index_start <- conds_conditioned[(step_ - 1)*stepsize + 1]
    index_end <- conds_conditioned[min(step_ * stepsize, nconds_conditioned)]
    condition_long_step <- condition_long[.(index_start:index_end),nomatch = NULL]
    
    # Store cat and cnt conditions separately
    cat_conds <- condition_long_step[variable %in% cat_cols,c("c_idx","variable","val")][, variable := factor(variable)]
    cnt_conds <- condition_long_step[variable %in% cnt_cols,c("c_idx","variable","min", "max","val")][,`:=` (variable = factor(variable),
                                                                                                             val = as.numeric(val))]
    # If cat conditions exist, calculate matching leaves
    if (nrow(cat_conds) != 0) {
      
      # Save leaf indices f_idx cat params in list column grouped by variable and val (value) and merge with cat conditions
      cat_relevant <- cat[, .(.(f_idx)), by=.(variable,val)]
      setkey(cat_relevant, variable, val)
      setkey(cat_conds, variable, val)
      cat_relevant <- cat_conds[cat_relevant, on = .(variable, val), nomatch = NULL]
      setkey(cat_relevant, c_idx)
      
      # or-combine different conditions on the same feature
      if (uniqueN(cat_relevant, by = c("c_idx", "variable")) != nrow(cat_relevant)) {
        cat_relevant <- cat_relevant[, .(.(Reduce(union, V1))), by = .(c_idx, variable)]
      }
      
      # Determine matching leaves for cat conditions
      relevant_leaves_changed_cat <- cat_relevant[, Reduce(intersect, V1), by = c_idx][, .(c_idx, f_idx = V1)]
      setorder(relevant_leaves_changed_cat)
      conditions_unchanged_cat <- setdiff(condition_long_step[, c_idx], cat_conds[, c_idx])
      relevant_leaves_unchanged_cat <- data.table(c_idx = rep(conditions_unchanged_cat, each = nrow(forest) ), f_idx = rep(forest[,f_idx],length(conditions_unchanged_cat)))
      relevant_leaves_cat <- rbind(relevant_leaves_changed_cat, relevant_leaves_unchanged_cat, fill = T)
      relevant_leaves_cat_list <- relevant_leaves_cat[,.(f_idx = .(f_idx)), by=c_idx]
    } else {
      relevant_leaves_cat <- data.table(c_idx = integer(), f_idx = integer())
    }
    
    # If cnt conditions exist, calculate matching leaves
    if (nrow(cnt_conds) != 0) {
      
      # Save min, max in cnt params in list columns grouped by variable and merge with cnt conditions
      cnt_relevant <- cnt[, .(min = .(min), max = .(max)), by = variable]
      cnt_conds_compact <- copy(cnt_conds)
      cnt_conds_compact[!is.na(val), `:=`(min = val, max = val)][, val := NULL]
      cnt_relevant <- cnt_conds_compact[cnt_relevant, on = .(variable), nomatch = NULL]
      setkey(cnt_relevant, c_idx)
      
      # If cat conds exist, use only matching subset of potentially relevant leaves for cnt conditions
      if (nrow(cat_conds) != 0) {
        cnt_relevant <- cnt_relevant[relevant_leaves_cat_list, on = .(c_idx)]
      } else {
        cnt_relevant[, f_idx := NA]
      }
      
      # Determine matching leaves for cnt conditions per row
      cnt_relevant <- cnt_relevant[, .(
        c_idx,
        variable,
        f_idx = Map(\(f_idx, min, max, i.min, i.max) {
          if (!inherits(f_idx, "logical")) {
            rel_cnt_min <- i.min[f_idx]
            rel_cnt_max <- i.max[f_idx]
            rel_min <- f_idx[which(max > rel_cnt_min)]
            rel_max <- f_idx[which(min <= rel_cnt_max)]
          } else {
            rel_cnt_min <- i.min
            rel_cnt_max <- i.max
            rel_min <- which(max > rel_cnt_min)
            rel_max <- which(min <= rel_cnt_max)
          }
          intersect(rel_min,rel_max) 
        }, f_idx = f_idx, min = min, max = max, i.min = i.min, i.max = i.max))]
      
      # or-combine different conditions on the same feature
      or_within_row_cnt <- uniqueN(cnt_relevant, by = c("c_idx", "variable")) != nrow(cnt_relevant)
      if (or_within_row_cnt) {
        cnt_relevant <- cnt_relevant[, .(f_idx = .(Reduce(union, f_idx))), by = .(c_idx, variable)]
      }
      
      # Determine matching leaves for cnt conditions
      relevant_leaves_changed_cnt <- cnt_relevant[, Reduce(intersect, f_idx),by = c_idx][,.(c_idx, f_idx = V1)]
      conditions_unchanged_cnt <- setdiff(condition_long_step[, c_idx], cnt_conds[, c_idx])
      relevant_leaves_unchanged_cnt <- data.table(c_idx = rep(conditions_unchanged_cnt, each = nrow(forest)), f_idx = rep(forest[,f_idx],length(conditions_unchanged_cnt)))
      relevant_leaves_cnt <- rbind(relevant_leaves_changed_cnt, relevant_leaves_unchanged_cnt)
      
      # Calculate updates for cnt params matching cnt conditions
      cnt_new <- merge(merge(relevant_leaves_cnt, cnt_conds, by = "c_idx", allow.cartesian = T, sort = F), cnt, by = c("f_idx", "variable"), all.x = T, allow.cartesian = T, sort = F)
      cnt_new[!is.na(val),`:=` (min = min.y,
                                max = max.y)]
      cnt_new[is.na(val),`:=` (min = pmax(min.x, min.y, na.rm = T),
                               max = pmin(max.x, max.y, na.rm = T))]
      cnt_new <- cnt_new[min <= max & sum(max.x, val, na.rm = T) != min.y, ]
      cnt_new[, prob := NA_real_]
      if (family == "truncnorm") {
        cnt_new[!is.na(val), prob := dtruncnorm(val, a=min.y, b=max.y, mean=mu, sd=sigma)*(val != min.y)]
        cnt_new[is.na(val) & (min == min.y) & (max == max.y), prob := 1]
        cnt_new[is.na(val) & is.na(prob), prob := ptruncnorm(max, a=min.y, b=max.y, mean=mu, sd=sigma) - ptruncnorm(min, a=min.y, max.y, mean=mu,sd=sigma)] 
      } else if (family == "unif") {
        cnt_new[!is.na(val), prob := dunif(val, min=min.y, max=max.y)*(val != min.y)]
        cnt_new[is.na(val) & (min == min.y) & (max == max.y), prob := 1]
        cnt_new[is.na(val) & is.na(prob), prob := punif(max, min=min.y, max=max.y) - punif(min, min=min.y, max.y)]  
      }
      cnt_new[, c("min.x","max.x","min.y","max.y") := NULL]
      
      # If or-combined cnt condition within rows exist, calculate likelihoods for ranges within leaves and norm to 1
      if (or_within_row_cnt) {
        cnt_new[, cvg_factor := sum(prob), by = .(f_idx, c_idx, variable)]
        cnt_new[, prob := prob/cvg_factor]
      } else {
        cnt_new[, `:=` (cvg_factor = prob, prob = 1)]
      }
      
      # Calculate final set of matching leaves
      if (nrow(cat_conds) > 0) {
        relevant_leaves <- merge(relevant_leaves_cnt, relevant_leaves_cat, by = c("c_idx", "f_idx"))[,.(c_idx, f_idx)]
      } else {
        relevant_leaves <- relevant_leaves_cnt[,.(c_idx, f_idx)]
      }
      
      # If no cnt conditions exist, output empty update table cnt_new for cnt params
    } else {
      relevant_leaves <- relevant_leaves_cat[,.(c_idx, f_idx)]
      cnt_new <- cbind(cnt[F,], data.table(cvg_factor = numeric(), c_idx = integer(), val = numeric(), prob = numeric()))
    }
    
    # Calculate updates for cat params matching cat conditions
    cat_new <- merge(merge(relevant_leaves, cat_conds, by = "c_idx", allow.cartesian = T), cat, by = c("f_idx","variable", "val")) 

    # Ensure probabilities sum to 1
    cat_new[, cvg_factor := sum(prob), by = .(f_idx, c_idx, variable)]
    cat_new[, prob := prob/cvg_factor]
    
    list(cnt_new = cnt_new, cat_new = cat_new, relevant_leaves = relevant_leaves)
  }
  if (isTRUE(parallel)) {
    updates_relevant_leaves <- foreach(step = 1:step_no, .combine = "rbind") %dopar% update_fun(step)
  } else {
    updates_relevant_leaves <- foreach(step = 1:step_no, .combine = "rbind") %do% update_fun(step)
  }
  
  # Combine results
  if (is.matrix(updates_relevant_leaves)) {
    updates_relevant_leaves <- lapply(as.data.table(updates_relevant_leaves), rbindlist) 
  }
  
  # Re-index matching leaves
  relevant_leaves <- updates_relevant_leaves$relevant_leaves[,`:=` (f_idx = .I, f_idx_uncond = f_idx)][]
  cnt_new <- setcolorder(merge(relevant_leaves, updates_relevant_leaves$cnt_new, by.x = c("c_idx", "f_idx_uncond"), by.y = c("c_idx", "f_idx"), sort = F), c("f_idx","c_idx","variable","min","max","val","cvg_factor"))[]
  cat_new <- setcolorder(merge(relevant_leaves, updates_relevant_leaves$cat_new, by.x = c("c_idx", "f_idx_uncond"), by.y = c("c_idx", "f_idx"), sort = F), c("f_idx","c_idx","variable","val","prob","cvg_factor"))[]

  # Check for conditions with no matching leaves and handle this according to row_mode
  if (relevant_leaves[,uniqueN(c_idx)] < nconds_conditioned) {
    if (relevant_leaves[,uniqueN(c_idx)] == 0 & row_mode == "or") {
      stop("For all entered evidence rows, no matching leaves could be found. This is probably because evidence lies outside of the distribution calculated by FORDE. For continuous data, consider setting epsilon>0 or finite_bounds='no' in forde(). For categorical data, consider setting alpha>0 in forde()")
    } else {
      warning("For some entered evidence rows, no matching leaves could be found. This is probably because evidence lies outside of the distribution calculated by FORDE. For continuous data, consider setting epsilon>0 or finite_bounds='no' in forde(). For categorical data, consider setting alpha>0 in forde()")
      conds_impossible <- conds_conditioned[!(conds_conditioned %in% relevant_leaves[,unique(c_idx)])]
      relevant_leaves <- setorder(rbind(relevant_leaves, data.table(c_idx = conds_impossible, f_idx = NA_integer_, f_idx_uncond = NA_integer_)))
    }
  }
  
  # Calculate new forest (set of leaves and weights)
  forest_new <- merge(relevant_leaves, forest, by.x = "f_idx_uncond", by.y = "f_idx", all.x = T, sort = F)
  setnames(forest_new, "cvg", "cvg_arf")
  
  cvg_new <- unique(rbind(cat_new[, .(f_idx, c_idx, variable, cvg_factor)],
               cnt_new[, .(f_idx, c_idx, variable, cvg_factor)]),
         by = c("f_idx", "variable"))[,-"variable"]
  
  if (nrow(cvg_new) > 0) {
    # Use log transformation to avoid overflow
    cvg_new[, cvg_factor := log(cvg_factor)]
    cvg_new <- cvg_new[, .(cvg_factor = sum(cvg_factor)), keyby = f_idx]
    cvg_new <- cbind(cvg_new, forest_new[!is.na(cvg_arf), .(c_idx, cvg_arf = log(cvg_arf))])
    cvg_new[,`:=` (cvg = cvg_factor + cvg_arf, cvg_factor = NULL, cvg_arf = NULL)]
    
    # Re-calculate weights and transform back from log scale, handle (numerically) impossible cases
    if (row_mode == "or") {
      if (cvg_new[,all(cvg == -Inf)]) {
        warning("All leaves have zero likelihood. This is probably because evidence contains an (almost) impossible combination.")
        cvg_new[, cvg := 1/.N]
      } else {
        cvg_new[, cvg := exp(cvg - max(cvg))]
        cvg_new <- cvg_new[, cvg := cvg / sum(cvg)]
      }
    } else {
      cvg_new[, leaf_zero_lik := all(cvg == -Inf), by = c_idx]
      if (any(cvg_new[, leaf_zero_lik])) {
        warning("All leaves have zero likelihood for some entered evidence rows. This is probably because evidence contains an (almost) impossible combination.")
        cvg_new[leaf_zero_lik == TRUE, cvg := 1/.N, by = c_idx]
      }
      cvg_new[leaf_zero_lik == FALSE, scale := max(cvg), by = c_idx]
      cvg_new[leaf_zero_lik == FALSE, cvg := exp(cvg - scale)]
      cvg_new[leaf_zero_lik == FALSE, scale := sum(cvg), by = c_idx]
      cvg_new[leaf_zero_lik == FALSE, cvg := cvg / scale]
      cvg_new[, `:=` (leaf_zero_lik = NULL, scale = NULL)]
    }
  }
  
  # Add conditions with no matching leaves to forest
  forest_new_noleaf <- data.table(c_idx = setdiff(unique(forest_new[,c_idx]), unique(cvg_new[,c_idx])))[,f_idx := NA_integer_]
  forest_new <- merge(forest_new, cvg_new[,.(f_idx, cvg)], by = "f_idx")
  forest_new <- forest_new[cvg > 0,]
  forest_new <- rbind(forest_new, forest_new_noleaf, fill = TRUE)
  if (row_mode == "or") {
    if (forest_new[,all(is.na(f_idx))]) {
      forest_new[is.na(f_idx), cvg := 1/.N]
    } else {
      forest_new[is.na(f_idx), cvg := 0]
    }
  } else {
    forest_new[is.na(f_idx), cvg := 1]
  }
  
  # Add all leaves for all-NA conditions to forest
  if (row_mode == "separate" & (nconds != nconds_conditioned)) {
    conds_unconditioned <- (1:nconds)[!(1:nconds) %in% conds_conditioned]
    forest_new_unconditioned <- copy(forest)
    forest_new_unconditioned <- rbindlist(replicate(length(conds_unconditioned), forest, simplify = F))
    forest_new_unconditioned[, `:=` (c_idx = rep(conds_unconditioned,each = nrow(forest)), f_idx_uncond = f_idx, cvg_arf = cvg)]
    forest_new <- rbind(forest_new, forest_new_unconditioned)
  }

  setorder(setcolorder(forest_new,c("f_idx","c_idx","f_idx_uncond","tree","leaf","cvg_arf","cvg")), c_idx, f_idx, f_idx_uncond, tree, leaf)
  
  list(evidence_input = evidence, evidence_prepped = condition_long, cnt = cnt_new, cat = cat_new, forest = forest_new)
}


#' Preprocess conditions
#' 
#' This function prepares conditions for computing conditional circuit paramaters via cforde
#' 
#' @param params Circuit parameters learned via \code{\link{forde}}. 
#' @param evidence Optional set of conditioning events.
#' @param row_mode Interpretation of rows in multi-row conditions.
#' 
#' @import data.table
#' @import stringr
#' 

prep_cond <- function(evidence, params, row_mode) {
  
  # To avoid data.table check issues
  c_idx <- family <- val <- variable <- val._x <- . <- NULL
  
  # If condition already in correct long format, do nothing
  if(ncol(evidence) == 5 && all(names(evidence) == c("c_idx", "variable", "min", "max", "val"))) {
    return(evidence)
  }
  
  n_row_cond <- nrow(evidence)
  meta <- params$meta
  cat <- params$cat
  cnt_cols <- intersect(meta[family != "multinom", variable], colnames(evidence))
  cat_cols <- intersect(meta[family == "multinom", variable], colnames(evidence))
  
  cond <- copy(evidence)
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
  }
  
  cond[, c_idx := .I]
  
  suppressWarnings(
    condition_long <- melt(cond,id.vars = "c_idx", value.name = "val")[!is.na(val),]
  )
  
  # handle logical or within rows
  if (any(cols_check_or > 0, na.rm = TRUE)) {
    condition_long <- condition_long[, .(val = unlist(str_split(val,"\\|"))), by = .(c_idx, variable)]
  }
  
  # Logical not
  cond_lnot <- condition_long[(variable %in% cat_cols) & str_detect(val, "^!"), ]
  if (nrow(cond_lnot) > 0) {
    cond_lnot[, val := str_remove(val, "^!")]
    cond_lnot <- merge(cond_lnot, params$levels, by = "variable", suffixes = c("._x", ""), allow.cartesian = TRUE)
    cond_lnot <- cond_lnot[val != val._x, ][, val._x := NULL]
    condition_long <- rbind(condition_long[!((variable %in% cat_cols) & str_detect(val, "^!")), ], 
                            cond_lnot)
  }
  
  # Interval syntax, e.g. (X,Inf)
  condition_long[(variable %in% cnt_cols) & str_detect(val, "\\("), 
    c("val", "min", "max") := cbind(c(NA_real_, transpose(strsplit(substr(val, 2, nchar(val) - 1), split = ","))))]
  
  # >, < syntax
  condition_long[(variable %in% cnt_cols) & str_detect(val, "<"), 
                 c("val", "min", "max") := list(NA_real_, -Inf, as.numeric(str_remove_all(str_remove_all(val, "\\s"), "<")))]
  condition_long[(variable %in% cnt_cols) & str_detect(val, ">"), 
                 c("val", "min", "max") := list(NA_real_, as.numeric(str_remove_all(str_remove_all(val, "\\s"), ">")), Inf)]
  
  condition_long[, c("min", "max") := lapply(.SD, as.numeric), .SDcols = c("min", "max")]
  setcolorder(condition_long, c("c_idx", "variable", "min", "max"))
  
  setorder(condition_long, c_idx)
  condition_long[]
}
