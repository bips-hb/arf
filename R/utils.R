#' Adaptive column renaming
#' 
#' This function renames columns in case the input colnames includes any
#' colnames required by internal functions (e.g., \code{"y"}).
#' 
#' @param cn Column names.
#' @param old_name Name of column to be renamed.
#' @keywords internal

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
#' @param cn Old column names.
#'
#' @return New columns names.
#' @keywords internal

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
#' @keywords internal

resample <- function(x, ...) {
  x[sample.int(length(x), ...)]
}

#' which.max() with random at ties
#'
#' @param x A numeric vector.
#'
#' @return Index of maximum value in x, with random tie-breaking.
#' @keywords internal

which.max.random <- function(x) {
  if (all(is.na(x))) {
    return(NA)
  }
  which(rank(x, ties.method = "random", na.last = FALSE) == length(x))
}

#' Preprocess input data
#' 
#' This function prepares input data for ARFs.
#' 
#' @param x Input data.frame.
#' @param verbose Show warning if recoding integers?
#' @keywords internal

prep_x <- function(x, verbose = TRUE) {
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
      if (verbose) {
        warning('Recoding integers with more than 5 unique values as numeric. ', 
              'To override this behavior, explicitly code these variables as factors.')
      }
      x[, to_numeric] <- lapply(x[, to_numeric, drop = FALSE], as.numeric)
    }
    to_factor <- sapply(seq_len(ncol(x)), function(j) {
      idx_integer[j] & length(unique(x[[j]])) < 6
    })
    if (any(to_factor)) {
      if (verbose) {
        warning('Recoding integers with fewer than 6 unique values as ordered factors. ', 
              'To override this behavior, explicitly code these variables as numeric.')
      }
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
#' @param round Round continuous variables to their respective maximum precision in the real data set?
#' 
#' @import data.table
#' @keywords internal

post_x <- function(x, params, round = TRUE) {
  
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
  if (sum(idx_numeric) > 0L & round) {
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
         if (round) {
           as.integer(round(x[[j]]))
         } else {
           x[[j]]
         }
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
#' @param nomatch What to do if no leaf matches a condition in \code{evidence}?
#'   Options are to force sampling from a random leaf (\code{"force"}) or return 
#'   \code{NA} (\code{"na"}). The default is \code{"force"}.
#' @param verbose Show warnings, e.g. when no leaf matches a condition?   
#' @param stepsize Stepsize defining number of condition rows handled in one for each step.
#' @param parallel Compute in parallel? Requires a registered \code{foreach}
#'   backend (\code{doParallel}, \code{doFuture}) or active \code{mirai}
#'   daemons. See \code{\link{arf-options}}.
#'   
#' @return List with conditions (\code{evidence_input}), prepared conditions (\code{evidence_prepped})
#'   and leaves that match the conditions in evidence with continuous data (\code{cnt}) 
#'   and categorical data (\code{cat}) as well as leaf info (\code{forest}).
#' 
#' @import data.table
#' @importFrom foreach foreach %dopar%
#' @importFrom truncnorm dtruncnorm ptruncnorm 
#' @importFrom stats dunif punif
#' @keywords internal

cforde <- function(params, 
                   evidence, 
                   row_mode = c("separate", "or"), 
                   nomatch = c("force", "na"),
                   verbose = TRUE,
                   stepsize = 0, 
                   parallel = TRUE) {
  
  row_mode <- match.arg(row_mode)
  nomatch <- match.arg(nomatch)
  
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
      stepsize <- ceiling(nconds_conditioned/arf_n_workers())
    } else {
      stepsize <- nconds_conditioned
    }
  }
  step_no <- ceiling(nconds_conditioned/stepsize)
  
  # Per-step work lives once in arf_cforde_step() (cforde_workers.R); this closure
  # adapts it to foreach's one-argument iteration.
  par_fun <- function(step_) {
    arf_cforde_step(step_, condition_long, conds_conditioned, nconds_conditioned,
                    stepsize, cat_cols, cnt_cols, params, family)
  }
  # Parallelism is across condition steps: 1 step is inherently serial. mirai only
  # for step_no > 1. Under "separate" evidence_row_mode the caller (forge/expct)
  # already parallelizes and passes parallel = FALSE, so no daemon-in-daemon.
  use_mirai <- FALSE
  if (step_no > 1) {
    backend <- arf_select_backend(parallel)
    use_mirai <- identical(backend, "mirai")
  }
  if (use_mirai) {
    arf_load_on_daemons()  # daemons call arf/data.table internals via bare names
    params_shared <- mori::share(params)
    condition_long_shared <- mori::share(condition_long)
    # arf_cforde_step returns list(cnt_new, cat_new, relevant_leaves), so combine
    # each component across steps (arf_mirai_tree_map assumes data.table results).
    res <- mirai::mirai_map(
      seq_len(step_no),
      arf_cforde_step,
      .args = list(condition_long = condition_long_shared,
                   conds_conditioned = conds_conditioned,
                   nconds_conditioned = nconds_conditioned, stepsize = stepsize,
                   cat_cols = cat_cols, cnt_cols = cnt_cols,
                   params = params_shared, family = family))[]
    arf_stop_on_mirai_error(res)
    updates_relevant_leaves <- list(
      cnt_new = rbindlist(lapply(res, `[[`, "cnt_new")),
      cat_new = rbindlist(lapply(res, `[[`, "cat_new")),
      relevant_leaves = rbindlist(lapply(res, `[[`, "relevant_leaves")))
  } else if (isTRUE(parallel) && step_no > 1) {
    updates_relevant_leaves <- foreach(step = 1:step_no, .combine = "rbind") %dopar% par_fun(step)
  } else {
    updates_relevant_leaves <- foreach(step = 1:step_no, .combine = "rbind") %do% par_fun(step)
  }
  
  # Combine results
  if (is.matrix(updates_relevant_leaves)) {
    updates_relevant_leaves <- lapply(as.data.table(updates_relevant_leaves), rbindlist) 
  }
  
  # Re-index matching leaves
  relevant_leaves <- updates_relevant_leaves$relevant_leaves[,`:=` (f_idx = .I, f_idx_uncond = f_idx)][]
  cnt_new <- setcolorder(merge(relevant_leaves, updates_relevant_leaves$cnt_new, by.x = c("c_idx", "f_idx_uncond"), by.y = c("c_idx", "f_idx"), sort = FALSE), c("f_idx","c_idx","variable","min","max","val","cvg_factor"))[]
  cat_new <- setcolorder(merge(relevant_leaves, updates_relevant_leaves$cat_new, by.x = c("c_idx", "f_idx_uncond"), by.y = c("c_idx", "f_idx"), sort = FALSE), c("f_idx","c_idx","variable","val","prob","cvg_factor"))[]
  
  # Check for conditions with no matching leaves and handle this according to row_mode
  conds_impossible <- conds_conditioned[!(conds_conditioned %in% relevant_leaves[,unique(c_idx)])]
  if (relevant_leaves[,uniqueN(c_idx)] < nconds_conditioned) {
    if (relevant_leaves[,uniqueN(c_idx)] == 0 & row_mode == "or") {
      stop("For all entered evidence rows, no matching leaves could be found. This is probably because evidence lies outside of the distribution calculated by FORDE. For continuous data, consider setting epsilon>0 or finite_bounds='no' in forde(). For categorical data, consider setting alpha>0 in forde().")
    } else {
      if (verbose) {
        wrn <- "For some entered evidence rows, no matching leaves could be found. This is probably because evidence lies outside of the distribution calculated by FORDE. For continuous data, consider setting epsilon>0 or finite_bounds='no' in forde(). For categorical data, consider setting alpha>0 in forde()."
        if (nomatch == "force") {
          warning(paste(wrn, "Sampling from all leaves with equal probability (can be changed with 'nomatch' argument)."))
        } else {
          warning(paste(wrn, "Returning NA for those rows (can be changed with 'nomatch' argument)."))
        }
      }
      impossible_leaves <- data.table(c_idx = conds_impossible, f_idx = NA_integer_, f_idx_uncond = NA_integer_)
      relevant_leaves <- setorder(rbind(relevant_leaves, impossible_leaves))
    }
  }
  
  # Calculate new forest (set of leaves and weights)
  forest_new <- merge(relevant_leaves, forest, by.x = "f_idx_uncond", by.y = "f_idx", all.x = TRUE, sort = FALSE)
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
        if (nomatch == "force") {
          cvg_new[, cvg := 1/.N]
        } else {
          cvg_new[, cvg := NA]
        }
        if (verbose) {
          wrn <- "All leaves have zero likelihood. This is probably because evidence contains an (almost) impossible combination."
          if (nomatch == "force") {
            warning(paste(wrn, "Sampling from all possible leaves with equal probability."))
          } else {
            warning(paste(wrn, "Returning NA."))
          }
        }
      } else {
        cvg_new[, cvg := exp(cvg - max(cvg))]
        cvg_new <- cvg_new[, cvg := cvg / sum(cvg)]
      }
    } else {
      cvg_new[, leaf_zero_lik := all(cvg == -Inf), by = c_idx]
      if (any(cvg_new[, leaf_zero_lik])) {
        if (nomatch == "force") {
          cvg_new[leaf_zero_lik == TRUE, cvg := 1/.N, by = c_idx]
        } else {
          cvg_new <- cvg_new[leaf_zero_lik == FALSE, ]
        }
        if (verbose) {
          wrn <- "All leaves have zero likelihood for some entered evidence rows. This is probably because evidence contains an (almost) impossible combination."
          if (nomatch == "force") {
            warning(paste(wrn, "Sampling from all possible leaves with equal probability (can be changed with 'nomatch' argument)."))
          } else {
            warning(paste(wrn, "Returning NA for those rows (can be changed with 'nomatch' argument)."))
          }
        }
      }
      if (any(cvg_new[, !leaf_zero_lik])) {
        cvg_new[leaf_zero_lik == FALSE, scale := max(cvg), by = c_idx]
        cvg_new[leaf_zero_lik == FALSE, cvg := exp(cvg - scale)]
        cvg_new[leaf_zero_lik == FALSE, scale := sum(cvg), by = c_idx]
        cvg_new[leaf_zero_lik == FALSE, cvg := cvg / scale]
        cvg_new[, scale := NULL]
      }
      cvg_new[, leaf_zero_lik := NULL]
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
  if (row_mode == "separate" & nconds != nconds_conditioned) {
    conds_unconditioned <- (1:nconds)[!(1:nconds) %in% conds_conditioned]
    forest_new_unconditioned <- copy(forest)
    forest_new_unconditioned <- rbindlist(replicate(length(conds_unconditioned), forest, simplify = FALSE))
    forest_new_unconditioned[, `:=` (c_idx = rep(conds_unconditioned,each = nrow(forest)), f_idx_uncond = f_idx, cvg_arf = cvg)]
    forest_new <- rbind(forest_new, forest_new_unconditioned)
  }
  
  setorder(setcolorder(forest_new,c("f_idx","c_idx","f_idx_uncond","tree","leaf","cvg_arf","cvg")), c_idx, f_idx, f_idx_uncond, tree, leaf)
  
  list(evidence_input = evidence, evidence_prepped = condition_long, cnt = cnt_new, cat = cat_new, forest = forest_new)
}


#' Preprocess conditions
#' 

#' 
#' @param params Circuit parameters learned via \code{\link{forde}}. 
#' @param evidence Optional set of conditioning events.
#' @param row_mode Interpretation of rows in multi-row conditions.
#' 
#' @import data.table
#' @import stringr
#' @keywords internal

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
  
  cols_check_range <- cond[,sapply(.SD, function(x) sum((str_sub(x,,1) == "(") | is.na(x))), .SDcols = cnt_cols]
  cols_check_or <- cond[,sapply(.SD, function(x) sum(str_detect(x, "\\|")))]
  
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
