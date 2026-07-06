#' @keywords internal
#' @noRd

# Per-step worker for cforde(), extracted so foreach and mirai backends share one
# definition (params + condition_long mori-shared under mirai). Calls arf/data.table
# internals via bare names, so mirai daemons must have arf loaded.
arf_cforde_step <- function(step_, condition_long, conds_conditioned,
                            nconds_conditioned, stepsize, cat_cols, cnt_cols,
                            params, family) {
  # data.table NSE silencing
  . <- c_idx <- cvg <- cvg_arf <- cvg_factor <- f_idx <- f_idx_uncond <- i.max <-
    i.min <- leaf <- max.x <- max.y <- min.x <- min.y <- mu <- prob <- sigma <-
    tree <- V1 <- val <- variable <- min <- max <- NULL

  forest <- params$forest
  cat <- params$cat
  cnt <- params$cnt

  # Define subset of conditions for step_
  index_start <- conds_conditioned[(step_ - 1) * stepsize + 1]
  index_end <- conds_conditioned[min(step_ * stepsize, nconds_conditioned)]
  condition_long_step <- condition_long[.(index_start:index_end), nomatch = NULL]

  # Store cat and cnt conditions separately
  cat_conds <- condition_long_step[variable %in% cat_cols, c("c_idx", "variable", "val")][, variable := factor(variable)]
  cnt_conds <- condition_long_step[variable %in% cnt_cols, c("c_idx", "variable", "min", "max", "val")][, `:=`(variable = factor(variable),
                                                                                                              val = as.numeric(val))]
  # If cat conditions exist, calculate matching leaves
  if (nrow(cat_conds) != 0) {

    # Save leaf indices f_idx cat params in list column grouped by variable and val (value) and merge with cat conditions
    cat_relevant <- cat[, .(.(f_idx)), by = .(variable, val)]
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
    relevant_leaves_unchanged_cat <- data.table(c_idx = rep(conditions_unchanged_cat, each = nrow(forest)), f_idx = rep(forest[, f_idx], length(conditions_unchanged_cat)))
    relevant_leaves_cat <- rbind(relevant_leaves_changed_cat, relevant_leaves_unchanged_cat, fill = TRUE)
    relevant_leaves_cat_list <- relevant_leaves_cat[, .(f_idx = .(f_idx)), by = c_idx]
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
      f_idx = Map(function(f_idx, min, max, i.min, i.max) {
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
        intersect(rel_min, rel_max)
      }, f_idx = f_idx, min = min, max = max, i.min = i.min, i.max = i.max))]

    # or-combine different conditions on the same feature
    or_within_row_cnt <- uniqueN(cnt_relevant, by = c("c_idx", "variable")) != nrow(cnt_relevant)
    if (or_within_row_cnt) {
      cnt_relevant <- cnt_relevant[, .(f_idx = .(Reduce(union, f_idx))), by = .(c_idx, variable)]
    }

    # Determine matching leaves for cnt conditions
    relevant_leaves_changed_cnt <- cnt_relevant[, Reduce(intersect, f_idx), by = c_idx][, .(c_idx, f_idx = V1)]
    conditions_unchanged_cnt <- setdiff(condition_long_step[, c_idx], cnt_conds[, c_idx])
    relevant_leaves_unchanged_cnt <- data.table(c_idx = rep(conditions_unchanged_cnt, each = nrow(forest)), f_idx = rep(forest[, f_idx], length(conditions_unchanged_cnt)))
    relevant_leaves_cnt <- rbind(relevant_leaves_changed_cnt, relevant_leaves_unchanged_cnt)

    # Calculate updates for cnt params matching cnt conditions
    cnt_new <- merge(merge(relevant_leaves_cnt, cnt_conds, by = "c_idx", allow.cartesian = TRUE, sort = FALSE), cnt, by = c("f_idx", "variable"), all.x = TRUE, allow.cartesian = TRUE, sort = FALSE)
    cnt_new[!is.na(val), `:=`(min = min.y,
                              max = max.y)]
    cnt_new[is.na(val), `:=`(min = pmax(min.x, min.y, na.rm = TRUE),
                             max = pmin(max.x, max.y, na.rm = TRUE))]
    cnt_new <- cnt_new[min <= max & sum(max.x, val, na.rm = TRUE) != min.y, ]
    cnt_new[, prob := NA_real_]
    if (family == "truncnorm") {
      cnt_new[!is.na(val), prob := dtruncnorm(val, a = min.y, b = max.y, mean = mu, sd = sigma) * (val != min.y)]
      cnt_new[is.na(val) & (min == min.y) & (max == max.y), prob := 1]
      cnt_new[is.na(val) & is.na(prob), prob := ptruncnorm(max, a = min.y, b = max.y, mean = mu, sd = sigma) - ptruncnorm(min, a = min.y, max.y, mean = mu, sd = sigma)]
    } else if (family == "unif") {
      cnt_new[!is.na(val), prob := dunif(val, min = min.y, max = max.y) * (val != min.y)]
      cnt_new[is.na(val) & (min == min.y) & (max == max.y), prob := 1]
      cnt_new[is.na(val) & is.na(prob), prob := punif(max, min = min.y, max = max.y) - punif(min, min = min.y, max.y)]
    }
    cnt_new[, c("min.x", "max.x", "min.y", "max.y") := NULL]

    # If or-combined cnt condition within rows exist, calculate likelihoods for ranges within leaves and norm to 1
    if (or_within_row_cnt) {
      cnt_new[, cvg_factor := sum(prob), by = .(f_idx, c_idx, variable)]
      cnt_new[, prob := prob / cvg_factor]
    } else {
      cnt_new[, `:=`(cvg_factor = prob, prob = 1)]
    }

    # Calculate final set of matching leaves
    if (nrow(cat_conds) > 0) {
      relevant_leaves <- merge(relevant_leaves_cnt, relevant_leaves_cat, by = c("c_idx", "f_idx"))[, .(c_idx, f_idx)]
    } else {
      relevant_leaves <- relevant_leaves_cnt[, .(c_idx, f_idx)]
    }

    # If no cnt conditions exist, output empty update table cnt_new for cnt params
  } else {
    relevant_leaves <- relevant_leaves_cat[, .(c_idx, f_idx)]
    cnt_new <- cbind(cnt[FALSE, ], data.table(cvg_factor = numeric(), c_idx = integer(), val = numeric(), prob = numeric()))
  }

  # Calculate updates for cat params matching cat conditions
  cat_new <- merge(merge(relevant_leaves, cat_conds, by = "c_idx", allow.cartesian = TRUE), cat, by = c("f_idx", "variable", "val"))

  # Ensure probabilities sum to 1
  cat_new[, cvg_factor := sum(prob), by = .(f_idx, c_idx, variable)]
  cat_new[, prob := prob / cvg_factor]

  list(cnt_new = cnt_new, cat_new = cat_new, relevant_leaves = relevant_leaves)
}
