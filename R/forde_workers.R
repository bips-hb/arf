#' @keywords internal
#' @noRd

# Per-tree worker bodies extracted from forde() so foreach and mirai
# backends share one definition. All dependencies are explicit
# arguments; nothing captured from enclosing scope. This is what makes
# the mirai backend possible — daemons run in clean environments.

# Worker 1: compute leaf bounds for one tree.
arf_bnd_fn <- function(tree, forest, d, finite_bounds, factor_cols,
                       x, epsilon, colnames_x) {
  # data.table NSE silencing
  variable <- NULL  # nolint

  num_nodes <- length(forest$split.varIDs[[tree]])
  lb <- matrix(-Inf, nrow = num_nodes, ncol = d)
  ub <- matrix(Inf, nrow = num_nodes, ncol = d)
  if (finite_bounds == "global" && any(!factor_cols)) {
    for (j in which(!factor_cols)) {
      min_j <- min(x[[j]], na.rm = TRUE)
      max_j <- max(x[[j]], na.rm = TRUE)
      gap <- max_j - min_j
      lb[, j] <- min_j - epsilon / 2 * gap
      ub[, j] <- max_j + epsilon / 2 * gap
    }
  }
  for (i in seq_len(num_nodes)) {
    left_child <- forest$child.nodeIDs[[tree]][[1]][i] + 1L
    right_child <- forest$child.nodeIDs[[tree]][[2]][i] + 1L
    splitvarID <- forest$split.varIDs[[tree]][i] + 1L
    splitval <- forest$split.values[[tree]][i]
    if (left_child > 1) {
      ub[left_child, ] <- ub[right_child, ] <- ub[i, ]
      lb[left_child, ] <- lb[right_child, ] <- lb[i, ]
      if (left_child != right_child) {
        ub[left_child, splitvarID] <- lb[right_child, splitvarID] <- splitval
      }
    }
  }
  leaves <- which(forest$child.nodeIDs[[tree]][[1]] == 0L)
  colnames(lb) <- colnames(ub) <- colnames_x
  data.table::merge.data.table(
    data.table::melt(
      data.table::data.table(tree = tree, leaf = leaves,
                             lb[leaves, , drop = FALSE]),
      id.vars = c("tree", "leaf"), value.name = "min"),
    data.table::melt(
      data.table::data.table(tree = tree, leaf = leaves,
                             ub[leaves, , drop = FALSE]),
      id.vars = c("tree", "leaf"), value.name = "max"),
    by = c("tree", "leaf", "variable"), sort = FALSE)
}

# Worker 2: compute continuous-variable distribution params for one tree.
arf_psi_cnt_fn <- function(tree, x, factor_cols, pred, inbag.counts, n,
                           oob, bnds, finite_bounds, epsilon, family) {
  # data.table NSE silencing
  leaf <- variable <- value <- min_emp <- max_emp <- length_emp <- mu <-
    sigma <- NA_share <- new_min <- new_max <- mid <- sigma0 <- f_idx <- . <-
    NULL  # nolint

  dt <- data.table::data.table(x[, !factor_cols, drop = FALSE],
                               leaf = pred[, tree])
  if (isTRUE(oob)) {
    dt <- dt[inbag.counts[[tree]][1:n] == 0L, ]
    dt <- dt[!is.na(leaf)]
  } else if (identical(oob, "inbag")) {
    dt <- dt[inbag.counts[[tree]][1:n] > 0L, ]
    dt <- dt[!is.na(leaf)]
  }
  dt <- data.table::melt(dt, id.vars = "leaf",
                         variable.factor = FALSE)[, tree := tree]
  dt <- data.table::merge.data.table(
    dt,
    bnds[, .(tree, leaf, variable, min, max, f_idx)],
    by = c("tree", "leaf", "variable"), sort = FALSE)
  if (finite_bounds == "local") {
    dt[, c("min_emp", "max_emp") := .(min(value, na.rm = TRUE),
                                       max(value, na.rm = TRUE)),
       by = .(leaf, variable)]
    dt[, length_emp := max_emp - min_emp]
    length_emp_0_replace <- min(
      dt[length_emp > 0, min(length_emp, na.rm = TRUE)],
      max(epsilon, 1e-12))
    dt[length_emp == 0,
       c("min_emp", "max_emp", "length_emp") :=
         .(min_emp - length_emp_0_replace / 2,
           max_emp + length_emp_0_replace / 2,
           length_emp_0_replace)]
    dt[, c("min", "max", "min_emp", "max_emp", "length_emp") :=
         .(data.table::fifelse(
             !is.finite(min) & !is.na(min_emp),
             min_emp - length_emp * (epsilon / 2), min),
           data.table::fifelse(
             !is.finite(max) & !is.na(max_emp),
             max_emp + length_emp * (epsilon / 2), max),
           NULL, NULL, NULL)]
  }
  if (family == "truncnorm") {
    dt[, c("mu", "sigma", "NA_share") :=
         .(mean(value, na.rm = TRUE),
           stats::sd(value, na.rm = TRUE),
           sum(is.na(value)) / .N),
       by = .(leaf, variable)]
    dt[, c("min_emp", "max_emp") :=
         .(min(value, na.rm = TRUE), max(value, na.rm = TRUE)),
       by = variable]
    dt[NA_share == 1,
       c("min", "max") :=
         .(data.table::fifelse(is.infinite(min), min_emp, min),
           data.table::fifelse(is.infinite(max), max_emp, max))]
    dt[, c("min_emp", "max_emp") := NULL]
    dt[NA_share == 1, mu := (max + min) / 2]
    dt[is.na(sigma), sigma := 0]
    if (any(dt[, sigma == 0])) {
      dt[, new_min := data.table::fifelse(
              !is.finite(min), min(value, na.rm = TRUE), min),
         by = variable]
      dt[, new_max := data.table::fifelse(
              !is.finite(max), max(value, na.rm = TRUE), max),
         by = variable]
      dt[, mid := (new_min + new_max) / 2]
      dt[, sigma0 := (new_max - mid) / stats::qnorm(0.975)]
      # Bayesian prior: 95% mass within bounding box, df = 2.
      # See forde.R history for derivation.
      dt[sigma == 0, sigma := sqrt(2 / .N * sigma0^2),
         by = .(variable, leaf)]
      dt[, c("new_min", "new_max", "mid", "sigma0") := NULL]
    }
  } else if (family == "unif") {
    dt[, NA_share := sum(is.na(value)) / .N, by = .(leaf, variable)]
  }
  unique(dt[, c("tree", "leaf", "value") := NULL])
}

# Worker 3: compute categorical-variable distribution params for one tree.
arf_psi_cat_fn <- function(tree, x, factor_cols, pred, bnds, oob,
                           lvl_df_rf, alpha) {
  # data.table NSE silencing
  leaf <- variable <- val <- NA_share <- count <- val_count <- k <-
    level <- prob <- f_idx <- . <- NULL  # nolint

  dt <- data.table::data.table(x[, factor_cols, drop = FALSE],
                               leaf = pred[, tree])
  if (isTRUE(oob)) {
    dt <- dt[!is.na(leaf)]
  }
  dt <- data.table::melt(dt, id.vars = "leaf",
                         variable.factor = FALSE,
                         value.factor = FALSE,
                         value.name = "val")[, tree := tree]
  dt[, NA_share := sum(is.na(val)) / .N, by = .(leaf, variable)]
  dt <- dt[!(is.na(val) & NA_share != 1)]
  if (dt[, any(NA_share == 1)]) {
    all_na <- unique(dt[NA_share == 1, ])
    dt <- dt[NA_share != 1, ]
    all_na <- data.table::merge.data.table(
      all_na, bnds[, .(tree, leaf, variable, min, max, f_idx)],
      by = c("tree", "leaf", "variable"), sort = FALSE)
    all_na[!is.finite(min), min := 0.5]
    for (j in names(which(factor_cols))) {
      all_na[!is.finite(max) & variable == j,
             max := lvl_df_rf[variable == j, max(level)]]
    }
    all_na[!grepl("\\.5", min), min := min + 0.5]
    all_na[!grepl("\\.5", max), max := max + 0.5]
    all_na[, min := min + 0.5][, max := max - 0.5]
    all_na <- all_na[, .(level = seq(min, max), NA_share),
                     by = .(leaf, variable)]
    all_na <- data.table::merge.data.table(
      all_na, lvl_df_rf, by = c("variable", "level"))
    all_na[, level := NULL][, tree := tree]
    data.table::setcolorder(all_na, colnames(dt))
    dt <- rbind(dt, all_na)
  }
  dt[, count := .N, by = .(leaf, variable)]
  dt <- data.table::merge.data.table(
    dt, bnds[, .(tree, leaf, variable, min, max, f_idx)],
    by = c("tree", "leaf", "variable"), sort = FALSE)
  dt[, c("tree", "leaf") := NULL]
  if (alpha == 0) {
    dt <- unique(dt[, prob := .N / count, by = .(f_idx, variable, val)])
  } else {
    dt <- unique(dt[, val_count := .N, by = .(f_idx, variable, val)])
    dt <- data.table::merge.data.table(
      dt, lvl_df_rf[, .(k = .N), by = variable], by = "variable")
    dt[!is.finite(min), min := 0.5][!is.finite(max), max := k + 0.5]
    dt[!grepl("\\.5", min), min := min + 0.5][
      !grepl("\\.5", max), max := max + 0.5]
    dt[, k := max - min]
    tmp <- dt[, seq(min[1] + 0.5, max[1] - 0.5),
              by = .(f_idx, variable)]
    data.table::setnames(tmp, "V1", "level")
    tmp <- data.table::merge.data.table(
      tmp, lvl_df_rf, by = c("variable", "level"),
      sort = FALSE)[, level := NULL]
    tmp <- data.table::merge.data.table(
      tmp, unique(dt[, .(f_idx, variable, count, k)]),
      by = c("f_idx", "variable"), sort = FALSE)
    dt <- data.table::merge.data.table(
      tmp, dt, by = c("f_idx", "variable", "val", "count", "k"),
      all.x = TRUE, sort = FALSE)
    dt[is.na(val_count), val_count := 0]
    dt[, NA_share := mean(NA_share, na.rm = TRUE),
       by = .(f_idx, variable)]
    dt[, prob := (val_count + alpha) / (count + alpha * k),
       by = .(f_idx, variable, val)]
    dt[, c("val_count", "k") := NULL]
  }
  dt[, c("count", "min", "max") := NULL]
  data.table::setcolorder(dt, c("f_idx", "variable", "val", "prob",
                                 "NA_share"))
  dt
}
