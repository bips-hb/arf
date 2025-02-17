
#' Generate synthetic data by sampling from leaves of a random forest
#'
#' @param arf Random forest to sample from
#' @param x_real Data to sample from
#' @param factor_cols Logical vector indicating which columns are factors (optional)
#' @param lvls List of factor levels (optional)
#'
#' @returns A data.table of synthetic data
#' @export
#'
#' @examples
#' arf <- adversarial_rf(iris)
#' sample_from_leaves(arf, iris)
sample_from_leaves <- function(arf, x_real, factor_cols = NULL, lvls = NULL, prep = TRUE) {
  if (prep) {
    x_real <- prep_x(x_real, verbose = FALSE)
  }
  n <- nrow(x_real)
  d <- ncol(x_real)
  if (is.null(factor_cols)) {
    factor_cols <- sapply(x_real, is.factor)
  }
  if (is.null(lvls)) {
    lvls <- lapply(x_real[factor_cols], levels)
  }
  # Sample leaves and get values from other observations in the same leaf
  nodeIDs <- stats::predict(arf, x_real, type = 'terminalNodes')$predictions
  tmp <- data.table('tree' = rep(seq_len(arf$num.trees), each = n), 
                    'leaf' = as.integer(nodeIDs))
  tmp2 <- tmp[sample(.N, n, replace = TRUE)]
  tmp2 <- unique(tmp2[, cnt := .N, by = .(tree, leaf)])
  draw_from <- rbindlist(lapply(seq_len(arf$num.trees), function(b) {
    x_real_b <- cbind(x_real, tmp[tree == b])
    x_real_b[, factor_cols] <- lapply(x_real_b[, factor_cols, drop = FALSE], as.numeric)
    merge(tmp2, x_real_b, by = c('tree', 'leaf'), 
          sort = FALSE)[, N := .N, by = .(tree, leaf)]
  }))
  # Draw new observations by sampling marginally from those leaves
  draw_params_within <- unique(draw_from, by = c('tree','leaf'))[, .(cnt, N)]
  adj_absolut_col <- rep(c(0, draw_params_within[-.N, cumsum(N)]), 
                         times = draw_params_within$cnt)
  adj_absolut <- rep(adj_absolut_col, d) + rep(seq(0, d - 1) * nrow(draw_from), each = n)
  idx_drawn_within <- ceiling(runif(n * d, 0, rep(draw_params_within$N, draw_params_within$cnt)))
  idx_drawn <- idx_drawn_within + adj_absolut
  draw_from_stacked <- unlist(draw_from[, -c('tree', 'leaf', 'cnt', 'N')], 
                              use.names = FALSE)
  values_drawn_stacked <- data.table('col_id' = rep(seq_len(d), each = n), 
                                     'values' = draw_from_stacked[idx_drawn])
  # Return synthetic data
  x_synth <- as.data.table(split(values_drawn_stacked, by = 'col_id', keep.by = FALSE))
  setnames(x_synth, names(x_real))
  if (any(factor_cols)) {
    x_synth[, names(which(factor_cols))] <- lapply(names(which(factor_cols)), function(j) {
      lvls[[j]][x_synth[[j]]]
    })
  }
  x_synth
}