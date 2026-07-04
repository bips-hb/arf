#' Generate synthetic data by sampling from the leaves of a random forest
#'
#' Draws synthetic samples by sampling, for each observation, a leaf from the
#' forest and then drawing each feature independently (marginally) from the real
#' observations that fall into that leaf. This is the intra-leaf sampling step
#' used internally by \code{\link{adversarial_rf}} to generate synthetic data
#' during the adversarial loop, exposed here as a stand-alone function.
#'
#' @param arf A trained ARF, as returned by \code{\link{adversarial_rf}} (a
#'   \code{ranger} object).
#' @param x_real Data whose intra-leaf structure is used for sampling, typically
#'   the data the forest was trained on.
#' @param params Optional circuit parameters as returned by \code{\link{forde}}.
#'   If supplied, the synthetic data is post-processed with the same routine used
#'   by \code{\link{forge}}: variable types and factor levels are restored,
#'   continuous variables are rounded to their observed precision (see
#'   \code{round}), and the class of the original input is reinstated. If
#'   \code{NULL}, a minimally processed \code{data.table} is returned with factor
#'   columns encoded as character, matching the representation used internally by
#'   \code{adversarial_rf}.
#' @param round Round continuous variables to their maximum precision in the real
#'   data? Only relevant when \code{params} is supplied.
#' @param factor_cols Optional logical vector flagging the factor columns of
#'   \code{x_real}. Computed from \code{x_real} if \code{NULL}. Mainly for
#'   internal use to avoid recomputation.
#' @param lvls Optional list of factor levels for the factor columns of
#'   \code{x_real}. Computed from \code{x_real} if \code{NULL}. Mainly for
#'   internal use.
#' @param prep Prepare \code{x_real} with the internal pre-processing routine
#'   before sampling? Set to \code{FALSE} if \code{x_real} is already prepared
#'   (internal use).
#'
#' @return A dataset of \code{nrow(x_real)} synthetic samples. When \code{params}
#'   is supplied, its class and column types match the original data; otherwise a
#'   \code{data.table} with factor columns encoded as character.
#'
#' @references
#' Watson, D., Blesch, K., Kapar, J., & Wright, M. (2023). Adversarial random
#' forests for density estimation and generative modeling. In \emph{Proceedings
#' of the 26th International Conference on Artificial Intelligence and
#' Statistics}, pp. 5357-5375.
#'
#' @seealso
#' \code{\link{adversarial_rf}}, \code{\link{forde}}, \code{\link{forge}}
#'
#' @examples
#' arf <- adversarial_rf(iris)
#'
#' # Minimally processed output (factors as character)
#' x_synth <- sample_from_leaves(arf, iris)
#'
#' # Fully post-processed output, consistent with forge()
#' psi <- forde(arf, iris)
#' x_synth <- sample_from_leaves(arf, iris, params = psi)
#'
#' @export
#'
sample_from_leaves <- function(arf, x_real, params = NULL, round = TRUE,
                               factor_cols = NULL, lvls = NULL, prep = TRUE) {
  # To avoid data.table check issues
  i <- b <- cnt <- obs <- tree <- leaf <- N <- . <- NULL
  # Prep data
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
  # Assemble synthetic data
  x_synth <- as.data.table(split(values_drawn_stacked, by = 'col_id', keep.by = FALSE))
  setnames(x_synth, names(x_real))
  if (any(factor_cols)) {
    x_synth[, names(which(factor_cols))] <- lapply(names(which(factor_cols)), function(j) {
      lvls[[j]][x_synth[[j]]]
    })
  }
  # Optionally post-process with the same routine used by forge()
  if (!is.null(params)) {
    x_synth <- post_x(x_synth, params, round = round)
  }
  x_synth
}
