#' @keywords internal
#' @noRd

# Difference of two values represented on the log scale, where log_x >= log_y.
arf_logspace_sub <- function(log_x, log_y) {
  x <- log_y - log_x
  out <- rep(-Inf, length(x))
  positive <- log_y < log_x
  small <- positive & x < -log(2)
  large <- positive & !small
  out[small] <- log1p(-exp(x[small]))
  out[large] <- log(-expm1(x[large]))
  log_x + out
}

# Log probability that a normal variable falls between lower and upper.
arf_log_pnorm_diff <- function(lower, upper, mean, sd) {
  lower_z <- (lower - mean) / sd
  upper_z <- (upper - mean) / sd
  out <- numeric(length(lower_z))
  left <- upper_z <= 0
  right <- !left & lower_z >= 0
  middle <- !left & !right

  if (any(left)) {
    out[left] <- arf_logspace_sub(
      stats::pnorm(upper_z[left], log.p = TRUE),
      stats::pnorm(lower_z[left], log.p = TRUE)
    )
  }
  if (any(right)) {
    out[right] <- arf_logspace_sub(
      stats::pnorm(lower_z[right], lower.tail = FALSE, log.p = TRUE),
      stats::pnorm(upper_z[right], lower.tail = FALSE, log.p = TRUE)
    )
  }
  if (any(middle)) {
    out[middle] <- arf_logspace_sub(
      stats::pnorm(upper_z[middle], log.p = TRUE),
      stats::pnorm(lower_z[middle], log.p = TRUE)
    )
  }
  out[upper_z <= lower_z] <- -Inf
  out
}

# Truncated-normal log density without truncnorm's probability-underflow fallback.
arf_log_dtruncnorm <- function(x, a, b, mean, sd) {
  out <- stats::dnorm(x, mean = mean, sd = sd, log = TRUE) -
    arf_log_pnorm_diff(a, b, mean, sd)
  outside <- !is.na(x) & (x < a | x > b)
  out[outside] <- -Inf
  out
}

# Log of a sum of positive values represented on the log scale.
arf_log_sum_exp <- function(x) {
  max_x <- max(x)
  if (!is.finite(max_x)) {
    return(max_x)
  }
  max_x + log(sum(exp(x - max_x)))
}

# Per-fold worker for lik(), extracted so foreach and mirai backends share one
# definition (params/x/preds/omega mori-shared under mirai). `has_arf` replaces
# the closed-over arf object (only used for is.null() branching; preds carries
# the arf-derived leaf assignments). Uses bare data.table verbs, so mirai
# daemons must have arf loaded (arf imports data.table).
arf_lik_fold <- function(fold, params, x, factor_cols, leaves, omega, preds, batch_idx, pure, has_arf) {
  # To avoid data.table check issues
  tree <- cvg <- leaf <- variable <- mu <- sigma <- value <- obs <- prob <-
    V1 <- relation <- f_idx <- wt <- val <- family <- f_idx_uncond <- . <-
      lik <- lik_cnt <- lik_cat <- s_idx <- min <- max <- NULL

  psi_cnt <- psi_cat <- NULL

  # Continuous data
  if (!all(factor_cols)) {
    fam <- params$meta[class == 'numeric', unique(family)]
    x_long <- melt(
      data.table(obs = batch_idx[[fold]], x[batch_idx[[fold]], !factor_cols, drop = FALSE]),
      id.vars = 'obs',
      variable.factor = FALSE
    )
    if (!has_arf) {
      psi_cnt <- merge(params$cnt[f_idx %in% leaves], x_long, by = 'variable', sort = FALSE, allow.cartesian = TRUE)
      rm(x_long)
    } else {
      preds_cnt <- merge(preds[f_idx %in% leaves], x_long, by = 'obs', sort = FALSE, allow.cartesian = TRUE)
      rm(x_long)
      psi_cnt <- merge(params$cnt[f_idx %in% leaves], preds_cnt, by = c('f_idx', 'variable'), sort = FALSE)
      rm(preds_cnt)
    }
    if (fam == 'truncnorm') {
      psi_cnt[, lik := arf_log_dtruncnorm(value, a = min, b = max, mean = mu, sd = sigma)]
    } else if (fam == 'unif') {
      psi_cnt[, lik := stats::dunif(value, min = min, max = max, log = TRUE)]
    }
    psi_cnt[value == min, lik := -Inf]
    psi_cnt[, lik := sum(lik), by = .(f_idx, obs)]
    psi_cnt <- unique(psi_cnt[!is.na(lik) & lik != -Inf, .(f_idx, obs, lik)])
    if (!has_arf & !isTRUE(pure)) {
      # Leaves with zero continuous density need no categorical grid
      leaves <- psi_cnt[, unique(f_idx)]
    }
  }

  # Categorical data
  if (any(factor_cols)) {
    x_tmp <- x[batch_idx[[fold]], factor_cols, drop = FALSE]
    x_long <- melt(
      data.table(obs = batch_idx[[fold]], x_tmp),
      id.vars = 'obs',
      value.name = 'val',
      variable.factor = FALSE
    )
    # Speedups are possible if there are many duplicates
    is_unique <- !duplicated(x_tmp)
    if (all(is_unique)) {
      x_unique <- x_long
      colnames(x_unique)[1] <- 's_idx'
    } else {
      x_dt <- as.data.table(x_tmp)
      x_pattern <- unique(x_dt)
      x_unique <- melt(
        data.table(s_idx = seq_len(nrow(x_pattern)), x_pattern),
        id.vars = 's_idx',
        value.name = 'val',
        variable.factor = FALSE
      )
      # Each row maps to the pattern it matches, duplicates included
      idx_dt <- data.table(
        obs = batch_idx[[fold]],
        s_idx = x_pattern[x_dt, on = names(x_dt), which = TRUE]
      )
    }
    if (!has_arf) {
      grd <- rbindlist(lapply(which(factor_cols), function(j) {
        expand.grid(
          'f_idx' = leaves,
          'variable' = colnames(x)[j],
          'val' = x_long[variable == colnames(x)[j], unique(val)],
          stringsAsFactors = FALSE
        )
      }))
      rm(x_long)
      psi_cat <- merge(
        params$cat[f_idx %in% leaves],
        grd,
        by = c('f_idx', 'variable', 'val'),
        sort = FALSE,
        all.y = TRUE
      )
      rm(grd)
      psi_cat[is.na(prob), prob := 0]
      psi_cat <- merge(psi_cat, x_unique, by = c('variable', 'val'), sort = FALSE, allow.cartesian = TRUE)
      psi_cat[, lik := sum(log(prob)), by = .(f_idx, s_idx)]
      psi_cat <- unique(psi_cat[!is.na(lik) & lik != -Inf, .(f_idx, s_idx, lik)])
      if (all(is_unique)) {
        setnames(psi_cat, 's_idx', 'obs')
      } else {
        psi_cat <- merge(psi_cat, idx_dt, by = 's_idx', sort = FALSE, allow.cartesian = TRUE)[, s_idx := NULL]
        setcolorder(psi_cat, c('f_idx', 'obs', 'lik'))
      }
    } else {
      preds_cat <- merge(preds[f_idx %in% leaves], x_long, by = 'obs', sort = FALSE, allow.cartesian = TRUE)
      rm(x_long)
      psi_cat <- merge(
        params$cat,
        preds_cat,
        by = c('f_idx', 'variable', 'val'),
        sort = FALSE,
        allow.cartesian = TRUE,
        all.y = TRUE
      )
      rm(preds_cat)
      psi_cat[is.na(prob), prob := 0]
      psi_cat[, lik := sum(log(prob)), by = .(f_idx, obs)]
      psi_cat <- unique(psi_cat[!is.na(lik) & lik != -Inf, .(f_idx, obs, lik)])
    }
  }

  # Put it together. Both blocks are filtered to lik > 0, so a leaf where one
  # block is zero is missing from that table: inner join, else its product is
  # the other block's density alone.
  if (isTRUE(pure)) {
    psi_x <- rbind(psi_cnt, psi_cat)
  } else {
    psi_x <- merge(psi_cnt, psi_cat, by = c('f_idx', 'obs'), sort = FALSE, suffixes = c('_cnt', '_cat'))
    psi_x <- psi_x[, .(f_idx, obs, lik = lik_cnt + lik_cat)]
  }

  # Reduce to per-observation log-likelihoods here rather than on the calling
  # process: folds cover disjoint obs and omega is a worker argument, so the
  # reduction is fold-local. This shrinks the returned object from one row per
  # (obs, leaf) to one per obs (less to serialize back under mirai) and
  # parallelizes what used to be a serial post-pass over every obs x leaf pair.
  psi_x <- merge(psi_x, omega, by = 'f_idx', sort = FALSE)
  psi_x[, lik := log(wt) + lik]
  psi_x <- psi_x[, .(lik = arf_log_sum_exp(lik)), by = obs]
  psi_x
}
