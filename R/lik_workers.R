#' @keywords internal
#' @noRd

# Per-fold worker for lik(), extracted so foreach and mirai backends share one
# definition (params/x/preds/omega mori-shared under mirai). `has_arf` replaces
# the closed-over arf object (only used for is.null() branching; preds carries
# the arf-derived leaf assignments). Uses bare data.table verbs, so mirai
# daemons must have arf loaded (arf imports data.table).
arf_lik_fold <- function(fold, params, x, factor_cols, leaves, omega, preds,
                         batch_idx, pure, has_arf) {
  # To avoid data.table check issues
  tree <- cvg <- leaf <- variable <- mu <- sigma <- value <- obs <- prob <-
    V1 <- relation <- f_idx <- wt <- val <- family <- f_idx_uncond <- . <-
    lik <- s_idx <- min <- max <- NULL

  # Prep work
  psi_cnt <- psi_cat <- NULL
  if (!has_arf & !isTRUE(pure)) {
    omega_tmp <- rbindlist(lapply(batch_idx[[fold]], function(i) {
      omega$obs <- i
      omega$wt <- NULL
      return(omega)
    }))
  }

  # Continuous data
  if (any(!factor_cols)) {
    fam <- params$meta[class == 'numeric', unique(family)]
    x_long <- melt(
      data.table(obs = batch_idx[[fold]],
                 x[batch_idx[[fold]], !factor_cols, drop = FALSE]),
      id.vars = 'obs', variable.factor = FALSE
    )
    if (!has_arf) {
      psi_cnt <- merge(params$cnt[f_idx %in% leaves], x_long, by = 'variable',
                       sort = FALSE, allow.cartesian = TRUE)
      rm(x_long)
    } else {
      preds_cnt <- merge(preds[f_idx %in% leaves], x_long, by = 'obs',
                         sort = FALSE, allow.cartesian = TRUE)
      rm(x_long)
      psi_cnt <- merge(params$cnt[f_idx %in% leaves], preds_cnt,
                       by = c('f_idx', 'variable'), sort = FALSE)
      rm(preds_cnt)
    }
    if (fam == 'truncnorm') {
      psi_cnt[, lik := truncnorm::dtruncnorm(value, a = min, b = max,
                                             mean = mu, sd = sigma)]
    } else if (fam == 'unif') {
      psi_cnt[, lik := stats::dunif(value, min = min, max = max)]
    }
    psi_cnt[value == min, lik := 0]
    psi_cnt[, lik := prod(lik), by = .(f_idx, obs)]
    psi_cnt <- unique(psi_cnt[lik > 0, .(f_idx, obs, lik)])
    if (!has_arf & !isTRUE(pure)) {
      omega_tmp <- merge(omega_tmp, psi_cnt[, .(f_idx, obs)],
                         by = c('f_idx', 'obs'), sort = FALSE)
      leaves <- omega_tmp[, unique(f_idx)]
    }
  }

  # Categorical data
  if (any(factor_cols)) {
    x_tmp <- x[batch_idx[[fold]], factor_cols, drop = FALSE]
    n_tmp <- nrow(x_tmp)
    x_long <- melt(
      data.table(obs = batch_idx[[fold]], x_tmp),
      id.vars = 'obs', value.name = 'val', variable.factor = FALSE
    )
    # Speedups are possible if there are many duplicates
    is_unique <- !duplicated(x_tmp)
    if (all(is_unique)) {
      x_unique <- x_long
      colnames(x_unique)[1] <- 's_idx'
    } else {
      x_unique <- unique(x_tmp)
      x_unique <- melt(
        data.table(s_idx = seq_len(nrow(x_unique)), x_unique),
        id.vars = 's_idx', value.name = 'val', variable.factor = FALSE
      )
      s_idx <- integer(length = n_tmp)
      s_idx[is_unique] <- seq_len(sum(is_unique))
      for (i in 2:n_tmp) {
        if (s_idx[i] == 0L) {
          s_idx[i] <- s_idx[i - 1L]
        }
      }
      idx_dt <- data.table(obs = batch_idx[[fold]], s_idx = s_idx)
    }
    if (!has_arf) {
      grd <- rbindlist(lapply(which(factor_cols), function(j) {
        expand.grid('f_idx' = leaves, 'variable' = colnames(x)[j],
                    'val' = x_long[variable == colnames(x)[j], unique(val)],
                    stringsAsFactors = FALSE)
      }))
      rm(x_long)
      psi_cat <- merge(params$cat[f_idx %in% leaves], grd,
                       by = c('f_idx', 'variable', 'val'),
                       sort = FALSE, all.y = TRUE)
      rm(grd)
      psi_cat[is.na(prob), prob := 0]
      psi_cat <- merge(psi_cat, x_unique, by = c('variable', 'val'),
                       sort = FALSE, allow.cartesian = TRUE)
      psi_cat[, lik := prod(prob), by = .(f_idx, s_idx)]
      psi_cat <- unique(psi_cat[lik > 0, .(f_idx, s_idx, lik)])
      if (all(is_unique)) {
        setnames(psi_cat, 's_idx', 'obs')
      } else {
        if (!isTRUE(pure)) {
          omega_tmp <- merge(idx_dt, omega_tmp, by = 'obs', sort = FALSE)
          psi_cat <- merge(psi_cat, omega_tmp, by = c('f_idx', 's_idx'),
                           sort = FALSE)[, s_idx := NULL]
          rm(omega_tmp)
          setcolorder(psi_cat, c('f_idx', 'obs', 'lik'))
          psi_cnt <- merge(psi_cnt, psi_cat[, .(f_idx, obs)],
                           by = c('f_idx', 'obs'), sort = FALSE)
        }
      }
    } else {
      preds_cat <- merge(preds[f_idx %in% leaves], x_long, by = 'obs',
                         sort = FALSE, allow.cartesian = TRUE)
      rm(x_long)
      psi_cat <- merge(params$cat, preds_cat, by = c('f_idx', 'variable', 'val'),
                       sort = FALSE, allow.cartesian = TRUE, all.y = TRUE)
      rm(preds_cat)
      psi_cat[is.na(prob), prob := 0]
      psi_cat[, lik := prod(prob), by = .(f_idx, obs)]
      psi_cat <- unique(psi_cat[lik > 0, .(f_idx, obs, lik)])
    }
  }

  # Put it together
  psi_x <- rbind(psi_cnt, psi_cat)
  if (!isTRUE(pure)) {
    psi_x <- psi_x[, prod(lik), by = .(f_idx, obs)]
    setnames(psi_x, 'V1', 'lik')
  }

  # Reduce to per-observation log-likelihoods here rather than on the calling
  # process: folds cover disjoint obs and omega is a worker argument, so the
  # reduction is fold-local. This shrinks the returned object from one row per
  # (obs, leaf) to one per obs (less to serialize back under mirai) and
  # parallelizes what used to be a serial post-pass over every obs x leaf pair.
  psi_x <- merge(psi_x, omega, by = 'f_idx', sort = FALSE)
  psi_x <- psi_x[, log(crossprod(wt, lik)), by = obs]
  setnames(psi_x, 'V1', 'lik')
  psi_x
}
