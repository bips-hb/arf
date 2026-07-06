#' @keywords internal
#' @noRd

# Per-step worker for forge(), extracted so the foreach and mirai backends share
# one definition. All dependencies are explicit arguments (params/evidence are
# mori-shared read-only under the mirai backend). Unlike the forde workers this
# calls arf internals (cforde, resample, post_x), so mirai daemons must have arf
# loaded (see the everywhere() call in forge()).
arf_forge_step <- function(step_, params, evidence, n_synth, factor_cols,
                           evidence_row_mode, nomatch, verbose, round,
                           sample_NAs, stepsize, stepsize_cforde,
                           parallel_cforde) {
  # data.table NSE silencing
  tree <- cvg <- leaf <- idx <- family <- mu <- sigma <- prob <- dat <-
    variable <- relation <- wt <- j <- f_idx <- val <- . <- step_x <- c_idx <-
    f_idx_uncond <- N <- I <- V1 <- min <- max <- NA_share <- NULL

  # Prepare the event space
  if (is.null(evidence) || (ncol(evidence) == 2 && all(colnames(evidence) == c("f_idx", "wt")))) {
    cparams <- NULL
  } else {
    index_start <- (step_ - 1) * stepsize + 1
    index_end <- min(step_ * stepsize, nrow(evidence))
    evidence_part <- evidence[index_start:index_end, ]
    cparams <- cforde(params, evidence_part, evidence_row_mode, nomatch, verbose,
                      stepsize_cforde, parallel_cforde)
    if (is.null(cparams)) {
      n_synth <- n_synth * nrow(evidence_part)
    }
  }

  if (is.null(cparams)) {
    if (is.null(evidence)) {
      num_trees <- params$forest[, max(tree)]
      omega <- params$forest[, .(f_idx, f_idx_uncond = f_idx, cvg)]
      omega[, `:=`(c_idx = 1, wt = cvg / num_trees)]
      omega[, cvg := NULL]
    } else {
      omega <- data.table::copy(evidence)
      omega[, f_idx_uncond := f_idx]
      omega[, c_idx := 1]
    }
  } else {
    omega <- cparams$forest[, .(c_idx, f_idx, f_idx_uncond, wt = cvg)]
  }
  omega <- omega[wt > 0, ]

  if (nrow(omega) == 1) {
    omega <- omega[rep(1, n_synth), ][, idx := .I]
  } else {
    if (evidence_row_mode == "or") {
      draws <- omega[, .(f_idx = resample(f_idx, size = n_synth, replace = TRUE, prob = wt))]
      omega <- merge(draws, omega, by = "f_idx", sort = FALSE)[, idx := .I]
    } else {
      draws <- omega[, .(f_idx = resample(f_idx, size = n_synth, replace = TRUE, prob = wt)), by = c_idx]
      omega <- merge(draws, omega, by = c("c_idx", "f_idx"), sort = FALSE)[, idx := .I]
    }
    data.table::setcolorder(omega, "idx")
  }

  # Simulate continuous data
  synth_cnt <- synth_cat <- NULL
  if (any(!factor_cols)) {
    fam <- params$meta[family != 'multinom', unique(family)]
    if (is.null(cparams)) {
      psi_cond <- data.table::data.table()
    } else {
      psi_cond <- merge(omega, cparams$cnt[, -c("cvg_factor", "f_idx_uncond")], by = c('c_idx', 'f_idx'),
                        sort = FALSE, allow.cartesian = TRUE)[prob > 0, ]
      if (any(psi_cond[, prob != 1])) {
        psi_cond[, I := .I]
        psi_cond <- psi_cond[sort(c(psi_cond[prob == 1, I],
                        psi_cond[prob > 0 & prob < 1, data.table::fifelse(.N > 1, resample(I, 1, prob = prob), 0), by = .(variable, idx)][, V1])), -"I"]
      }
      psi_cond[, prob := NULL]
    }
    psi <- unique(rbind(psi_cond,
                        merge(omega, params$cnt, by.x = 'f_idx_uncond', by.y = 'f_idx',
                              sort = FALSE, allow.cartesian = TRUE)[, val := NA_real_]),
                  by = c("idx", "variable"))
    if (fam == 'truncnorm') {
      psi[is.na(val), val := truncnorm::rtruncnorm(.N, a = min, b = max, mean = mu, sd = sigma)]
      psi[is.na(val), val := mu]
    } else if (fam == 'unif') {
      psi[is.na(val), val := stats::runif(.N, min = min, max = max)]
    }
    NA_share_cnt <- psi[, .(idx, variable, NA_share)]
    synth_cnt <- data.table::dcast(psi, idx ~ variable, value.var = 'val')[, idx := NULL]
  }

  # Simulate categorical data
  if (any(factor_cols)) {
    if (is.null(cparams)) {
      psi <- merge(omega, params$cat, by.x = 'f_idx_uncond', by.y = 'f_idx', sort = FALSE, allow.cartesian = TRUE)
    } else {
      psi_cond <- merge(omega, cparams$cat[, -c("cvg_factor", "f_idx_uncond")], by = c('c_idx', 'f_idx'),
                        sort = FALSE, allow.cartesian = TRUE)
      psi_uncond <- merge(omega, params$cat, by.x = 'f_idx_uncond', by.y = 'f_idx',
                          sort = FALSE, allow.cartesian = TRUE)
      psi_uncond_relevant <- psi_uncond[!psi_cond, on = .(idx, variable)]
      psi <- rbind(psi_cond, psi_uncond_relevant)
    }
    psi[prob < 1, val := sample(val, 1, prob = prob), by = .(variable, idx)]
    psi <- unique(psi[, .(idx, variable, val, NA_share)])
    NA_share_cat <- psi[, .(idx, variable, NA_share)]
    synth_cat <- data.table::dcast(psi, idx ~ variable, value.var = 'val')[, idx := NULL]
  }

  # Combine, optionally impose constraint(s)
  x_synth <- cbind(synth_cnt, synth_cat)
  if (length(x_synth) == 0) {
    x_synth <- evidence_part[FALSE, ]
  }
  x_synth <- post_x(x_synth, params, round)

  if (sample_NAs) {
    data.table::setDT(x_synth)
    NA_share <- rbind(NA_share_cnt, NA_share_cat)
    data.table::setorder(NA_share[, variable := factor(variable, levels = params$meta[, variable])], variable, idx)
    NA_share[, dat := stats::rbinom(.N, 1, prob = NA_share)]
    x_synth[data.table::dcast(NA_share, formula = idx ~ variable, value.var = "dat")[, -"idx"] == 1] <- NA
    x_synth <- post_x(x_synth, params, round)
  }
  if (evidence_row_mode == "separate" & any(omega[, is.na(f_idx)])) {
    data.table::setDT(x_synth)
    indices_na <- cparams$forest[is.na(f_idx), c_idx]
    indices_sampled <- cparams$forest[!is.na(f_idx), unique(c_idx)]
    rows_na <- data.table::dcast(rbind(data.table::data.table(c_idx = 0, variable = params$meta[, variable]),
                                       cparams$evidence_prepped[c_idx %in% indices_na, ],
                                       fill = TRUE),
                                 c_idx ~ variable, value.var = "val")[c_idx != 0, ]
    rows_na <- data.table::rbindlist(replicate(n_synth, rows_na, simplify = FALSE))
    if (nomatch == "force") {
      # nested recovery draw runs serial (a worker must not spawn its own backend)
      rows_na_sampled <- forge(params, n_synth = nrow(rows_na), sample_NAs = sample_NAs, parallel = FALSE)
      rows_na[is.na(rows_na)] <- rows_na_sampled[is.na(rows_na[, -1])]
    }
    x_synth[, c_idx := rep(indices_sampled, each = n_synth)]
    x_synth <- rbind(x_synth, rows_na, fill = TRUE)
    data.table::setorder(x_synth, c_idx)[, c_idx := NULL]
    x_synth <- post_x(x_synth, params, round)
  }
  x_synth
}
