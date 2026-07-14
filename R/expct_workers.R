#' @keywords internal
#' @noRd

# Per-step worker for expct(), extracted so foreach and mirai backends share one
# definition (params/evidence mori-shared under mirai). Calls arf internals
# (cforde, post_x, which.max.random) so mirai daemons must have arf loaded.
arf_expct_step <- function(step_, params, evidence, query, factor_cols,
                           evidence_row_mode, nomatch, verbose, round,
                           stepsize, stepsize_cforde, parallel_cforde) {
  # To avoid data.table check issues
  variable <- tree <- f_idx <- cvg <- wt <- V1 <- value <- val <- family <-
    mu <- sigma <- obs <- prob <- f_idx_uncond <- c_idx <- idx <- NA_share <-
    . <- I <- NULL

  # Prepare the event space
  if (is.null(evidence) || (ncol(evidence) == 2 && all(colnames(evidence) == c("f_idx", "wt")))) {
    cparams <- NULL
  } else {
    # Call cforde with part of the evidence for this step
    index_start <- (step_ - 1) * stepsize + 1
    index_end <- min(step_ * stepsize, nrow(evidence))
    evidence_part <- evidence[index_start:index_end, ]
    cparams <- cforde(params, evidence_part, evidence_row_mode, nomatch, verbose,
                      stepsize_cforde, parallel_cforde)
  }

  # omega contains the weight (wt) for each leaf (f_idx) for each condition (c_idx)
  if (is.null(cparams)) {
    if (is.null(evidence)) {
      num_trees <- params$forest[, max(tree)]
      omega <- params$forest[, .(f_idx, f_idx_uncond = f_idx, cvg)]
      omega[, `:=`(c_idx = 1, wt = cvg / num_trees)]
      omega[, cvg := NULL]
    } else {
      omega <- copy(evidence)
      omega[, f_idx_uncond := f_idx]
      omega[, c_idx := 1]
    }
  } else {
    omega <- cparams$forest[, .(c_idx, f_idx, f_idx_uncond, wt = cvg)]
  }
  omega <- omega[wt > 0, ]
  omega[, idx := .I]

  # Synthesize expectations for one block of conditions (subset of omega
  # rows; cparams/params merges below self-restrict via the c_idx join keys).
  synth_block <- function(omega_) {
    synth_cnt <- synth_cat <- NULL
    # Continuous data
    if (any(!factor_cols)) {
      if (is.null(cparams) || nrow(cparams$cnt) == 0) {
        psi_cond <- data.table()
      } else {
        psi_cond <- merge(omega_, cparams$cnt[variable %in% query, -c("cvg_factor", "f_idx_uncond")], by = c('c_idx', 'f_idx'),
                          sort = FALSE, allow.cartesian = TRUE)[prob > 0, ]
        # calculate absolute weights for sub-leaf areas (resulting from within-row or-conditions)
        if (any(psi_cond[, prob != 1])) {
          psi_cond[, wt := wt * prob]
          psi_cond[, I := seq_len(.N), by = .(variable, idx)]
        } else {
          psi_cond[, I := 1]
        }
        psi_cond[, prob := NULL]
      }
      psi <- unique(rbind(psi_cond,
                          merge(omega_, params$cnt[variable %in% query, ], by.x = 'f_idx_uncond', by.y = 'f_idx',
                                sort = FALSE, allow.cartesian = TRUE)[, `:=`(val = NA_real_, I = 1)]), by = c("c_idx", "f_idx", "variable", "I"))[, I := NULL]
      psi[NA_share == 1, wt := 0]
      cnt <- psi[is.na(val), val := sum(wt * mu) / sum(wt), by = .(c_idx, variable)]
      cnt <- unique(cnt[, .(c_idx, variable, val)])
      synth_cnt <- dcast(cnt, c_idx ~ variable, value.var = 'val')[, c_idx := NULL]
    }

    # Categorical data
    if (any(factor_cols)) {
      if (is.null(cparams) || nrow(cparams$cat) == 0) {
        psi <- merge(omega_, params$cat[variable %in% query, ], by.x = 'f_idx_uncond', by.y = 'f_idx', sort = FALSE, allow.cartesian = TRUE)
      } else {
        psi_cond <- merge(omega_, cparams$cat[variable %in% query, -c("cvg_factor", "f_idx_uncond")], by = c('c_idx', 'f_idx'),
                          sort = FALSE, allow.cartesian = TRUE)
        psi_uncond <- merge(omega_, params$cat[variable %in% query, ], by.x = 'f_idx_uncond', by.y = 'f_idx',
                            sort = FALSE, allow.cartesian = TRUE)
        psi_uncond_relevant <- psi_uncond[!psi_cond, on = .(idx, variable)]
        psi <- rbind(psi_cond, psi_uncond_relevant)
      }
      psi[NA_share == 1, wt := 0]
      cat <- psi[, sum(wt * prob), by = .(c_idx, variable, val)]
      cat <- setDT(cat)[, .SD[which.max.random(V1)], by = .(c_idx, variable)]
      synth_cat <- dcast(cat, c_idx ~ variable, value.var = 'val')[, c_idx := NULL]
    }
    cbind(synth_cnt, synth_cat)
  }

  # The merges in synth_block materialize (#matched leaves x #query variables)
  # rows PER CONDITION -- all at once for the step, three times over via
  # merge/rbind/unique, on every backend alike; at large forest x condition
  # counts this reaches tens of GB. Conditions are independent here
  # (every aggregation groups by c_idx, omega arrives sorted by c_idx, and the
  # per-group RNG order of which.max.random is ascending c_idx either way), so
  # when the estimated join size exceeds block_cap rows, process conditions in
  # blocks that keep each materialization bounded. A single condition cannot be
  # split; its leaves x variables product is the floor of this algorithm.
  # Lower via options(arf.block_rows) for tight-memory runs; see ?arf-options.
  block_cap <- max(1, as.numeric(getOption("arf.block_rows", 5e6)))
  n_vars <- max(1L, length(query))
  # as.double: nrow * n_vars overflows integer arithmetic exactly in the
  # large-forest regime the blocking exists for
  if (as.double(nrow(omega)) * n_vars <= block_cap || omega[, uniqueN(c_idx)] == 1L) {
    x_synth <- synth_block(omega)
  } else {
    sizes <- omega[, .N, by = c_idx]  # ascending c_idx
    g <- integer(nrow(sizes)); gi <- 1L; acc <- 0
    for (i in seq_len(nrow(sizes))) {
      r <- sizes$N[i] * n_vars
      if (acc > 0 && acc + r > block_cap) { gi <- gi + 1L; acc <- 0 }
      g[i] <- gi; acc <- acc + r
    }
    x_synth <- rbindlist(lapply(split(sizes$c_idx, g), function(cs) {
      synth_block(omega[c_idx %in% cs])
    }))
  }

  # Create dataset with expectations
  x_synth <- post_x(x_synth, params, round)

  if (evidence_row_mode == "separate" & any(omega[, is.na(f_idx)])) {
    setDT(x_synth)
    indices_na <- cparams$forest[is.na(f_idx), c_idx]
    indices_sampled <- cparams$forest[!is.na(f_idx), unique(c_idx)]
    rows_na <- dcast(rbind(data.table(c_idx = 0, variable = params$meta[, variable]),
                           cparams$evidence_prepped[c_idx %in% indices_na, ],
                           fill = TRUE),
                     c_idx ~ variable, value.var = "val")[c_idx != 0, ]
    if (nomatch == "force") {
      # nested recovery runs serial (a worker must not spawn its own backend)
      rows_na_sampled <- expct(params, parallel = FALSE)
      rows_na[is.na(rows_na)] <- rows_na_sampled[is.na(rows_na[, -1])]
    }
    x_synth[, c_idx := indices_sampled]
    x_synth <- rbind(x_synth, rows_na, fill = TRUE)
    setorder(x_synth, c_idx)[, c_idx := NULL]
    x_synth <- post_x(x_synth, params, round)
  }

  x_synth
}
