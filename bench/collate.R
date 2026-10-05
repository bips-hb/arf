# Deltas, the PR report, and the anchor history. Thresholds come from
# bench/DESIGN.md: memory above 3% is real, time below 10% is inconclusive.

BENCH_MEM_THRESHOLD <- 3
BENCH_TIME_THRESHOLD <- 10

.bench_cell_key <- function(rows) {
  paste(rows$op, rows$backend, rows$workers, rows$n, rows$p, rows$trees, rows$rowmode, sep = "|")
}

.bench_pct <- function(new, old) {
  # A baseline of 0 or NA means the baseline measurement failed, which is not a
  # 100% improvement.
  # Both sides: a zero on the NEW side is a sampler that caught nothing, which
  # would otherwise render as -100% with verdict "real".
  ifelse(
    is.na(new) | is.na(old) | old <= 0 | new <= 0,
    NA_real_,
    (new - old) / old * 100
  )
}

.bench_verdict <- function(pct, threshold) {
  ifelse(
    is.na(pct),
    "no baseline",
    ifelse(
      pct == 0,
      "unchanged",
      ifelse(abs(pct) >= threshold, "real", "inconclusive")
    )
  )
}

bench_deltas <- function(rows, baseline = "main") {
  rows$key <- .bench_cell_key(rows)
  base <- rows[rows$ref == baseline, ]
  if (!nrow(base)) {
    stop("baseline ref '", baseline, "' is not in the results", call. = FALSE)
  }
  idx <- match(rows$key, base$key)
  # On the marginal, not the floor-inflated total (bench/DESIGN.md).
  rows$delta_mem_pct <- .bench_pct(rows$peak_delta_mb, base$peak_delta_mb[idx])
  rows$delta_time_pct <- .bench_pct(rows$time_median, base$time_median[idx])
  rows$mem_verdict <- .bench_verdict(rows$delta_mem_pct, BENCH_MEM_THRESHOLD)
  rows$time_verdict <- .bench_verdict(rows$delta_time_pct, BENCH_TIME_THRESHOLD)
  rows$digest_differs <- !is.na(rows$digest) & !is.na(base$digest[idx]) & rows$digest != base$digest[idx]
  # One peak sample cannot support a percentage claim: peak memory is a
  # stochastic function of GC scheduling, and identical code re-run differs by
  # tens of percent on the marginal (measured same-commit: -2.7% to -11.6% at
  # n=1e3, up to -38.6% at n=1e4). Refuse to call such a delta a finding.
  single <- !is.na(rows$mem_reps) & rows$mem_reps < 2L & rows$mem_verdict == "real"
  rows$mem_verdict[single] <- "single sample, not a finding"
  rows$key <- NULL
  rows
}


.bench_fmt_pct <- function(x) ifelse(is.na(x), "   n/a", sprintf("%+6.1f%%", x))

bench_report <- function(rows, baseline = "main") {
  d <- bench_deltas(rows, baseline)
  d <- d[d$ref != baseline, ]
  metrics <- unique(rows$metric)
  if (length(metrics) > 1L) {
    # DESIGN: rows with different metrics must never be compared.
    warning(
      "results mix memory metrics (",
      paste(metrics, collapse = ", "),
      "), which are not comparable; split the run by metric",
      call. = FALSE
    )
  }
  out <- c(
    "# arf benchmark report",
    "",
    sprintf(
      "Baseline: `%s`. Metric: %s. Host(s): %s, kernel %s, R %s.",
      baseline,
      paste(metrics, collapse = " / "),
      paste(unique(rows$host), collapse = ", "),
      paste(unique(rows$kernel), collapse = ", "),
      paste(unique(rows$r_version), collapse = ", ")
    ),
    sprintf(
      "Memory deltas at or above %g%% are real; time deltas below %g%% are inconclusive.",
      BENCH_MEM_THRESHOLD,
      BENCH_TIME_THRESHOLD
    ),
    ""
  )
  for (i in seq_len(nrow(d))) {
    r <- d[i, ]
    out <- c(
      out,
      sprintf(
        "- %-14s %-10s w=%-3s n=%-6g trees=%-4g %-12s mem %s (%s)  time %s (%s)",
        r$op,
        r$backend,
        if (is.na(r$workers)) "-" else as.character(r$workers),
        r$n,
        r$trees,
        r$ref,
        .bench_fmt_pct(r$delta_mem_pct),
        r$mem_verdict,
        .bench_fmt_pct(r$delta_time_pct),
        r$time_verdict
      )
    )
    if (isTRUE(r$digest_differs)) {
      out <- c(
        out,
        sprintf(
          paste0(
            "    !! OUTPUT DIFFERS (%s digest) against %s: this delta is ",
            "not apples-to-apples, confirm the change is intended"
          ),
          r$digest_kind,
          baseline
        )
      )
    }
  }
  c(out, "")
}

# Anchors only: HEAD and main move, so recording them would be noise. Committed
# by hand afterwards, for the runs worth keeping.
bench_history_append <- function(rows, path = "bench/history.csv") {
  keep <- rows[!rows$ref %in% c("HEAD", "main"), bench_schema()]
  if (!nrow(keep)) {
    return(invisible(keep))
  }
  if (file.exists(path) && length(readLines(path, warn = FALSE)) > 1L) {
    old <- utils::read.csv(path, stringsAsFactors = FALSE)
    keep <- rbind(old, keep)
  }
  utils::write.csv(keep, path, row.names = FALSE)
  invisible(keep)
}
