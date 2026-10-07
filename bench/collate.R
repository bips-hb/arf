# Deltas, the PR report, and the anchor history. Thresholds come from
# bench/DESIGN.md: memory above 3% is real, time below 10% is inconclusive.

BENCH_MEM_THRESHOLD <- 3
# A marginal below this fraction of its floor is below the instrument's
# resolution; see the note in bench_deltas().
BENCH_MEM_RESOLUTION <- 1
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
  # Resolution guard. The marginal is a difference against a floor of an R
  # interpreter plus arf plus the fixture, so when it is small relative to that
  # floor the shared-cgroup delta cannot resolve it and no replicate count
  # helps. Measured same-commit spread against the marginal-to-floor ratio:
  # 0.29 -> 14-23%, 0.73 -> 4.2%, 2.0 -> <=3.6%, so the op must allocate at
  # least as much as the fixed floor before a 3% verdict means anything.
  # Suppress it rather than report noise as a finding.
  # Per CELL, not per row: a delta has two sides, so if either the ref or its
  # baseline is below the instrument's resolution then so is the comparison.
  # Judging each row on its own ratio let one ref in a cell read "real" while
  # the other read "too small", which is how a +10.3% noise reading got a
  # verdict on the cluster's first run.
  thin_row <- !is.na(rows$peak_delta_mb) &
    !is.na(rows$floor_mb) &
    rows$peak_delta_mb < BENCH_MEM_RESOLUTION * rows$floor_mb
  thin_cell <- rows$key %in% unique(rows$key[thin_row])
  rows$mem_verdict[thin_cell & rows$mem_verdict %in% c("real", "inconclusive")] <-
    "cell too small to resolve memory"
  rows$key <- NULL
  rows
}


.bench_fmt_pct <- function(x) ifelse(is.na(x), "   n/a", sprintf("%+6.1f%%", x))

bench_report <- function(rows, baseline = "main") {
  d <- bench_deltas(rows, baseline)
  base <- d[d$ref == baseline, ]
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
  idx <- match(.bench_cell_key(d), .bench_cell_key(base))

  num <- function(x, digits = 1) ifelse(is.na(x), "n/a", formatC(x, format = "f", digits = digits))
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
      paste0(
        "Memory is the marginal (peak minus the measured per-cell floor), median of ",
        "%s paired rounds. Deltas at or above %g%% are real; time deltas below %g%% ",
        "are inconclusive. A marginal below its floor is below the instrument's ",
        "resolution and gets no verdict."
      ),
      paste(unique(stats::na.omit(rows$mem_reps)), collapse = "/"),
      BENCH_MEM_THRESHOLD,
      BENCH_TIME_THRESHOLD
    ),
    "",
    paste(
      "| op | backend | w | n | trees | ref | marginal MB | floor MB |",
      "mem | time s | time |"
    ),
    "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"
  )
  for (i in seq_len(nrow(d))) {
    r <- d[i, ]
    b <- base[idx[i], ]
    out <- c(
      out,
      sprintf(
        "| %s | %s | %s | %g | %g | `%s` | %s vs %s | %s | %s %s | %s vs %s | %s %s |",
        r$op,
        r$backend,
        if (is.na(r$workers)) "-" else r$workers,
        r$n,
        r$trees,
        r$ref,
        num(r$peak_delta_mb),
        num(b$peak_delta_mb),
        num(r$floor_mb),
        .bench_fmt_pct(r$delta_mem_pct),
        r$mem_verdict,
        num(r$time_median, 3),
        num(b$time_median, 3),
        .bench_fmt_pct(r$delta_time_pct),
        r$time_verdict
      )
    )
    if (isTRUE(r$digest_differs)) {
      out <- c(
        out,
        sprintf(
          paste0(
            "| | | | | | | | | **OUTPUT DIFFERS** (%s digest) against %s: ",
            "this delta is not apples-to-apples, confirm the change is intended | | |"
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

# Shared tail of a run and a collate: one place that writes the CSV, appends
# anchor rows to the history, and renders the report.
bench_write_results <- function(rows) {
  dir.create("bench/results", showWarnings = FALSE, recursive = TRUE)
  out <- file.path(
    "bench/results",
    sprintf("bench-%s-%s.csv", rows$tier[1], format(Sys.time(), "%Y%m%d-%H%M%S"))
  )
  utils::write.csv(rows, out, row.names = FALSE)
  message("Written to ", out)
  # Anchor rows are appended automatically; only the commit is manual.
  bench_history_append(rows)
  report <- bench_report(rows)
  cat(report, sep = "\n")
  writeLines(report, "bench/results/report.md")
  message("Report at bench/results/report.md")
  invisible(rows)
}
