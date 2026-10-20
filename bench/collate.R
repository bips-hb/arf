# Deltas, verdicts, the report, and the two committed CSVs.

BENCH_MEM_THRESHOLD <- 3
# A marginal below this fraction of its floor is below the instrument's
# resolution: same-commit spread at marginal/floor 0.29 was 14-23%, at 2.0 under 4%.
BENCH_MEM_RESOLUTION <- 1
# Above this many comparison rows the report summarises instead of tabulating.
BENCH_REPORT_ROWS <- 40
# Below this many seconds a percentage is jitter, not a finding.
BENCH_TIME_FLOOR_S <- 0.5
# A peak this close to the cgroup limit is the limit, not the workload.
BENCH_MEM_CAP_FRACTION <- 0.95
BENCH_TIME_THRESHOLD <- 10

# Evidence and synth counts are part of the cell: forge and expct carry
# variants at otherwise identical settings.
.bench_cell_key <- function(rows) {
  paste(
    rows$op,
    rows$backend,
    rows$workers,
    rows$n,
    rows$p,
    rows$trees,
    rows$rowmode,
    rows$n_evidence,
    rows$n_synth,
    sep = "|"
  )
}

# A zero or NA on either side is a failed measurement, not a 100% change.
.bench_pct <- function(new, old) {
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
  # On the marginal: the floor would dilute a raw-peak percentage several-fold.
  rows$delta_mem_pct <- .bench_pct(rows$peak_delta_mb, base$peak_delta_mb[idx])
  rows$delta_time_pct <- .bench_pct(rows$time_median, base$time_median[idx])
  rows$mem_verdict <- .bench_verdict(rows$delta_mem_pct, BENCH_MEM_THRESHOLD)
  rows$time_verdict <- .bench_verdict(rows$delta_time_pct, BENCH_TIME_THRESHOLD)
  # Results predating the column are assumed stable, so old runs keep their verdicts.
  stable <- function(x) if (is.null(x)) TRUE else is.na(x) | x
  # An incomplete result is not a changed result; checked first.
  rows$out_rows_differ <- if (is.null(rows$out_rows)) {
    rep(FALSE, nrow(rows))
  } else {
    !is.na(rows$out_rows) &
      !is.na(base$out_rows[idx]) &
      rows$out_rows != base$out_rows[idx]
  }
  rows$digest_differs <- !rows$out_rows_differ &
    !is.na(rows$digest) &
    !is.na(base$digest[idx]) &
    rows$digest != base$digest[idx] &
    stable(rows$digest_stable) &
    stable(base$digest_stable[idx])
  too_fast <- (!is.na(rows$time_median) & rows$time_median < BENCH_TIME_FLOOR_S) |
    (!is.na(base$time_median[idx]) & base$time_median[idx] < BENCH_TIME_FLOOR_S)
  rows$time_verdict[too_fast & rows$time_verdict %in% c("real", "inconclusive")] <-
    "too fast to time reliably"
  # Timings per ref are iters x mem_reps; one of either cannot carry a percentage.
  n_timings <- ifelse(is.na(rows$iters) | is.na(rows$mem_reps), NA_integer_, rows$iters * rows$mem_reps)
  one_timing <- !is.na(n_timings) & n_timings < 2L & rows$time_verdict %in% c("real", "inconclusive")
  rows$time_verdict[one_timing] <- "single sample, not a finding"
  single <- !is.na(rows$mem_reps) & rows$mem_reps < 2L & rows$mem_verdict == "real"
  rows$mem_verdict[single] <- "single sample, not a finding"
  # Every guard below is per comparison: a delta has two sides, so it is
  # invalid when the ref or its baseline is affected, and unaffected when some
  # third ref in the cell is.
  capped <- !is.na(rows$mem_limit_mb) &
    !is.na(rows$peak_mb) &
    rows$peak_mb >= BENCH_MEM_CAP_FRACTION * rows$mem_limit_mb
  capped_pair <- capped | (!is.na(idx) & capped[match(rows$key, base$key)])
  thin_row <- !is.na(rows$peak_delta_mb) &
    !is.na(rows$floor_mb) &
    rows$peak_delta_mb < BENCH_MEM_RESOLUTION * rows$floor_mb
  thin_pair <- thin_row | (!is.na(idx) & thin_row[match(rows$key, base$key)])
  rows$mem_verdict[thin_pair & rows$mem_verdict %in% c("real", "inconclusive")] <-
    "cell too small to resolve memory"
  rows$mem_verdict[capped_pair] <- "hit the memory cap, not a measurement"
  # memory.high sits at the cap, so a capped run was throttled: its time is not the op's.
  rows$time_verdict[capped_pair] <- "hit the memory cap, not a measurement"
  # Last, so it survives the guards above: a ref that cannot run a cell is a
  # finding, not a missing baseline.
  rows$mem_verdict[is.na(rows$peak_delta_mb) & !is.na(base$peak_delta_mb[idx])] <-
    "ref did not complete"
  rows$time_verdict[is.na(rows$time_median) & !is.na(base$time_median[idx])] <-
    "ref did not complete"
  neither <- is.na(rows$peak_delta_mb) & is.na(base$peak_delta_mb[idx])
  rows$mem_verdict[neither] <- "no ref completed this cell"
  rows$time_verdict[is.na(rows$time_median) & is.na(base$time_median[idx])] <-
    "no ref completed this cell"
  rows$key <- NULL
  rows
}


.bench_fmt_pct <- function(x) ifelse(is.na(x), "   n/a", sprintf("%+6.1f%%", x))

# forge/expct/lik cells differ only in these, so without them the table shows
# rows that look like duplicates.
.bench_fmt_variant <- function(r) {
  parts <- c(
    if (!is.na(r$n_evidence)) paste0("e", r$n_evidence),
    if (!is.na(r$rowmode)) r$rowmode,
    if (!is.na(r$n_synth)) paste0("s", r$n_synth),
    if (!is.na(r$n_folds)) paste0("k", r$n_folds)
  )
  if (length(parts)) paste(parts, collapse = " ") else "-"
}

.bench_report_header <- c(
  paste(
    "| op | backend | w | n | trees | variant | ref | marginal MB | floor MB |",
    "mem | time s | time |"
  ),
  "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"
)

bench_report <- function(rows, baseline = "main") {
  d <- bench_deltas(rows, baseline)
  base <- d[d$ref == baseline, ]
  d <- d[d$ref != baseline, ]
  metrics <- unique(rows$metric)
  if (length(metrics) > 1L) {
    warning(
      "results mix memory metrics (",
      paste(metrics, collapse = ", "),
      "), which are not comparable; split the run by metric",
      call. = FALSE
    )
  }
  idx <- match(.bench_cell_key(d), .bench_cell_key(base))
  full_table <- nrow(d) <= BENCH_REPORT_ROWS

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
    if (all(is.na(rows$mem_limit_mb))) {
      paste(
        "**No memory limit recorded in these results**, so the cap guard could",
        "not be applied: a peak that was clamped by its cgroup limit would",
        "still read as a measurement here. Results from before `mem_limit_mb`",
        "entered the schema."
      )
    } else {
      NULL
    },
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
    ""
  )
  shown <- if (full_table) {
    seq_len(nrow(d))
  } else {
    which(
      d$mem_verdict == "real" |
        d$time_verdict == "real" |
        d$digest_differs |
        d$out_rows_differ
    )
  }
  if (!full_table) {
    # A full-tier run has hundreds of rows: summarise by op, show the findings.
    counts <- table(d$op, d$mem_verdict)
    out <- c(
      out,
      sprintf(
        "%d comparisons across %d cells. Memory verdicts by op:",
        nrow(d),
        length(unique(.bench_cell_key(d)))
      ),
      "",
      "```",
      utils::capture.output(print(counts)),
      "```",
      "",
      sprintf(
        "Rows below are the %d comparison(s) with a verdict of `real`, or a digest mismatch. Full data in the CSV.",
        length(shown)
      ),
      ""
    )
  }
  out <- c(out, .bench_report_header)
  for (i in shown) {
    r <- d[i, ]
    b <- base[idx[i], ]
    out <- c(
      out,
      sprintf(
        "| %s | %s | %s | %g | %g | %s | `%s` | %s vs %s | %s | %s %s | %s vs %s | %s %s |",
        r$op,
        r$backend,
        if (is.na(r$workers)) "-" else r$workers,
        r$n,
        r$trees,
        .bench_fmt_variant(r),
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
    if (isTRUE(r$out_rows_differ)) {
      out <- c(
        out,
        sprintf(
          paste0(
            "| | | | | | | | | | **INCOMPLETE RESULT**: %s rows against %s's %s. ",
            "A parallel worker that dies drops its whole chunk and arf returns ",
            "the short result without erroring, so neither side's timing or ",
            "peak is for the full workload | | |"
          ),
          format(r$out_rows),
          baseline,
          format(b$out_rows)
        )
      )
    } else if (isTRUE(r$digest_differs)) {
      out <- c(
        out,
        sprintf(
          paste0(
            "| | | | | | | | | | **OUTPUT DIFFERS** (%s digest) against %s: ",
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

# Anchors only: HEAD and main move. Committed by hand, for runs worth keeping.
bench_history_append <- function(rows, path = "bench/history.csv") {
  anchors <- rows[!rows$ref %in% c("HEAD", "main"), bench_schema()]
  keep <- anchors[!is.na(anchors$peak_delta_mb) | !is.na(anchors$time_median), ]
  if (!nrow(keep)) {
    message(
      "Nothing appended to ",
      path,
      ": of ",
      nrow(rows),
      " row(s), ",
      nrow(anchors),
      " were anchors and ",
      nrow(anchors) - nrow(keep),
      " of those carried no measurement (HEAD and main are never history)."
    )
    return(invisible(keep))
  }
  if (nrow(keep) < nrow(anchors)) {
    message(
      nrow(anchors) - nrow(keep),
      " anchor row(s) carried no measurement and ",
      "are not in the history; see the report for which cells those were."
    )
  }
  added <- nrow(keep)
  keep <- .bench_bind_old(keep, path)
  # Newest row per configuration, or the trend panel draws a zig-zag.
  keep <- .bench_newest_by(
    keep,
    key = c(
      "ref",
      "arf_version",
      "op",
      "backend",
      "workers",
      "n",
      "p",
      "trees",
      "rowmode",
      "n_evidence",
      "n_synth",
      "host",
      "metric"
    ),
    stamp = "timestamp",
    sort_by = c("op", "n", "trees", "backend", "workers")
  )
  utils::write.csv(keep, path, row.names = FALSE)
  message("History at ", path, ": ", nrow(keep), " row(s) total, ", added, " new or replaced.")
  invisible(keep)
}

# The committed CSVs outlive the schema that wrote them: fill missing columns
# rather than refuse the file.
.bench_bind_old <- function(new, path) {
  if (!file.exists(path) || length(readLines(path, warn = FALSE)) < 2L) {
    return(new)
  }
  old <- utils::read.csv(path, stringsAsFactors = FALSE)
  for (nm in setdiff(names(new), names(old))) {
    old[[nm]] <- NA
  }
  rbind(old[, names(new)], new)
}

# Newest row per key, replacing rather than accumulating.
.bench_newest_by <- function(df, key, stamp, sort_by = key) {
  dt <- data.table::as.data.table(df)
  data.table::setorderv(dt, c(key, stamp))
  dt <- unique(dt, by = key, fromLast = TRUE)
  data.table::setorderv(dt, sort_by)
  as.data.frame(dt)
}

# Per-call seconds per cell, committed so slurm walltimes come from data.
# Separate from the history: scheduling needs the slowest ref in a cell,
# whichever that is, and cells no anchor runs.
bench_timings_append <- function(rows, path = "bench/timings.csv") {
  have <- rows[!is.na(rows$time_median), ]
  if (!nrow(have)) {
    return(invisible(NULL))
  }
  key <- c("op", "backend", "workers", "n", "p", "trees", "n_evidence", "n_synth", "rowmode")
  # data.table, not aggregate(): aggregate() drops groups with an NA by-column.
  agg <- as.data.frame(
    data.table::as.data.table(have)[,
      list(call_seconds = max(time_median, na.rm = TRUE)),
      by = key
    ]
  )
  agg$host <- rows$host[1]
  agg$measured <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  agg <- .bench_bind_old(agg, path)
  agg <- .bench_newest_by(
    agg,
    key = key,
    stamp = "measured",
    sort_by = c("op", "n", "n_evidence", "rowmode", "backend", "workers")
  )
  utils::write.csv(agg, path, row.names = FALSE)
  message("Timings at ", path, ": ", nrow(agg), " cell(s) known to the scheduler.")
  invisible(agg)
}

# Shared tail of a run and a collate.
bench_write_results <- function(rows) {
  dir.create("bench/results", showWarnings = FALSE, recursive = TRUE)
  out <- file.path(
    "bench/results",
    sprintf("bench-%s-%s.csv", rows$tier[1], format(Sys.time(), "%Y%m%d-%H%M%S"))
  )
  utils::write.csv(rows, out, row.names = FALSE)
  message("Written to ", out)
  bench_history_append(rows)
  bench_timings_append(rows)
  report <- bench_report(rows)
  cat(report, sep = "\n")
  writeLines(report, "bench/results/report.md")
  message("Report at bench/results/report.md")
  invisible(rows)
}
