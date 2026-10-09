# Deltas, the PR report, and the anchor history. Thresholds come from
# Memory above 3% is real, time below 10% is inconclusive.

BENCH_MEM_THRESHOLD <- 3
# A marginal below this fraction of its floor is below the instrument's
# resolution; see the note in bench_deltas().
BENCH_MEM_RESOLUTION <- 1
# Above this many comparison rows the report summarises instead of tabulating.
BENCH_REPORT_ROWS <- 40
# Below this many seconds a percentage is not a finding: process and GC jitter
# at tens of milliseconds easily exceeds the 10% line. Measured on a real quick
# tier run, the ONLY "real" time verdict under a second was the fastest cell in
# the grid (62 ms, +12.7%), while every cell above 1 s read inconclusive.
BENCH_TIME_FLOOR_S <- 0.5
# A peak this close to the cgroup limit is the limit, not the workload: six
# cells of a real cluster run sat at 99.75% of a 256 GB cap, identical to the
# megabyte across a fivefold difference in n.
BENCH_MEM_CAP_FRACTION <- 0.95
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
  # On the marginal, not the floor-inflated total: a ~200 MB floor of
  # interpreter plus arf plus fixture dilutes a raw-peak percentage fivefold.
  rows$delta_mem_pct <- .bench_pct(rows$peak_delta_mb, base$peak_delta_mb[idx])
  rows$delta_time_pct <- .bench_pct(rows$time_median, base$time_median[idx])
  rows$mem_verdict <- .bench_verdict(rows$delta_mem_pct, BENCH_MEM_THRESHOLD)
  rows$time_verdict <- .bench_verdict(rows$delta_time_pct, BENCH_TIME_THRESHOLD)
  rows$digest_differs <- !is.na(rows$digest) & !is.na(base$digest[idx]) & rows$digest != base$digest[idx]
  # One peak sample cannot support a percentage claim: peak memory is a
  # stochastic function of GC scheduling, and identical code re-run differs by
  # tens of percent on the marginal (measured same-commit: -2.7% to -11.6% at
  # n=1e3, up to -38.6% at n=1e4). Refuse to call such a delta a finding.
  too_fast <- (!is.na(rows$time_median) & rows$time_median < BENCH_TIME_FLOOR_S) |
    (!is.na(base$time_median[idx]) & base$time_median[idx] < BENCH_TIME_FLOOR_S)
  rows$time_verdict[too_fast & rows$time_verdict %in% c("real", "inconclusive")] <-
    "too fast to time reliably"
  # Time has the same replication problem as memory and had no guard for it: a
  # validation run at --iters 1 --mem-reps 1 produced one timing per ref and
  # still printed "real" on a 13% delta. Timings per ref are iters x mem_reps.
  n_timings <- ifelse(is.na(rows$iters) | is.na(rows$mem_reps), NA_integer_, rows$iters * rows$mem_reps)
  one_timing <- !is.na(n_timings) & n_timings < 2L & rows$time_verdict %in% c("real", "inconclusive")
  rows$time_verdict[one_timing] <- "single sample, not a finding"
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
  capped <- !is.na(rows$mem_limit_mb) &
    !is.na(rows$peak_mb) &
    rows$peak_mb >= BENCH_MEM_CAP_FRACTION * rows$mem_limit_mb
  # Per COMPARISON, not per cell: a delta has two sides, so it is invalid
  # when the ref or its baseline capped, and unaffected when some third ref
  # in the same cell did. Where only the hungriest anchor caps, HEAD against
  # main is still a measurement and must not be discarded with it.
  capped_pair <- capped | (!is.na(idx) & capped[match(rows$key, base$key)])
  rows$mem_verdict[capped_pair] <- "hit the memory cap, not a measurement"
  thin_row <- !is.na(rows$peak_delta_mb) &
    !is.na(rows$floor_mb) &
    rows$peak_delta_mb < BENCH_MEM_RESOLUTION * rows$floor_mb
  # Same pairing. Usually identical to a cell-wide rule anyway, since both
  # sides share the cell's floor.
  thin_pair <- thin_row | (!is.na(idx) & thin_row[match(rows$key, base$key)])
  rows$mem_verdict[thin_pair & rows$mem_verdict %in% c("real", "inconclusive")] <-
    "cell too small to resolve memory"
  rows$mem_verdict[capped_pair] <- "hit the memory cap, not a measurement"
  # Last, so it survives the guards above. "no baseline" blames the wrong side
  # when the baseline measured fine and the ref is the one that died: the
  # evidence=1000 expct cell needs ~1.8 TB on cran-0.2.5, more than a node has,
  # so its child is OOM-killed and returns NA while main returns 267 GB. That
  # a ref cannot run a cell at all is a finding, not a missing baseline.
  rows$mem_verdict[is.na(rows$peak_delta_mb) & !is.na(base$peak_delta_mb[idx])] <-
    "ref did not complete"
  rows$time_verdict[is.na(rows$time_median) & !is.na(base$time_median[idx])] <-
    "ref did not complete"
  # Both sides dead means the cell is beyond this hardware for every ref, not
  # missing a reference point.
  neither <- is.na(rows$peak_delta_mb) & is.na(base$peak_delta_mb[idx])
  rows$mem_verdict[neither] <- "no ref completed this cell"
  rows$time_verdict[is.na(rows$time_median) & is.na(base$time_median[idx])] <-
    "no ref completed this cell"
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
  # A full-tier run is 704 non-baseline rows, which is not something anyone
  # reads in a pull request. Past a readable size, summarise by op and show
  # only the rows that carry a verdict, with the CSV as the full record.
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
    "",
    paste(
      "| op | backend | w | n | trees | ref | marginal MB | floor MB |",
      "mem | time s | time |"
    ),
    "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"
  )
  shown <- if (full_table) {
    seq_len(nrow(d))
  } else {
    which(d$mem_verdict == "real" | d$time_verdict == "real" | d$digest_differs)
  }
  if (!full_table) {
    counts <- table(d$op, d$mem_verdict)
    out <- c(
      out[seq_len(length(out) - 2L)],
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
      "",
      paste(
        "| op | backend | w | n | trees | ref | marginal MB | floor MB |",
        "mem | time s | time |"
      ),
      "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"
    )
  }
  for (i in shown) {
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
  anchors <- rows[!rows$ref %in% c("HEAD", "main"), bench_schema()]
  # A row with no measurement is not history. Say so: a cell whose anchor was
  # OOM-killed leaves nothing to append, and a silent return looks exactly like
  # a collate that never ran.
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
  if (file.exists(path) && length(readLines(path, warn = FALSE)) > 1L) {
    old <- utils::read.csv(path, stringsAsFactors = FALSE)
    keep <- rbind(old, keep)
  }
  # One row per anchor per cell per host and metric, newest kept: this file is
  # committed and plotted, so re-running a configuration must replace its row
  # rather than accumulate duplicates the trend panel would draw as a zig-zag.
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

# Newest row per key, replacing rather than accumulating. Both committed CSVs
# need it: re-running a configuration must replace its row, or the trend panel
# draws a zig-zag and the scheduler averages stale timings.
.bench_newest_by <- function(df, key, stamp, sort_by = key) {
  dt <- data.table::as.data.table(df)
  data.table::setorderv(dt, c(key, stamp))
  dt <- unique(dt, by = key, fromLast = TRUE)
  data.table::setorderv(dt, sort_by)
  as.data.frame(dt)
}

# Measured per-call seconds per cell, committed so the scheduler can size
# slurm walltimes from data instead of one global figure. Separate from the
# history, which keeps anchors only: scheduling cares about the SLOWEST ref in
# a cell, whichever that is, and about cells no anchor ever runs.
bench_timings_append <- function(rows, path = "bench/timings.csv") {
  have <- rows[!is.na(rows$time_median), ]
  if (!nrow(have)) {
    return(invisible(NULL))
  }
  key <- c("op", "backend", "workers", "n", "p", "trees", "n_evidence", "rowmode")
  # Worst ref in the cell, not the mean: a walltime that fits the fast ref and
  # not the slow one loses the whole cell.
  # data.table, not aggregate(): aggregate() drops every group with an NA in a
  # `by` column, which here is every sequential cell and every op with no
  # evidence arguments.
  agg <- as.data.frame(
    data.table::as.data.table(have)[,
      list(call_seconds = max(time_median, na.rm = TRUE)),
      by = key
    ]
  )
  agg$host <- rows$host[1]
  agg$measured <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  if (file.exists(path) && length(readLines(path, warn = FALSE)) > 1L) {
    old <- utils::read.csv(path, stringsAsFactors = FALSE)
    agg <- rbind(old[, names(agg)], agg)
  }
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
  bench_timings_append(rows)
  report <- bench_report(rows)
  cat(report, sep = "\n")
  writeLines(report, "bench/results/report.md")
  message("Report at bench/results/report.md")
  invisible(rows)
}
