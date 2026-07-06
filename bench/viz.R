#!/usr/bin/env Rscript
# Plot the latest benchmark results: runtime and peak memory across backends
# and worker counts, faceted by benchmark setting (op, rowmode, evidence/synth
# scale) and problem size. Writes a multi-page PDF to bench/results/.
#
# Usage:
#   Rscript bench/viz.R              # newest sweep-ops CSV per config label
#   Rscript bench/viz.R file1.csv... # explicit files (old forde-only sweep-*.csv work too)
suppressMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args)) {
  files <- args
} else {
  files <- list.files("bench/results", pattern = "^sweep-ops-.*\\.csv$",
                      full.names = TRUE)
  if (!length(files)) stop("no sweep-ops CSVs in bench/results/")
  # newest file per config label (filename: sweep-ops[-<label>]-<stamp>.csv)
  info <- data.table(file = files)
  info[, base := sub("\\.csv$", "", sub("^sweep-ops-?", "", basename(file)))]
  info[, stamp := sub("^.*?([0-9]{8}-[0-9]{6})$", "\\1", base)]
  info[, label := sub("-?[0-9]{8}-[0-9]{6}$", "", base)]
  info <- info[order(-stamp), .SD[1], by = label]
  files <- info$file
}
message("Reading: ", paste(basename(files), collapse = ", "))
res <- rbindlist(lapply(files, fread), fill = TRUE)
if (!"op" %in% names(res)) res[, op := NA_character_]
res[is.na(op), op := "forde"]  # old forde-only sweep.R CSVs

# One facet label per benchmark setting; forge/expct differ by rowmode/scale.
res[, setting := op]
if ("rowmode" %in% names(res)) {
  res[op %in% c("forge", "expct") & !is.na(rowmode),
      setting := paste0(op, " ", rowmode, " (evi=", n_evidence,
                        fifelse(!is.na(n_synth),
                                paste0(", synth=", n_synth), ""), ")")]
}
res[, minutes := seconds / 60]
res[, gib := peak_mb / 1024]

# sequential has one measurement (workers = NA) per facet: draw as a per-facet
# reference line, not a point in the workers series
seq_res <- res[backend == "sequential"]
par_res <- res[backend != "sequential"]

plot_metric <- function(yvar, ylab, title, log10 = FALSE) {
  p <- ggplot(par_res, aes(x = workers, y = .data[[yvar]],
                           color = backend, fill = backend)) +
    facet_wrap(vars(setting, n, trees), labeller = label_both,
               ncol = 4, scales = "free_y") +
    geom_hline(data = seq_res, aes(yintercept = .data[[yvar]]),
               color = "blue", linetype = "dashed") +
    geom_line() +
    geom_point(size = 2.5, key_glyph = "rect") +
    scale_x_continuous(breaks = sort(unique(par_res$workers))) +
    labs(title = title, subtitle = "Dashed: sequential run",
         x = "# of workers", y = ylab, color = NULL, fill = NULL) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", plot.title.position = "plot")
  if (log10) p <- p + scale_y_log10()
  p
}

out <- sprintf("bench/results/viz-%s.pdf", format(Sys.time(), "%Y%m%d-%H%M%S"))
n_facets <- nrow(unique(par_res[, .(setting, n, trees)]))
pdf(out, width = 16, height = max(9, 3 * ceiling(n_facets / 4)))
print(plot_metric("minutes", "Time (minutes)", "Runtime in total"))
print(plot_metric("minutes", "Time (minutes) log10", "Runtime in total", log10 = TRUE))
print(plot_metric("gib", "Memory (GiB)", "Peak (shared) memory consumption"))
print(plot_metric("gib", "Memory (GiB) log10", "Peak (shared) memory consumption", log10 = TRUE))
invisible(dev.off())
cat("Written to", out, "\n")
