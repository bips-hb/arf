#' @keywords internal
#' @noRd

# Internal helpers for the experimental mirai+mori backend.
# Activated by getOption("arf.backend", "foreach") == "mirai".
# Hidden from public API; mirai and mori in Suggests only.

arf_check_mirai_ready <- function() {
  if (!requireNamespace("mirai", quietly = TRUE)) {
    stop("arf.backend = 'mirai' requires the 'mirai' package. ",
         "Install it or set options(arf.backend = 'foreach').",
         call. = FALSE)
  }
  if (!requireNamespace("mori", quietly = TRUE)) {
    stop("arf.backend = 'mirai' requires the 'mori' package. ",
         "Install it or set options(arf.backend = 'foreach').",
         call. = FALSE)
  }
  status <- mirai::status()
  if (is.null(status$connections) || status$connections < 1L) {
    stop("arf.backend = 'mirai' requires daemons() to be set. ",
         "Call mirai::daemons(n) before forde().",
         call. = FALSE)
  }
  invisible(TRUE)
}

# Note: data.table objects come back from mori with truelength == 0
# (selfref dropped by the ALTREP round-trip), but read-only operations
# in the worker bodies tolerate this — see bench/compat/01-datatable.R.
# No setalloccol() repair is needed because workers only read the
# shared tables and build fresh data.tables from the results.

# Dispatch a per-tree worker over mirai daemons. Trees are split into
# one chunk per daemon; each task loops the per-tree worker over its
# block and rbinds locally (e.g. 100 trees -> 4 tasks on 4 daemons).
# At small scale this is timing-neutral vs one-task-per-tree, but it
# bounds the number of result objects serialized back to n_workers
# rather than num_trees, which matters for the large-grid cases.
arf_mirai_tree_map <- function(num_trees, worker_fn, shared_args) {
  st <- mirai::status()
  n_workers <- max(1L, as.integer(st$connections))
  n_chunks <- max(1L, min(n_workers, num_trees))
  chunks <- split(seq_len(num_trees),
                   rep(seq_len(n_chunks), length.out = num_trees))
  chunk_runner <- function(trees, worker_fn, shared_args) {
    parts <- lapply(trees, function(tr) {
      do.call(worker_fn, c(list(tr), shared_args))
    })
    data.table::rbindlist(parts)
  }
  res <- mirai::mirai_map(
    chunks,
    chunk_runner,
    .args = list(worker_fn = worker_fn, shared_args = shared_args)
  )[]
  data.table::rbindlist(res)
}
