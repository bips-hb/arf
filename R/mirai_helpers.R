#' @keywords internal
#' @noRd

# Internal helpers for the experimental mirai+mori backend.
# Activated by getOption("arf.backend", "foreach") == "mirai".
# Hidden from public API; mirai and mori in Suggests only.

# Session-scoped state for once-per-run notifications.
.arf_env <- new.env(parent = emptyenv())

# Emit a backend notification at most once per distinct state per session,
# and only when getOption("arf.verbose", TRUE) is TRUE (opt-out switch).
arf_backend_inform <- function(msg, key) {
  if (!isTRUE(getOption("arf.verbose", TRUE))) {
    return(invisible(FALSE))
  }
  shown <- get0("backend_shown", envir = .arf_env, ifnotfound = character(0))
  if (key %in% shown) {
    return(invisible(FALSE))
  }
  assign("backend_shown", c(shown, key), envir = .arf_env)
  message(msg)
  invisible(TRUE)
}

# Decide which parallel backend forde() should use, and tell the user.
# Only meaningful when parallel = TRUE. Precedence:
#   1. explicit options(arf.backend = "foreach" | "mirai")  (validated)
#   2. active mirai daemons                                  -> "mirai"
#   3. otherwise                                             -> "foreach"
# For the foreach path we additionally report whether a real parallel backend
# is registered (>1 worker) or whether it will fall back to sequential.
# Returns one of "sequential", "foreach", "mirai".
arf_select_backend <- function(parallel) {
  if (!isTRUE(parallel)) {
    return("sequential")
  }
  mirai_ready <- requireNamespace("mirai", quietly = TRUE) &&
    requireNamespace("mori", quietly = TRUE) &&
    {
      st <- mirai::status()
      !is.null(st$connections) && st$connections >= 1L
    }
  dopar_workers <- if (requireNamespace("foreach", quietly = TRUE)) {
    foreach::getDoParWorkers()
  } else {
    1L
  }

  opt <- getOption("arf.backend", NULL)
  if (!is.null(opt)) {
    # Explicit user choice wins; validate and hard-check mirai readiness.
    backend <- match.arg(opt, c("foreach", "mirai"))
    if (backend == "mirai") {
      arf_check_mirai_ready()
    }
  } else if (mirai_ready) {
    backend <- "mirai"
  } else {
    backend <- "foreach"
  }

  # Report the effective backend: parallel = TRUE only, at most once per
  # distinct state per session, suppressible via options(arf.verbose = FALSE).
  if (backend == "mirai") {
    n <- mirai::status()$connections
    arf_backend_inform(
      paste0("arf: using 'mirai' backend (", n, " daemon",
             if (n != 1L) "s" else "", ")."),
      key = paste0("mirai:", n))
  } else if (dopar_workers > 1L) {
    arf_backend_inform(
      paste0("arf: using 'foreach' backend (", foreach::getDoParName(),
             ", ", dopar_workers, " workers)."),
      key = paste0("foreach:", dopar_workers))
  } else {
    arf_backend_inform(
      paste0("arf: parallel = TRUE but no parallel backend is registered; ",
             "computing sequentially. Register a foreach backend (e.g. ",
             "doParallel) or start mirai daemons via mirai::daemons()."),
      key = "sequential-fallback")
  }
  backend
}

# Load arf on all mirai daemons so workers can call its internals (forge/expct/
# lik workers use cforde/resample/post_x, unlike forde's self-contained per-tree
# workers). No-op on daemons where arf is already loaded (e.g. dev-loaded in
# tests; see tests/testthat/helper-mirai.R).
arf_load_on_daemons <- function() {
  mirai::everywhere(suppressMessages(loadNamespace("arf")))
  invisible(TRUE)
}

# Worker count for sizing step/fold chunks. Whichever parallel backend is active
# reports >1; the inactive one reports 1, so max() picks the right pool without
# re-deriving backend-selection precedence. Needed because stepsize was sized via
# foreach::getDoParWorkers() alone, which is 1 under a pure mirai backend (daemons
# set, no foreach registered) -> step_no 1 -> mirai never engages.
# ponytail: sizing hint only; the dispatch gate governs real parallelism, so a
# slight over/under-split when both backends are up is harmless.
arf_n_workers <- function() {
  mirai_conns <- if (requireNamespace("mirai", quietly = TRUE)) {
    st <- tryCatch(mirai::status(), error = function(e) NULL)
    if (!is.null(st$connections)) as.integer(st$connections) else 0L
  } else {
    0L
  }
  dopar <- if (requireNamespace("foreach", quietly = TRUE)) {
    foreach::getDoParWorkers()
  } else {
    1L
  }
  max(1L, mirai_conns, dopar)
}

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
# in the worker bodies tolerate this. No setalloccol() repair is needed
# because workers only read the shared tables and build fresh data.tables
# from the results.

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
