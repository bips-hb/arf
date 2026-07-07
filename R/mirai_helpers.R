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
#
# everywhere() round-trips all daemons on every call, which is dead cost for a
# workflow doing many small forge/expct/lik calls against one pool. Register the
# load once per pool and skip on repeat. The pool key is the dispatcher URL
# (status()$daemons), which is minted fresh by every daemons() call: a
# teardown+rebuild yields a new key so the cache self-invalidates, and daemons
# joining an existing pool auto-run the registered everywhere() expression, so
# one call per pool suffices. If status() is unavailable (key NULL) we fall back
# to the old always-load behavior.
arf_load_on_daemons <- function() {
  key <- tryCatch(mirai::status()$daemons, error = function(e) NULL)
  cached <- get0("arf_loaded_key", envir = .arf_env, ifnotfound = NULL)
  if (!is.null(key) && identical(key, cached)) {
    return(invisible(FALSE))
  }
  mirai::everywhere(suppressMessages(loadNamespace("arf")))
  assign("arf_loaded_key", key, envir = .arf_env)
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
# Chunks MUST be contiguous (sort()): rbindlist concatenates them in chunk
# order, so interleaved chunks would scramble positional output (forge/expct
# rows are per-evidence-row; forde is keyed so order-independent, but we can't
# rely on that here). See test "mirai preserves row order".
# `combine` stacks per-item results within a chunk and then across chunks. Default
# rbindlist returns a data.table (forde wants that; results are re-keyed). forge/
# expct pass a rbind-based combine so the object class matches their serial foreach
# .combine="rbind" output (data.table vs data.frame is otherwise a parity break;
# values are identical either way).
# CAUTION on closures: serializing a function serializes its enclosing
# environment. A closure defined inside a caller's frame drags that whole frame
# (params, evidence, training data, ...) into EVERY task, defeating mori
# sharing. Pass package-level functions (serialized as a namespace reference)
# or strip base-R-only closures to globalenv() before shipping.
# This is a documented mirai gotcha, see the Community FAQ:
# https://mirai.r-lib.org/articles/v07-questions.html
# The FAQ recommends carrier::crate() for the general case. We skip it here
# because our shipped functions are pure base R with all inputs as explicit
# arguments: there is nothing to crate, and environment(fn) <- globalenv()
# gets the same zero-payload serialization without adding a carrier dependency.
# Split trees into one contiguous block per worker. Contiguous (sort()) is
# load-bearing: rbindlist/c concatenate blocks in chunk order, so interleaved
# chunks would scramble positional output. Shared by arf_mirai_tree_map() and
# the prune dispatch in adversarial_rf(); the invariant lives here once.
arf_tree_chunks <- function(num_trees, n_workers) {
  n_chunks <- max(1L, min(as.integer(n_workers), num_trees))
  split(seq_len(num_trees),
        sort(rep(seq_len(n_chunks), length.out = num_trees)))
}

arf_mirai_tree_map <- function(num_trees, worker_fn, shared_args,
                               combine = data.table::rbindlist) {
  st <- mirai::status()
  n_workers <- max(1L, as.integer(st$connections))
  chunks <- arf_tree_chunks(num_trees, n_workers)
  chunk_runner <- function(trees, worker_fn, shared_args, combine) {
    parts <- lapply(trees, function(tr) {
      do.call(worker_fn, c(list(tr), shared_args))
    })
    combine(parts)
  }
  # base-R body, all inputs explicit args: strip so this frame (chunks,
  # shared_args, combine, ...) is not serialized into every task
  environment(chunk_runner) <- globalenv()
  res <- mirai::mirai_map(
    chunks,
    chunk_runner,
    .args = list(worker_fn = worker_fn, shared_args = shared_args,
                 combine = combine)
  )[]
  combine(res)
}

# Combine for forge/expct step results: matches serial foreach .combine="rbind"
# (same class, clean 1..n row.names). Package-level on purpose: an inline
# closure in forge()/expct() would serialize their whole frame, params
# included, into every task (see note above).
arf_rbind_steps <- function(parts) {
  r <- do.call(rbind, parts)
  rownames(r) <- NULL
  r
}
