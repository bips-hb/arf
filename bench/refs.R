# Ref resolution: a ref spec becomes an installed library. See bench/README.md
# for the spec syntax.

# "git:<ref>", "cran:<version>", or a bare git ref such as "HEAD" or "main".
bench_parse_ref <- function(spec) {
  if (!nzchar(spec)) {
    stop("empty ref spec")
  }
  if (!grepl(":", spec, fixed = TRUE)) {
    return(list(kind = "git", value = spec, label = spec))
  }
  kind <- sub(":.*$", "", spec)
  value <- sub("^[^:]*:", "", spec)
  if (!kind %in% c("git", "cran")) {
    stop("unknown ref kind: ", kind)
  }
  if (!nzchar(value)) {
    stop(kind, " ref needs a version or revision")
  }
  list(
    kind = kind,
    value = value,
    label = if (kind == "cran") paste0("cran-", value) else value
  )
}

# Modified tracked paths only; untracked files are expected under bench/.
bench_tree_dirty <- function(dir = ".") {
  out <- system2("git", c("-C", dir, "status", "--porcelain", "--untracked-files=no"), stdout = TRUE)
  if (!length(out)) {
    return(character())
  }
  trimws(substring(out, 4))
}

# A number from a tree nobody can reconstruct cannot be cited. The two CSVs
# are the suite's own output, written by every collect.
bench_assert_clean_tree <- function(dir = ".", ignore = c("bench/history.csv", "bench/timings.csv")) {
  dirty <- setdiff(bench_tree_dirty(dir), ignore)
  if (length(dirty)) {
    stop(
      "working tree has modified tracked files, so results would not be reproducible:\n  ",
      paste(dirty, collapse = "\n  "),
      "\nCommit them, including to a throwaway branch, then rerun.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# One shared dependency layer, the session's own libpaths: per-ref dependency
# libraries would let versions diverge and confound the comparison.
bench_assert_deps <- function() {
  deps <- c(
    "data.table",
    "ranger",
    "stringr",
    "truncnorm",
    "foreach",
    "doParallel",
    "mirai",
    "mori",
    "digest",
    "callr",
    "batchtools"
  )
  missing <- deps[!vapply(deps, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("missing benchmark dependencies: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

# Identical for every ref. Never add --no-byte-compile: it changes what is timed.
.bench_install_flags <- c("--no-docs", "--no-help", "--no-multiarch")

.bench_install_tarball <- function(tarball, lib) {
  args <- c("CMD", "INSTALL", .bench_install_flags, "-l", shQuote(lib), shQuote(tarball))
  status <- system2(file.path(R.home("bin"), "R"), args, stdout = FALSE, stderr = FALSE)
  if (status != 0L) {
    stop("R CMD INSTALL failed for ", basename(tarball), call. = FALSE)
  }
  invisible(TRUE)
}

.bench_install_git <- function(rev, lib) {
  wt <- file.path(tempdir(), paste0("arf-wt-", substr(rev, 1, 12), "-", Sys.getpid()))
  unlink(wt, recursive = TRUE)
  status <- system2("git", c("worktree", "add", "--detach", shQuote(wt), shQuote(rev)), stdout = FALSE, stderr = FALSE)
  if (status != 0L) {
    stop("git worktree add failed for ", rev, call. = FALSE)
  }
  on.exit(
    {
      system2("git", c("worktree", "remove", "--force", shQuote(wt)), stdout = FALSE, stderr = FALSE)
    },
    add = TRUE
  )
  commit <- system2("git", c("-C", shQuote(wt), "rev-parse", "--short", "HEAD"), stdout = TRUE)
  .bench_install_tarball(wt, lib)
  commit
}

# The tarball CRAN shipped is the faithful anchor; current and archived
# releases live at different URLs.
.bench_install_cran <- function(version, lib) {
  urls <- sprintf(
    c(
      "https://cran.r-project.org/src/contrib/arf_%s.tar.gz",
      "https://cran.r-project.org/src/contrib/Archive/arf/arf_%s.tar.gz"
    ),
    version
  )
  dest <- file.path(tempdir(), sprintf("arf_%s.tar.gz", version))
  ok <- FALSE
  for (u in urls) {
    ok <- isTRUE(tryCatch(
      utils::download.file(u, dest, quiet = TRUE, mode = "wb") == 0L,
      error = function(e) FALSE,
      warning = function(w) FALSE
    ))
    if (ok) break
  }
  if (!ok) {
    stop("could not download arf ", version, " from CRAN or its archive", call. = FALSE)
  }
  .bench_install_tarball(dest, lib)
  paste0("cran-", version)
}

.bench_ref_id <- function(ref) {
  if (ref$kind == "cran") {
    return(paste0("cran-", ref$value))
  }
  id <- suppressWarnings(system2("git", c("rev-parse", "--short", shQuote(ref$value)), stdout = TRUE, stderr = FALSE))
  if (!length(id)) {
    stop("cannot resolve git ref: ", ref$value, call. = FALSE)
  }
  id[1]
}

# Idempotent via a stamp file, so one install serves every cell.
bench_install_ref <- function(spec, root, force = FALSE) {
  ref <- bench_parse_ref(spec)
  lib <- file.path(root, ref$label)
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  id <- .bench_ref_id(ref)
  stamp <- file.path(lib, ".bench-stamp")
  cached <- !force &&
    dir.exists(file.path(lib, "arf")) &&
    file.exists(stamp) &&
    identical(readLines(stamp, warn = FALSE)[1], id)
  if (!cached) {
    if (ref$kind == "git") {
      .bench_install_git(ref$value, lib)
    } else {
      .bench_install_cran(ref$value, lib)
    }
    writeLines(id, stamp)
  }
  list(
    label = ref$label,
    lib = lib,
    arf_version = as.character(utils::packageVersion("arf", lib.loc = lib)),
    commit = id,
    backends = bench_ref_backends(lib)
  )
}

# arf ignores an `arf.backend` it does not know and runs its default, so a ref
# is measured only on the backends its own DESCRIPTION supports.
BENCH_BACKEND_DEPS <- c(
  sequential = NA_character_,
  foreach = "foreach",
  psock = "foreach",
  mirai = "mirai"
)

bench_ref_backends <- function(lib) {
  d <- read.dcf(file.path(lib, "arf", "DESCRIPTION"))
  cols <- intersect(c("Depends", "Imports", "Suggests"), colnames(d))
  deps <- trimws(sub("\\(.*", "", unlist(strsplit(d[1, cols], "[,\n]"))))
  names(BENCH_BACKEND_DEPS)[
    is.na(BENCH_BACKEND_DEPS) | BENCH_BACKEND_DEPS %in% deps
  ]
}

# A ref without the probe (results predating it) is kept.
.bench_usable_refs <- function(refs, backend) {
  can <- vapply(
    refs,
    function(r) is.null(r$backends) || backend %in% r$backends,
    logical(1)
  )
  if (!all(can)) {
    message(
      "  skipping ",
      paste(vapply(refs[!can], function(r) r$label, character(1)), collapse = ", "),
      " for backend ",
      backend,
      ": that version does not have it, so the cell would measure its default ",
      "path under another name"
    )
  }
  refs[can]
}

# R drops or partially matches an argument a ref lacks, so the call would run
# against that version's defaults. Checked in a clean child, because this
# session may have another arf attached.
bench_assert_args <- function(lib, calls) {
  have <- callr::r(
    function(lib, fns) {
      .libPaths(c(lib, .libPaths()))
      loadNamespace("arf", lib.loc = lib)
      stats::setNames(lapply(fns, function(f) names(formals(getExportedValue("arf", f)))), fns)
    },
    args = list(lib = lib, fns = names(calls))
  )
  for (fn in names(calls)) {
    missing <- setdiff(calls[[fn]], have[[fn]])
    if (length(missing)) {
      stop("arf in ", lib, " has no argument(s) ", paste(missing, collapse = ", "), " for ", fn, "()", call. = FALSE)
    }
  }
  invisible(TRUE)
}

# Identity for deduplication is package content, not the commit: a branch that
# only touches bench/ is the same arf as its base.
.bench_pkg_id <- function(ref) {
  if (ref$kind == "cran") {
    return(paste0("cran-", ref$value))
  }
  paths <- c("DESCRIPTION", "NAMESPACE", "R", "src", "inst", "man")
  tree <- suppressWarnings(system2(
    "git",
    c("ls-tree", "-r", shQuote(ref$value), "--", paths),
    stdout = TRUE,
    stderr = FALSE
  ))
  if (!length(tree)) {
    return(.bench_ref_id(ref))
  }
  paste(tree, collapse = "\n")
}

# Deduplicate by content; the baseline spec survives a collapse, since
# bench_deltas() needs it by name.
bench_resolve_refs <- function(specs, baseline = "main") {
  ids <- vapply(specs, function(s) .bench_pkg_id(bench_parse_ref(s)), character(1))
  prefer <- order(specs != baseline)
  specs <- specs[prefer]
  ids <- ids[prefer]
  first <- !duplicated(ids)
  for (i in which(!first)) {
    message(
      "Skipping ",
      specs[i],
      ": installs the same package content as ",
      specs[first][match(ids[i], ids[first])],
      ", so it would measure the same code twice."
    )
  }
  specs[first]
}
