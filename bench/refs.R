# Ref resolution for the version-comparing benchmark: a ref spec becomes an
# installed library. See bench/README.md for the spec syntax.

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

# Modified TRACKED paths only. Untracked files are expected here: results/,
# logs/ and local experiment scripts all live under bench/.
bench_tree_dirty <- function(dir = ".") {
  out <- system2("git", c("-C", dir, "status", "--porcelain", "--untracked-files=no"), stdout = TRUE)
  if (!length(out)) {
    return(character())
  }
  trimws(substring(out, 4))
}

# A number from a tree nobody can reconstruct cannot be cited, so refuse it
# rather than label it: a number from a tree nobody can reconstruct cannot be
# cited, so refusing beats labelling.
bench_assert_clean_tree <- function(dir = ".", ignore = "bench/history.csv") {
  # history.csv is benchmark OUTPUT that the suite appends to and the author
  # commits afterwards. Its modification cannot make a measurement
  # irreproducible, and gating on it would make `make bench` one-shot per
  # commit once the file is tracked.
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

# The shared dependency layer is simply the session's own libpaths: every dep
# already resolves there, and installing them per ref would let their versions
# diverge between refs and confound the comparison with somebody else's
# performance change. A per-ref library therefore holds only arf, and children
# prepend it with `.libPaths(c(ref_lib, .libPaths()))`.
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

# Identical flags for every ref. Never add --no-byte-compile: it would change
# the thing being timed.
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

.bench_install_cran <- function(version, lib) {
  # The tarball CRAN ships is the faithful anchor, and a current release lives
  # at a different URL than an archived one, so try both.
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

# Idempotent: one job per cell would otherwise reinstall every ref for every
# cell, which is most of the wall clock of a local run (20 quick cells x 3 refs
# = 60 installs). A stamp file records which commit the library holds, so a
# second request for the same ref is free. Slurm jobs land on their own
# tempdirs, so there is nothing to race against.
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
    commit = id
  )
}

# R drops or partially matches an argument a ref does not have, so the call
# would succeed against that version's defaults and produce a plausible number
# for a different computation. Check up front instead.
bench_assert_args <- function(lib, calls) {
  # In a clean child, because this session may already have arf attached from a
  # different library and loadNamespace() would hand back that cached copy.
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

# Two specs can name the same commit, which is the default state of every run
# on main ("HEAD" and "main"). Comparing a ref with itself pays two installs
# and two measurement sets per cell to print a verdict on sampling noise.
# Identity for deduplication: the PACKAGE content, not the commit. A branch
# that only touches bench/ or the Makefile is the same arf as its base, and
# benchmarking both doubles a cluster run to prove 0.0%.
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
    # No package paths resolved: fall back to the commit rather than treating
    # every such ref as identical.
    return(.bench_ref_id(ref))
  }
  paste(tree, collapse = "\n")
}

# `baseline` survives a collapse. Deduplication keeps the first spec, and the
# default order is HEAD before main, so on the main branch itself the two
# resolve alike and "main" was the one dropped -- leaving bench_deltas() to
# stop on "baseline ref 'main' is not in the results".
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
