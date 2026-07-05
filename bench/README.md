# Parallelism backend benchmarks

Scripts here compare the parallel backends available to `forde()`:

- **sequential** — `forde(parallel = FALSE)`
- **foreach** — `forde(parallel = TRUE)` with a registered `foreach` backend
  (e.g. `doParallel`). Data is serialized/copied to each worker.
- **mirai** — `options(arf.backend = "mirai")` with `mirai::daemons()` set. Data
  is shared read-only via `mori`, so it is not copied per worker.

The point of interest is **memory**: the foreach path copies the training data
(and forest) into every worker, whereas the mirai+mori path shares one copy. The
benchmark therefore reports peak memory as well as wall-clock time.

## Metric: PSS, not RSS

Memory is measured as the **peak total PSS (Proportional Set Size) across the R
process tree** on Linux (`/proc/<pid>/smaps_rollup`). PSS divides shared pages
among the processes that map them, so a page shared by N workers counts once
(split N ways) rather than N times. This is the fair way to compare a
copy-per-worker backend against a shared-memory one — RSS would double-count the
shared pages and understate mori's advantage's inverse (i.e. overstate its cost).

This is Linux-only. On other platforms the memory column is reported as `NA`.

## Running

```sh
Rscript bench/bench-backends.R                 # defaults
Rscript bench/bench-backends.R 5000 40 100 4   # n, p, num_trees, n_workers
```

Requires `doParallel` for the foreach path and `mirai` + `mori` for the mirai
path; each backend is skipped (with a note) if its packages are unavailable.

Run in an otherwise-idle session: the sampler sums PSS over the R process
subtree, so concurrent unrelated R work in the same tree would be attributed.

`bench/` is listed in `.Rbuildignore`, so nothing here ships in the package.
