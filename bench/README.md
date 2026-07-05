# Parallelism backend benchmarks

`bench-backends.R` compares the parallel backends available to `forde()` across
a grid of task sizes, using the [`bench`](https://bench.r-lib.org) package
(`bench::press` over configurations, `bench::mark` per backend):

- **sequential** — `forde(parallel = FALSE)`
- **foreach** — `forde(parallel = TRUE)` with a registered `foreach` backend
  (e.g. `doParallel`). Data is serialized/copied to each worker.
- **mirai** — `options(arf.backend = "mirai")` with `mirai::daemons()` set. Data
  is shared read-only via `mori`, so it is not copied per worker.

Each backend is included only if its packages are available (`doParallel` for
foreach; `mirai` + `mori` for mirai). `bench::mark`'s `check` runs a tolerant
`all.equal()` across backends, so the benchmark doubles as a correctness gate.

## Running

```sh
Rscript bench/bench-backends.R
ARF_BENCH_WORKERS=8 Rscript bench/bench-backends.R
```

Results (both the raw `bench_mark` object and a flat CSV summary) are written to
`bench/results/`, which is **gitignored** — nothing here is committed or shipped
(`bench/` is also in `.Rbuildignore`).

`bench-backends.R` reports **time only** (`bench::mark(memory = TRUE)` errors on
parallel code, and its `mem_alloc` sees only the main process anyway). Memory is
measured by a separate script.

## `mem-backends.R` — peak memory

```sh
Rscript bench/mem-backends.R
ARF_BENCH_N=20000 ARF_BENCH_TREES=200 ARF_BENCH_WORKERS=8 Rscript bench/mem-backends.R
```

Reports peak total memory across this user's R process tree while `forde()`
runs, per backend, and writes a CSV to `bench/results/`.

### Why PSS, and why manual sampling

The metric is **PSS (Proportional Set Size)** on Linux, read from
`/proc/<pid>/smaps_rollup`. PSS divides shared pages among the processes mapping
them, so mori's shared data counts once (split), not once per daemon — the fair
way to compare a copy-per-worker backend (foreach) against a shared-memory one
(mirai+mori).

We sample this manually rather than via a package because there is no better R
option for it: `bench` can't profile parallel code, and the `ps` package
(cleaner, cross-platform) only exposes **RSS**, which double-counts shared pages
and would *understate* mori's advantage. PSS is Linux-only; on other platforms
`mem-backends.R` falls back to RSS and labels the column accordingly.
