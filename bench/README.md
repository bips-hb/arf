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

## On the memory columns

`bench::mark` reports `mem_alloc` and gc counts, but these track allocations in
the **main R process only**. The parallel backends do their per-tree work in
separate worker/daemon processes, so `mem_alloc` does **not** capture the
per-worker data copies that motivate mirai+mori — it can even make the parallel
backends look *lighter* than sequential (the main process just collects
results). Use `mem_alloc`/gc for main-process signal and wall-clock time for
throughput.

To measure the cross-process memory footprint — the actual mirai/mori advantage
— you need OS-level RSS/PSS of the whole R process tree (e.g. sampling
`/proc/<pid>/smaps_rollup` on Linux). That is intentionally out of scope for
this `bench`-based script; add a separate memory harness if/when that comparison
is needed.
