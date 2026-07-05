# Parallelism backend benchmarks

Scripts here compare the parallel backends available to `forde()`:

- **sequential** — `forde(parallel = FALSE)`
- **foreach** — `forde(parallel = TRUE)` with a registered `foreach` backend
  (e.g. `doParallel`). Data is serialized/copied to each worker.
- **mirai** — `options(arf.backend = "mirai")` with `mirai::daemons()` set. Data
  is shared read-only via `mori`, so it is not copied per worker.

Shared logic (backend checks, data generation, PSS memory sampling, the
single-backend runner) lives in `bench-helpers.R`, sourced by the scripts below.
All backend packages (`doParallel`, `mirai`, `mori`) are required: every script
aborts up front — before the model fit — if any is missing, since a backend
comparison is meaningless without all of them.

Results are written to `bench/results/`, which is **gitignored** (and `bench/`
is in `.Rbuildignore`, so nothing here is committed or shipped).

## The scripts

| Script | Measures | Use for |
|---|---|---|
| `sweep.R` | time **and** peak memory over a grid of worker counts × sizes | **the cluster run** — backend scaling + the memory story |
| `mem-backends.R` | time + peak memory at a single config | a quick one-off memory snapshot |
| `bench-backends.R` | time only, via `bench::mark` (multi-iteration, `itr/sec`) | rigorous *steady-state* timing + a correctness gate |

```sh
# scaling sweep (the main cluster tool)
ARF_BENCH_WORKERS=1,2,4,8,16 ARF_BENCH_N=5000,20000 ARF_BENCH_TREES=100,200 \
  ARF_BENCH_DT_THREADS=1 Rscript bench/sweep.R

# single snapshot
ARF_BENCH_N=20000 ARF_BENCH_TREES=200 ARF_BENCH_WORKERS=8 Rscript bench/mem-backends.R

# rigorous timing + correctness check
ARF_BENCH_WORKERS=8 Rscript bench/bench-backends.R
```

`sweep.R` fits the ARF **once per (n, trees)** and reuses it across worker
counts. Grids are comma-separated env overrides (`ARF_BENCH_WORKERS`,
`ARF_BENCH_N`, `ARF_BENCH_TREES`, `ARF_BENCH_P`).

## Parallelism topology (important)

Run on **reserved/idle resources**, and mind the thread topology or the numbers
are confounded. `data.table` is multi-threaded by default, so without control
each worker spins up several threads and `N` workers oversubscribe the machine
(this also makes the "sequential" baseline secretly multi-threaded). The scripts
pin `data.table` to **one thread per worker** (main process, forked foreach
workers, and mirai daemons) via `ARF_BENCH_DT_THREADS` (default 1), so `N`
workers means `N` cores — the fair comparison.

The one-time `adversarial_rf()` fit happens *outside* the timed region, so
ranger's own threading affects only setup wall-time, not the reported backend
numbers — it is orthogonal to the backend comparison.

### Timing caveat: single-run vs steady-state

`sweep.R` and `mem-backends.R` time a **single** `forde()` call per cell, so
they include first-call/cold overhead (e.g. loading `data.table` on freshly
launched daemons, the first `mori::share`). `bench-backends.R` runs several
iterations on warm daemons and reports the median, i.e. **steady-state**
throughput. For "cost of one call" read the sweep; for "throughput once
warmed up" read `bench-backends.R`. They can legitimately disagree, especially
for mirai.

## Why PSS, and why manual sampling

Memory is measured as peak total **PSS (Proportional Set Size)** across this
user's R process tree, from `/proc/<pid>/smaps_rollup` on Linux. PSS divides
shared pages among the processes mapping them, so mori's shared data counts once
(split), not once per daemon — the fair way to compare a copy-per-worker backend
(foreach) against a shared-memory one (mirai+mori).

We sample manually because there is no better R option: `bench` can't profile
parallel code, and the `ps` package (cleaner, cross-platform) only exposes
**RSS**, which double-counts shared pages and would *understate* mori's
advantage. PSS is Linux-only; elsewhere the scripts fall back to RSS and label
the column accordingly.

Note the mori memory advantage is **scale-dependent**: at small data the fixed
per-daemon/mori overhead can exceed the per-worker copy it saves (mirai may use
*more* than foreach), while at large data it wins clearly. The sweep is designed
to expose that crossover.
