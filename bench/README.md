# Parallelism backend benchmarks

Scripts here compare the parallel backends available to `forde()`:

- **sequential** — `forde(parallel = FALSE)`
- **foreach** — `forde(parallel = TRUE)` with a registered `foreach` backend
  (e.g. `doParallel`). Data is serialized/copied to each worker.
- **mirai** — `options(arf.backend = "mirai")` with `mirai::daemons()` set. Data
  is shared read-only via `mori`, so it is not copied per worker.

Shared logic (backend checks, data generation, memory sampling, the
subprocess runner) lives in `bench-helpers.R`, sourced by the scripts below.
All backend packages (`doParallel`, `mirai`, `mori`) are required: every script
aborts up front — before the model fit — if any is missing, since a backend
comparison is meaningless without all of them.

Results are written to `bench/results/`, which is **gitignored** (and `bench/`
is in `.Rbuildignore`, so nothing here is committed or shipped).

## Process isolation (why, and how)

Each backend measurement runs in its **own fresh R process**, spawned with
`callr::r_bg` (which `posix_spawn`s a clean R, it does **not** fork). This is
required for correctness, not just tidiness: `mirai`/`nanonext` start background
threads, and forking afterwards is unsafe (SIGILL); forked workers can also
delete the shared session tempdir on exit. Running everything in one process —
and forking a sampler — is what crashed earlier. Now a child does exactly one
backend (mirai **or** fork, never both), and the parent samples the child's
`/proc` subtree without forking.

## The scripts

| Script | Measures | Use for |
|---|---|---|
| `sweep.R` | time **and** peak memory over a grid of worker counts × sizes | **the cluster run** — backend scaling + the memory story |
| `mem-backends.R` | time + peak memory at a single config | a quick one-off snapshot |

`run-sweep.sh` is the convenience launcher for a cluster node:

```sh
bench/run-sweep.sh                 # clean: 1 thread per worker (fair scaling baseline)
bench/run-sweep.sh realistic       # realistic: data.table + ranger at 10 threads each
bench/run-sweep.sh realistic 8     # ...with 8 threads
```

It just exports the `ARF_BENCH_*` variables and runs `sweep.R`; any variable you
set in the environment overrides its defaults. Or call the scripts directly:

```sh
ARF_BENCH_WORKERS=1,2,4,8,16 ARF_BENCH_N=5000,20000 ARF_BENCH_TREES=100,200 \
  ARF_BENCH_DT_THREADS=1 ARF_BENCH_RANGER_THREADS=1 ARF_BENCH_ITERS=1 Rscript bench/sweep.R

ARF_BENCH_N=20000 ARF_BENCH_TREES=200 ARF_BENCH_WORKERS=8 Rscript bench/mem-backends.R
```

`sweep.R` fits the ARF **once per (n, trees)**, caches it to disk, and each child
reads it (no refit). Grids are comma-separated env overrides (`ARF_BENCH_WORKERS`,
`ARF_BENCH_N`, `ARF_BENCH_TREES`, `ARF_BENCH_P`).

## Parallelism topology (important)

Run on **reserved/idle resources**. Two libraries also thread internally and
confound the backend comparison if left uncontrolled: `data.table` (the per-tree
work) and `ranger` (forde's `terminalNodes` prediction). Both are pinned via env
vars — `ARF_BENCH_DT_THREADS` (`data.table`, on main + forked workers + mirai
daemons) and `ARF_BENCH_RANGER_THREADS` (`ranger`, via `options(ranger.num.threads)`
which forde's `predict` inherits).

### Two regimes: clean vs realistic

- **clean** (`run-sweep.sh`, threads = 1): each worker is single-threaded, so
  `N` workers means `N` cores. This is the *unconfounded scaling baseline* —
  without it, `data.table`'s default multi-threading even makes "sequential"
  secretly parallel.
- **realistic** (`run-sweep.sh realistic`, threads = N): `data.table`/`ranger`
  multi-thread as in real-world use, so each worker uses several cores.

**Do not oversubscribe.** Keep `workers × threads` **under the core count, with
headroom** — the sweep prints an `OVERSUBSCRIBED` warning otherwise. Beyond the
obvious thrashing, mirai degrades *catastrophically* when oversubscribed: its
per-`forde` dispatch/collect over nanonext gets CPU-starved by the compute
threads (in one 16×10 run on 192 cores, mirai took 900s vs foreach's 53s). So in
realistic mode either lower the worker grid or the thread count; e.g. on C cores
pick `max(workers) × threads ≲ C`. `run-sweep.sh` also exports
`OMP_WAIT_POLICY=passive` so `data.table`'s idle OpenMP threads sleep instead of
spinning (spinning inflates load average and starves mirai's coordination).

Compare the two regimes: if mirai's memory advantage or the speed ordering flips,
that is a finding worth reporting. The `dt_threads`/`ranger_threads` columns in
the CSV record which regime produced each row.

### Timing caveat: cold vs steady-state

With `ARF_BENCH_ITERS=1` (default) each cell times a **single** `forde()` call,
which includes first-call/cold overhead (loading `data.table` on freshly
launched daemons, the first `mori::share`). Set `ARF_BENCH_ITERS` higher to run
several `forde()` calls per child and report the **median**, i.e. warmed-up
throughput — mirai in particular looks worse cold than warm. (The child is
timed only around `forde()`; the one-time daemon/cluster setup is excluded from
every iteration.)

## How memory is measured

Preferred metric: **cgroup-anon** — the `anon` counter of the cgroup v2
`memory.stat`, sampled at 5ms and reported as a delta against the pre-cell
baseline. The kernel charges each page once per cgroup, so COW pages shared by
fork workers and mori-shared pages count exactly once — the fair way to compare
a copy-per-worker backend (foreach) against a shared-memory one (mirai+mori),
and *exact* where summed PSS only approximates. Every process a cell spawns
inherits the cgroup, so short-lived fork workers are always counted, and `anon`
excludes page cache noise. Slurm gives each job its own cgroup, so cluster
measurements are isolated by construction.

Fallback (no cgroup v2 memory accounting): peak summed **PSS** across the
cell's process group, from `/proc/<pid>/smaps_rollup`, sampled at 50ms. Beware
its limits: each read forces a kernel VMA walk, so sweeps slow down as memory
grows, and workers that spawn and die between sweeps are missed. Last resort
without smaps_rollup is **RSS**, which double-counts shared pages and would
*understate* mori's advantage.

We sample manually because there is no better R option: `bench` can't profile
parallel code, and the `ps` package (cleaner, cross-platform) only exposes RSS.
The `metric` column in every CSV records which measurement was used — don't
compare peak_mb across CSVs with different metrics.

Note the mori memory advantage is **scale-dependent**: at small data the fixed
per-daemon/mori overhead can exceed the per-worker copy it saves (mirai may use
*more* than foreach), while at large data it wins clearly. The sweep is designed
to expose that crossover.
