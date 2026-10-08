# Parallelism backend benchmarks

Scripts here compare the parallel backends available to `forde()`:

- **sequential** — `forde(parallel = FALSE)`
- **psock** — foreach on a PSOCK cluster (`ARF_BENCH_BACKENDS=psock`): clean
  non-fork workers, useful to separate fork's parent-heap duplication from the
  per-worker input copies that mori removes
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

## Version comparison

`bench/arf-bench` compares refs rather than backends: it installs each ref into
its own library, measures the same grid against all of them in one allocation,
and reports deltas against `main`. See `DESIGN.md` for why the job unit is one
cell against every ref, and why a dirty tree refuses to run.

It is an [Rapp](https://github.com/r-lib/Rapp) script, so the defaults in
`--help` are the script's own top-level assignments and cannot drift from the
code. Install the launcher once:

    Rscript -e 'install.packages("Rapp"); Rapp::install_pkg_cli_apps("Rapp")'

Then:

    bench/arf-bench --help                  # commands
    bench/arf-bench run --help              # every option, with its real default
    bench/arf-bench plan                    # cells, executions and slurm sizing, nothing run

    make bench                              # quick tier, local
    bench/arf-bench run                     # the same thing

    bench/arf-bench -t full -c slurm plan   # what the cluster run would cost
    bench/arf-bench -t full -c slurm run    # submit and exit
    bench/arf-bench collect                 # collate the newest registry

Stage a big run with `--ops` and `--max-cells`, and cheapen a validation run
with `--iters 1 --mem-reps 1`:

    bench/arf-bench -t full -c slurm run -o forde -n 2 --iters 1 --mem-reps 1

Refs default to `HEAD`, `main`, and the anchors in `anchors.csv`, deduplicated
by resolved commit, so a run on `main` does not compare `HEAD` with itself.
Anchor rows are appended to `history.csv` automatically; committing them is
manual, for the runs worth keeping.

`collect` re-derives the CSV and report from the stored per-ref rows, so the
threshold and verdict rules can be revised and the report regenerated without
recomputing anything.

`viz.qmd` is the overview: render it with `bench/arf-bench viz` (or
`quarto render bench/viz.qmd`), which writes `bench/viz.html`. Sections: the
latest version comparison (sortable and filterable once it is large), a
**backend comparison** of `foreach` against `mirai` and `psock` by worker
count, per-cell measurements, which cells can resolve memory at all, the
release trend across anchors, and any earlier `sweep-ops-*.csv` runs. It shows the newest `results/bench-*.csv` as a delta table and
per-cell plots, a memory-resolution panel (which cells can carry a memory
verdict at all), the release trend from `history.csv`, and the legacy
`sweep-ops-*.csv` backend comparison when those files are present. Point it at
a specific run with `QUARTO_BENCH_CSV`, or change the baseline with
`QUARTO_BENCH_BASELINE`.

The cluster path needs no template argument: batchtools reads
`/etc/xdg/batchtools/config.R`, which already names the site template and sets
qos, partition and `max.concurrent.jobs`.

Cores and memory are requested together as matched fractions of a node
(1/16, 1/8, 1/4, 1/2, whole), per cell. A node is 192 threads and 1152 GB, so
6 GB per thread; requesting the dimensions separately strands whichever is left
over. `bench/arf-bench -t full -c slurm plan` prints the shapes:

    slurm:  30 cell(s) at 1/16  node =  12 cpus ( 6 cores),   72 GB
    slurm:   9 cell(s) at 1/8   node =  24 cpus (12 cores),  144 GB
    slurm: 185 cell(s) at 1/4   node =  48 cpus (24 cores),  288 GB
    slurm: 128 cell(s) at 1/2   node =  96 cpus (48 cores),  576 GB

The thread need is `2 x (workers + 2)`: the two spare cores are for mirai's
dispatcher and the 5 ms cgroup sampler, since a starved dispatcher makes mirai
lose comparisons it should win. The memory need comes from measured peaks,
because one figure cannot serve the grid: 64 GB OOM-killed an `n = 5e4` cell
while 256 GB was pinned at its own ceiling by `expct`. Every row records
`mem_limit_mb`, and a peak within 95% of it reads `hit the memory cap, not a
measurement`. `--cpus` and `--mem` override with one figure for every cell
(`--mem` is TOTAL MB, not `mem_per_cpu`, which the site defaults conflict
with).

`bench/arf-bench diag` prints what a machine looks like to the sampler: the
cgroup it reads, every limit up the hierarchy, the real `MemTotal`, and whether
the fixture lands on tmpfs. Run it under `sbatch` when a peak looks
implausible.

The three backend labels cover two arf backends: `foreach` and `psock` both set
`arf.backend = "foreach"` and differ in the foreach adapter, forked workers
versus a PSOCK cluster. That separates fork's parent-heap duplication from the
per-worker input copies `mori` removes.

Memory deltas are taken on `peak_delta_mb`, the peak minus a measured per-cell
floor (an R interpreter plus `arf` plus the fixture, about 200 MB here).
Peak and floor are measured in the same round and subtracted there, one round
per replicate interleaved and counterbalanced across refs, and the median of
those paired differences is reported. That matters more than the replicate
count: with the peak replicated but the floor measured once, same-commit spread
at `n = 1e3` stayed between 24% and 42% no matter how many replicates were
taken.

`--mem-reps` defaults to 3 in the quick tier and 5 in the full tier. A single
replicate is never reported as a finding, and neither is a cell whose marginal
is smaller than its floor: that is below the instrument's resolution and reads
`cell too small to resolve memory`. On current sizes that is every quick-tier
cell, so `make bench` is a time check that also records memory; memory verdicts
come from the full tier's large cells.

Time gets the same treatment at the other end of the scale: a comparison where
either side is faster than 0.5 s reads `too fast to time reliably`, because a
percentage on a 62 ms call is jitter. Measured on a real quick-tier run, the
only `real` time verdict under a second was the fastest cell in the grid.

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

Preferred metric: **cgroup-anon+shmem** — the `anon` plus `shmem` counters of
the cgroup v2 `memory.stat`, sampled at 5ms and reported as a delta against the
pre-cell baseline. The kernel charges each page once per cgroup, so COW pages
shared by fork workers count exactly once — the fair way to compare a
copy-per-worker backend (foreach) against a shared-memory one (mirai+mori), and
*exact* where summed PSS only approximates. `shmem` is required because mori
places its regions in `/dev/shm` (tmpfs), which the kernel files under `shmem`
rather than `anon`; results labelled plain `cgroup-anon` in older CSVs omit
that one shared copy of the inputs and understate mirai slightly. Every process a cell spawns
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
