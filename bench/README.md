# Parallelism backend benchmarks

Scripts here compare the parallel backends available to `forde()`:

- **sequential** — `forde(parallel = FALSE)`
- **psock** — foreach on a PSOCK cluster (`--backend psock` in the grid): clean
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
and reports deltas against `main`.

The job unit is one cell against **every** ref, with the ref loop inside the
job. Making `ref` a grid dimension would scatter refs across nodes and
reintroduce the cross-node drift that same-allocation A/B exists to remove;
cells still parallelise across the cluster, each comparison stays on one node.

A dirty tree refuses to run: the suite aborts unless no tracked file is
modified, so every row traces to a commit by construction. Untracked files are
fine, and `history.csv` is exempt because it is the suite's own output. To
benchmark work in progress, commit it, including to a throwaway branch.

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

Refs default to `HEAD` and `main` on the quick tier, plus the anchors in
`anchors.csv` on the full tier, deduplicated by resolved commit so a run on
`main` does not compare `HEAD` with itself. `--anchors` adds them to a quick
run; `--refs` overrides entirely. The anchors answer a release-trend question
the quick tier cannot (its cells are below the memory resolution limit), so
paying for them on every mid-work check would add a third to the cost for
nothing.
Anchor rows are appended to `history.csv` automatically; committing them is
manual, for the runs worth keeping.

A ref spec is a bare git ref, `git:<ref>`, or `cran:<version>`. A git ref
becomes a `git worktree` installed with `R CMD INSTALL -l`, which works for
unpushed branches and needs no network. `cran:` installs the published tarball
from the CRAN archive: a git tag and the tarball CRAN actually shipped are not
always identical, and for a release anchor the tarball is the more faithful
answer to "what did users have". Use `git:` inside the current development line
and `cran:` for historical releases.

Dependencies are deliberately **not** per ref. One shared library holds
`data.table`, `ranger`, `mirai` and the rest, each per-ref library holds only
`arf`, and `R_LIBS` is `lib/<ref>:lib/shared`. Per-ref dependencies would let
versions diverge between refs and confound the comparison with somebody else's
performance change.

Refs are deduplicated by **package content** (`DESCRIPTION`, `NAMESPACE`, `R/`,
`src/`, `inst/`, `man/`), so a branch that only touches `bench/` is not
measured twice against its own base. When every ref collapses to one, the run
says so: it will still fill `timings.csv` but no delta can come out of it.

`collect` re-derives the CSV and report from the stored per-ref rows, so the
threshold and verdict rules can be revised and the report regenerated without
recomputing anything.

`viz.qmd` is the overview: render it with `bench/arf-bench viz` (or
`quarto render bench/viz.qmd`), which writes `bench/viz.html`. Sections: the
newest `results/bench-*.csv` as a delta table, sortable and filterable once it
is large; a **backend comparison** of `foreach` against `mirai` and `psock` by
worker count; the per-cell measurements as plots; a memory-resolution panel
saying which cells can carry a memory verdict at all; and the release trend
across anchors from `history.csv`. Point it at a specific run with
`QUARTO_BENCH_CSV`, or change the baseline with `QUARTO_BENCH_BASELINE`.

The cluster path needs no template argument: batchtools reads
`/etc/xdg/batchtools/config.R`, which already names the site template and sets
qos, partition and `max.concurrent.jobs`.

Cores, memory and walltime are requested per cell. Run
`bench/arf-bench -t full -c slurm plan` to see every shape and the total cost
before submitting:

    slurm:  30 cell(s) at 1/16  node =  12 cpus ( 6 cores),   67 GB,   8.0 h
    slurm:   9 cell(s) at 1/8   node =  24 cpus (12 cores),  134 GB,   8.0 h
    slurm:  19 cell(s) at 1/4   node =  48 cpus (24 cores),  269 GB,   3.5 h
    slurm: 169 cell(s) at 1/4   node =  48 cpus (24 cores),  269 GB,   8.0 h
    slurm:  48 cell(s) at 1/2   node =  96 cpus (48 cores),  537 GB,   8.0 h
    slurm:  16 cell(s) at 1/1   node = 192 cpus (96 cores), 1074 GB,  16.0 h
    total: 1021 node-hours; 80 of 352 cells have a measured timing

Cores and memory go together as matched fractions of a node (1/16, 1/8, 1/4,
1/2, whole). A node is 192 threads and 1152 GB, so about 6 GB per thread, and
requesting the dimensions independently strands whichever is left over.

The thread need is `2 x (workers + 2)`: the two spare cores are for mirai's
dispatcher and the 5 ms cgroup sampler, since a starved dispatcher makes mirai
lose comparisons it should win. `ncpus` counts hyperthreads, not cores.

The memory need comes from measured peaks, keyed on op, size, evidence count
**and** row mode, because one figure cannot serve the grid: an `expct` cell
needs 183 GB at `evidence = 100, separate` and over 545 GB at
`evidence = 1000, "or"`. Every row records `mem_limit_mb`, and a peak within
95% of it reads `hit the memory cap, not a measurement`.

Walltime is estimated per cell from `timings.csv`, a committed record of
measured per-call times that every `collect` updates, so estimates sharpen as
the suite runs. The model is
`safety x mem_reps x n_refs x (iters x call_seconds + overhead)`. A cell with no
timing of its own inherits the slowest cell sharing its op, size, evidence
count and row mode, and never scales by worker count: on `expct`, `foreach`
scales as `1/w` while `mirai` gets *slower* with more workers (413 s at 1
worker against 1104 s at 16), so no per-worker rule is safe in either
direction. `plan` warns about any cell whose estimate already exceeds the QoS
ceiling, since such a cell will be killed at that ceiling whatever is
requested.

`--cpus`, `--mem` and `--walltime` override with one figure for every cell
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

Both guards are **per comparison, not per cell**: a delta has two sides, so it
is invalid when the ref or its baseline is affected and unaffected when some
third ref in the same cell is. A ref whose child died reads `ref did not
complete`, or `no ref completed this cell` when both sides did, rather than
blaming the baseline.

These labels are conservative per comparison and cannot aggregate. Several
independent cells each labelled `single sample` while agreeing closely is
stronger evidence than one cell labelled `real`, and no rule here will say so;
read the CSV.

### Output differences

Each cell hashes its own result, and the report prints `OUTPUT DIFFERS` beside
any comparison whose refs disagree, marking that delta as not apples-to-apples.
Nothing aborts: an intentional fix changes output by design, and only a human
knows which side is correct.

The digest differs by op class. Deterministic ops (`forde` parameters, `lik`
values) hash the rounded output directly. Stochastic ops (`forge`, `expct`) run
under a fixed seed and hash a rounded summary vector of column means and
standard deviations, because an exact hash there also moves when a refactor
merely consumes RNG draws in a different order. `adversarial_rf` is not hashed.
The report says which kind produced a flag, so a stochastic mismatch reads as
"look at this" rather than "this is broken".

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

| File | Role |
|---|---|
| `arf-bench` | the CLI: `run`, `collect`, `viz`, `plan`, `diag`. Everything goes through it. |
| `cells.R` | the grid as data, plus the per-cell slurm memory, cores and walltime |
| `refs.R` | a ref spec to an installed library, deduplicated by package content |
| `run-cell.R` | one cell against every ref, in callr children, interleaved and counterbalanced |
| `registry.R` | the batchtools registry, local or slurm, submitted per resource shape |
| `collate.R` | deltas, verdicts, `report.md`, `history.csv`, `timings.csv` |
| `bench-helpers.R` | the cgroup sampler, fixture generation, output digests |
| `tests.R` | the suite: `Rscript bench/tests.R` |

Run `bench/arf-bench --help`, or `bench/arf-bench <subcommand> --help`, for the
options and their defaults. `plan` shows what a run would submit without
submitting it.

## Parallelism topology (important)

Run on **reserved/idle resources**. Two libraries also thread internally and
confound the backend comparison if left uncontrolled: `data.table` (the per-tree
work) and `ranger` (forde's `terminalNodes` prediction).

Both are pinned to **one thread** for every measurement: `data.table` via
`setDTthreads()` on the main process, forked workers and mirai daemons alike,
and `ranger` via `options(ranger.num.threads)`, which forde's `predict`
inherits. That is the unconfounded scaling baseline, and without it
`data.table`'s default multi-threading makes even the `sequential` backend
secretly parallel. The `dt_threads` and `ranger_threads` columns record it per
row, so a future run at a different setting is never silently compared against
these.

A realistic multi-threaded regime is a separate question this grid does not
ask: every row here is single-threaded by construction. If you add one, mind
that mirai degrades *catastrophically* when oversubscribed, because its
per-`forde` dispatch over nanonext gets CPU-starved by the compute threads
(measured: 900 s against foreach's 53 s at 16 workers x 10 threads on 192
cores). Keep `workers x threads` under the core count with headroom, and export
`OMP_WAIT_POLICY=passive` so idle OpenMP threads sleep rather than spin.

Core allocation is handled for you on slurm: a cell asks for `2 x (workers + 2)`
hyperthreads, rounded up to a node fraction, so the worker pool always has
cores to spare for coordination.

### Timing caveat: cold vs steady-state

`--iters` sets how many op calls each child makes, and the cell reports their
median. At `--iters 1` that single call includes first-call cold overhead
(loading `data.table` on freshly launched daemons, the first `mori::share`), and
mirai in particular looks worse cold than warm. Higher values report warmed-up
throughput at proportionally more walltime. The child is timed only around the
op itself; one-time daemon and cluster setup is excluded from every iteration.

`--iters` also multiplies the walltime estimate, so `plan` is the cheap way to
see what raising it costs before submitting.

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
*more* than foreach), while at large data it wins clearly. The grid spans both
sizes at every worker count so the crossover shows up rather than being averaged
away, and the memory-resolution panel says which cells are large enough to carry
a verdict in the first place.
