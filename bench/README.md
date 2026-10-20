# arf benchmark suite

`bench/arf-bench` measures runtime and peak memory of arf git refs and CRAN releases on one grid of cells and reports deltas against `main`.
A cell is one op at one size, backend and worker count.
Every ref runs inside the same job, interleaved round by round, so a comparison never crosses nodes.
`bench/` is in `.Rbuildignore`; nothing here ships with the package.

## Requirements

- [Rapp](https://github.com/r-lib/Rapp) launcher, once: `Rscript -e 'install.packages("Rapp"); Rapp::install_pkg_cli_apps("Rapp")'`.
- R packages in the session library: arf's own deps (`data.table`, `ranger`, `stringr`, `truncnorm`, `foreach`), plus `doParallel`, `mirai`, `mori`, `digest`, `callr`, `batchtools`.
  `run` aborts up front if any is missing.
  These are shared by every ref on purpose: per-ref dependency libraries would let versions diverge and confound the comparison.
- `git` for ref worktrees; network only for `cran:` refs.
- cgroup v2 memory accounting for the process (any slurm job with `ConstrainRAMSpace`, or a systemd user slice); `run` refuses to start without it.
- Cluster runs need slurm with the site batchtools config at `/etc/xdg/batchtools/config.R` (template, partition, qos, `max.concurrent.jobs`).
- `quarto` on `PATH` for `viz`.

## Running

Always from the repository root.
A tree with modified tracked files refuses to run, so every row traces to a commit; commit work in progress to a throwaway branch.
`history.csv` is exempt, it is the suite's own output.

    bench/arf-bench --help                   # commands
    bench/arf-bench run --help               # every option with its real default
    bench/arf-bench plan                     # cells, resources, cost; nothing runs

    make bench                               # quick tier, local, HEAD vs main
    bench/arf-bench run                      # same thing

    bench/arf-bench -t full -c slurm plan    # what a cluster run costs
    bench/arf-bench -t full -c slurm run     # submit and exit
    bench/arf-bench collect                  # newest registry to CSV + report.md
    bench/arf-bench resubmit --dry-run       # which cells the run lost
    bench/arf-bench viz                      # bench/viz.html
    bench/arf-bench diag                     # what this machine looks like to the sampler
    Rscript bench/tests.R                    # the suite's own tests

Narrow a run with `--ops`, `--sizes`, `--backends`, `--rowmode`, `--evidence`, `--max-cells`; cheapen a validation run with `--iters 1 --mem-reps 1`.
`--refs` overrides the ref list, `--anchors` adds the CRAN anchors to a quick run.

Output: `results/bench-<tier>-<timestamp>.csv` and `results/report.md`, gitignored.
`collect` also appends to two committed files: `timings.csv` (per-call times, used for walltime estimates) and `history.csv` (anchor rows for the release trend).
Commit them for runs worth keeping.
`collect` re-derives everything from the stored per-ref rows, so verdict rules can change and the report be regenerated without rerunning.

### Cluster

One job per cell, resources per cell.
Cores and memory are matched fractions of a node (1/16 to 1/1, 192 threads and 1.1 TB per node, so about 6 GB per thread); requesting them independently strands whichever is left over.
Threads are `2 x (workers + 2)`: two spare cores keep mirai's dispatcher and the sampler from starving.
Memory comes from measured peaks keyed on op, size, evidence count and row mode.
Walltime comes from `timings.csv` as `1.3 x mem_reps x n_refs x (iters x call_seconds + 120 s)`, 8 h where a cell has no timing yet; a cell without its own timing inherits the slowest cell sharing op, size, evidence and row mode, never scaled by worker count.
`--cpus`, `--mem` (total MB) and `--walltime` force one figure for every cell.

`collect` warns about jobs that neither finished nor errored.
`resubmit` sends them back with the request each carried and twice the walltime; walltime kills and OOM kills look identical in the registry, so pass `--mem-factor 2` when `sacct -j <batch.id> --format=State,ExitCode,MaxRSS` says memory.

## The grid

Ops: `adversarial_rf`, `forde`, `lik`, `forge`, `expct`.

| tier | n | trees | backends | workers | mem_reps |
|---|---|---|---|---|---|
| quick | 1e3, 1e4 | 10, 50 | sequential | - | 3 |
| full | 1e4, 5e4 | 200 | sequential, foreach, psock, mirai | 1, 2, 4, 8, 16 | 5 |

`forge` and `expct` appear in four variants: evidence 100 or 200 rows, row mode `separate` or `or`, with `n_synth` 1 or 100.
Evidence stays at 200 because the 0.2.5 anchor needs more than a node for `expct` at 1000.

Backends:

- `sequential`: `parallel = FALSE`.
- `foreach`: `arf.backend = "foreach"` with `doParallel::registerDoParallel(cores = w)`, forked workers.
- `psock`: same arf backend on a PSOCK cluster, non-fork workers, so fork's copy-on-write sharing is removed and every worker holds its own input copy.
- `mirai`: `arf.backend = "mirai"` with `mirai::daemons(w)`, inputs shared read-only through `mori`.

`data.table` and `ranger` are pinned to one thread in every process (parent, fork, PSOCK, daemon), recorded per row as `dt_threads` and `ranger_threads`.
Without that the sequential backend is secretly parallel.
A multi-threaded regime is a separate question this grid does not ask; if you add one, keep `workers x threads` well under the core count, mirai's dispatcher degrades badly when starved.

With `stepsize = 0` the unit of parallel work is `ceiling(n_evidence / workers)`, and parallelism is skipped when that yields one step.
So a worker sweep is "the default configuration at w workers", not the same work on more workers, and `w = 1` on forge/expct runs in-process on every backend.

## Refs

Default refs are `HEAD` and `main`, plus `anchors.csv` on the full tier.
A spec is a bare git ref, `git:<ref>`, or `cran:<version>`.
Git refs become worktrees installed with `R CMD INSTALL -l`; `cran:` installs the published tarball, which is what users actually had.
Each ref library holds only arf, prepended to the session libpaths.

Refs are deduplicated by package content (`DESCRIPTION`, `NAMESPACE`, `R/`, `src/`, `inst/`, `man/`), so a branch that only touches `bench/` is not measured twice.
A ref runs only on the backends its own `DESCRIPTION` supports; arf ignores an unknown `arf.backend` and runs its default, which would otherwise be recorded under the wrong name.
The cell log says which refs it dropped.

## Measurement

Each ref measurement runs in a fresh R process via `callr::r_bg`, never a fork: `mirai` starts background threads and forking after that is unsafe, and forked workers can delete the shared tempdir.
A child does one backend, loads the fixture, pins threads, times `--iters` calls of the op (median reported, result dropped between calls) and samples memory throughout.
Setup is untimed, so one-time daemon or cluster start never enters a timing, but `--iters 1` still carries first-call cold costs.

Memory metric: `cgroup-anon+shmem`, the `anon` plus `shmem` counters of the job's cgroup v2 `memory.stat`, sampled every 5 ms.
The kernel charges each page once per cgroup, so pages shared by fork workers count once, and `shmem` catches mori's segments in `/dev/shm`.
Every process a cell spawns inherits the cgroup.
There is no fallback metric; PSS and RSS count shared pages differently and would not be comparable.

Reported memory is `peak_delta_mb`: the peak minus a per-cell floor (R, arf, the fixture; about 480 MB at n=1e4 and 1.7 GB at n=5e4 on the cluster).
Peak and floor are measured in the same round, rounds interleaved and counterbalanced across refs, and the median of the paired differences is reported.
That pairing matters more than the replicate count.
Every row records `mem_limit_mb`; a peak within 95% of it gets `hit the memory cap, not a measurement` on both memory and time, since memory.high throttling distorts runtime too.

## Reading the report

Deltas are the ref relative to the baseline: `+50%` on a `cran-0.2.5` row means 0.2.5 uses 50% more than `main`.
Memory deltas at or above 3% are `real`; time deltas below 10% are `inconclusive`.
Guards, per comparison not per cell:

- `single sample, not a finding`: one memory replicate or one timing.
- `cell too small to resolve memory`: marginal smaller than the floor, which is every quick-tier cell, so `make bench` is a time check that also records memory.
- `too fast to time reliably`: either side under 0.5 s.
- `ref did not complete` / `no ref completed this cell`: a child died; the log carries the error.

Labels do not aggregate: several agreeing `single sample` cells are stronger evidence than one `real` cell, and no rule says so.
Read the CSV.

Each cell digests its result and the report prints `OUTPUT DIFFERS` where refs disagree.
Deterministic ops hash the rounded output; `forge` and `expct` run under a fixed seed and hash a rounded summary, because an exact hash also moves when a refactor merely consumes RNG draws in a different order; `adversarial_rf` is not hashed.
Parallel `forge` in `separate` mode ignores the seed, so its digest is unstable at `w >= 2` and no mismatch is flagged there.

`viz.qmd` renders the newest CSV (or `--csv`) as a delta table, a backend comparison by worker count, per-cell plots, a memory-resolution panel, and the release trend from `history.csv`.

## Files

| File | Role |
|---|---|
| `arf-bench` | the CLI; everything goes through it |
| `cells.R` | the grid, plus per-cell slurm memory, cores and walltime |
| `refs.R` | ref spec to installed library, dedup by content, dependency and clean-tree checks |
| `run-cell.R` | one cell against every ref, callr children, interleaved rounds |
| `registry.R` | batchtools registry, local or slurm, submit and resubmit per resource shape |
| `collate.R` | deltas, verdicts, `report.md`, `history.csv`, `timings.csv` |
| `bench-helpers.R` | cgroup sampler, fixture, backend setup in the child, digests |
| `viz.qmd` | the HTML overview |
| `anchors.csv` | CRAN releases added on the full tier |
| `tests.R` | the suite's tests |
