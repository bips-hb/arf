# Benchmark suite design

Why this exists and why it is shaped this way.
`README.md` next to it covers how to run the thing, the memory metric, and the parallelism topology.
This file covers the version-comparison architecture on top of that.

## Purpose

Make claims like "this PR reduces peak memory 10% against main" and "runtime has not regressed since v0.2.5" reproducible, reviewable, and public.
Two audiences: a PR reviewer who wants one table, and the maintainers who want to see runtime and memory drift across releases.
Correctness outranks both: a performance number from a ref that computes something different is worse than no number.

## Decisions

### Comparison happens twice, from one producer

A/B in one allocation answers "is this PR better than main", because the deltas are same node, same kernel, same R, same dependency versions.
A committed history answers "what has happened across releases", which no single allocation can show.
Both read the same schema: every row carries a `ref` column, and `history.csv` is an append of the anchor rows from ordinary A/B runs.
One producer, two consumers, no second code path to keep aligned.

### Refs are anchors, not commits

Comparing every commit is wasteful and says little, so the ref set is deliberate and small: the working ref (`HEAD`), `main`, and a declared list of **anchors** in `anchors.csv`.
Anchors are git tags and CRAN releases, the points where the package's behaviour was actually published and where a regression matters.
The list is curated rather than derived from every tag, so adding an anchor is a decision with a cost someone chose to pay.

This has a consequence worth stating, because it is better than the usual arrangement.
Anchors are re-runnable at any time, so the whole release trend can be regenerated **on one node in one allocation** instead of stitched together from runs taken months apart on whatever hardware was free.
`history.csv` then stops being the source of truth and becomes a cache: it saves re-measuring old anchors on every run, and can be rebuilt from scratch whenever the numbers look suspect or the metric changes.

### batchtools orchestrates, the cgroup sampler measures

batchtools replaces `submit-ops.sh`: per-cell resources, retries, job status, and `reduceResultsDataTable()` for collation.
The quick tier uses its local cluster functions and the full tier uses slurm, so both tiers run the same code and cannot drift apart.
It does not supply the memory metric.
The `cgroup-anon+shmem` sampler in `bench-helpers.R` stays exactly as it is, inside each cell.
The earlier note that a batchtools migration "needs cgroup memory.peak, not MaxRSS" conflated scheduling with accounting: slurm's MaxRSS is indeed unusable here, which is why the metric never comes from the scheduler.

A constraint follows from this.
batchtools' natural unit is one job per grid row, but if `ref` is a grid dimension the refs scatter across nodes and reintroduce the cross-node drift that A/B exists to avoid.
So the job unit is `(op, cell)` and the ref loop runs *inside* the job.
Cells still parallelize across the cluster; each comparison stays on one node.

The ref *libraries*, however, are built once by the orchestrator into `bench/lib` on the shared filesystem, not inside each job.
Installing per job would have every job run `git worktree add` against the single shared `.git`, which at 352 cells means hundreds of concurrent mutations of `.git/worktrees`, plus a thousand redundant installs.
This is safe only because arf is pure R, with no `src/` and no `NeedsCompilation`, so an installed tree is architecture-independent.
A package with compiled code, or a cluster with heterogeneous nodes, would have to install inside the job and accept the cost.

### Rejected: targets with crew.cluster

targets sells a dependency graph plus change-based invalidation.
A benchmark grid has no interesting graph, and invalidation is actively hostile: inputs are identical by design across refs and across machines, and deliberately re-running an unchanged cell is the entire point.
Adopting it means setting `cue = "always"` everywhere, which switches off the feature it was adopted for.
crew workers are also persistent and reused, which collides with the fresh-process-per-cell requirement that the mirai-then-fork SIGILL, leftover daemons, and the shared session tempdir impose.
It would pay for the collate-and-report stage, which `viz.qmd` already covers.

### Rejected: the bench package

It cannot profile parallel code, and `mem_alloc` reports R-level allocation from the profiler: single process, blind to fork copies, to daemons, and to mori's `/dev/shm` regions.
Those are the quantities under study.
Its timing features duplicate the existing `ARF_BENCH_ITERS` medians.

### A dirty tree refuses to run

The benchmark aborts unless the working tree has no modifications to tracked files; untracked files are fine, since `results/`, `logs/` and local notes live there.
`bench/history.csv` is exempt: it is benchmark output that the suite appends to and the author commits afterwards, so its modification cannot make a measurement irreproducible, and gating on it would make `make bench` succeed once per commit and refuse thereafter.
Recording a dirtiness flag was the alternative and is worse: a number from a tree nobody can reconstruct cannot be cited, so the useful move is to refuse it rather than to label it.
Every row therefore traces to a commit by construction.
To benchmark work in progress, commit it, including to a throwaway branch.

### Replication is asymmetric

Wall time is noisy on shared nodes and peak memory turned out to be noisy too, so the original "memory is near-deterministic, one sample is enough" assumption does not survive measurement (see "What measurement showed" below).
Wall time on shared nodes is not, so the report takes a median over 5 iterations with min and max, and prints time deltas below 10% as inconclusive rather than as findings.

Timings come from one cell process: a cell spawns one `callr` child, calls the op 5 times inside it (`ARF_BENCH_ITERS`, which already works this way), and reports the median of those 5 timings.
Peak memory is the cgroup maximum over that whole child, which is one sample per child however many iterations run inside it.

### What measurement showed, and what changed

Two things were measured on 2026-10-05 and both changed the design.

Peak memory is dominated by a floor of about 200 MB: an R interpreter, `arf` and its dependencies, and the fixture, none of which any PR will change.
A percentage taken on the raw peak is therefore diluted four- to fivefold, and a genuine 10% regression in an op's own allocation renders as under 2%, below the 3% line.
So every cell now also measures its floor, with a child that loads the ref and reads the fixture but never runs the op, and `peak_delta_mb = peak_mb - floor_mb` is what deltas are computed on.

Peak memory is then not reproducible from one sample, because R's peak depends on when the garbage collector decides to grow the heap rather than on the algorithm alone.
Re-running identical code against itself gave marginal deltas of -2.7% to -11.6% at `n = 1e3` and as much as -38.6% at `n = 1e4`.
The raw peak looked stable only because the 200 MB floor was hiding this.
What made the measurement usable was the estimator, not the replicate count.
Replicating the cell peak alone did nothing: at `n = 1e3` the same-commit spread stayed between 24% and 42% whether one, two, three or five replicates were taken.
Three mistakes were stacked in that first version.
The floor was measured once while the peak was replicated, so un-replicated floor noise landed undiluted in a difference of two nearly equal numbers.
All of a ref's replicates ran before the next ref's, so any drift over the cell's lifetime fell entirely on the later ref, which showed up as a one-sided outlier.
And taking a minimum separately on each end of a subtraction inflates the spread rather than reducing it.

The version in the code measures peak and floor in the same round, subtracts them there, interleaves one round per replicate across all refs, and takes the median of the per-round marginals.
Same-commit spread at `n = 1e3` then fell to 1.4%, under the 3% threshold, so a 10% claim clears the noise by about sevenfold.

`mem_reps` is therefore a tier property: 3 for quick, which keeps it affordable on bertha, and 5 for full, where the numbers get published and the cluster has the budget.
`ARF_BENCH_MEM_REPS` overrides both for a one-off.
A single replicate is never reported as a finding: its verdict reads `single sample, not a finding`.

Replication has a resolution limit that no replicate count overcomes.
The marginal is a difference against a floor of roughly 200 MB, so when an op allocates much less than the floor the shared-cgroup delta cannot resolve it.
Same-commit spread against the ratio of marginal to floor: 0.29 gave 14-23%, 0.73 gave 4.2%, and 2.0 gave at most 3.6%.
So a cell whose marginal is below its floor carries no memory verdict at all, and reads `cell too small to resolve memory`.
Every current quick-tier cell falls under that line, which is consistent with the tier split: quick catches major regressions and reports time, and published memory numbers come from the full tier's large cells.
Raising the quick tier's `n` until the marginal clears its floor, or measuring per-child PSS instead of a shared-cgroup delta, would both lift the limit and neither is done here.

### Output differences are flagged, never fatal

Each cell hashes its own result and the report prints `OUTPUT DIFFERS` beside any comparison whose refs disagree, annotating that delta as not apples-to-apples.
Nothing aborts.
An intentional fix changes output by design, as #72 did for `lik()` on mixed data, and only a human knows which side is correct.
A hard failure would need an allowlist entry for every such fix, which is churn that buys nothing a loud line in the report does not.

The digest differs by op class.
Deterministic ops (`forde` parameters, `lik` values) hash the rounded output directly.
Stochastic ops (`forge`, `expct`) run under a fixed seed and hash a rounded summary vector (column means and standard deviations), because an exact hash there also moves when a refactor merely consumes RNG draws in a different order.
The report says which kind of digest produced a flag, so a stochastic mismatch is read as "look at this" rather than "this is broken".

## Layout

```
bench/
  arf-bench         the CLI: run / collect / plan subcommands   (new, Rapp)
  README.md         how to run, memory metric, topology        (exists, extend)
  DESIGN.md         this file                                   (new)
  bench-helpers.R   sampler, data generation, git stamping      (exists, keep)
  cells.R           the grid as data: quick and full tiers      (new)
  refs.R            git ref to installed libpath                (new)
  run-cell.R        one cell against all refs, callr children   (new)
  registry.R        batchtools registry, local or slurm         (new)
  collate.R         reduceResultsDataTable to csv and report.md (new)
  viz.qmd           plots, plus the release trend panel         (exists, extend)
  anchors.csv       curated ref set: git tags and CRAN releases  (new)
  history.csv       committed cache of anchor rows               (new)
  lib/              per-ref and shared libraries, gitignored     (new)
  results/ logs/    raw runs, gitignored                        (exists)
```

`bench` is already in `.Rbuildignore`, so none of this reaches the tarball.
`docs/` is gitignored pkgdown output, which is why this file lives here instead.

## Schema

One row per `(ref, op, cell)`:

```
ref arf_version commit op backend workers n p trees iters
dt_threads ranger_threads metric peak_mb floor_mb peak_delta_mb mem_reps
time_median time_min time_max
digest digest_kind
host kernel r_version job_id timestamp
```

`peak_delta_mb` is the op's own marginal and the only memory column deltas are taken on; `peak_mb` and `floor_mb` are kept so the dilution stays visible.
It is the median of the per-round paired differences, so it is close to but not identical to the difference of the two reported medians.
`mem_reps` records how many rounds produced the row.

`bench_git_commit()` already produces the commit hash, and its existing dirty detection becomes the precondition check rather than a column.
`arf_version` comes from the installed ref's `DESCRIPTION`, which is what makes an anchor row self-describing.
`metric` records which memory measurement was used, and rows with different metrics must never be compared, as `README.md` already warns.
`history.csv` uses these columns unchanged.

## Ref resolution

Refs come in two kinds, both named in `anchors.csv` and both resolving to an installed library.

`git:<ref>` becomes a `git worktree`, installed with `R CMD INSTALL -l`.
This works for unpushed branches, needs no network, and guarantees the installed tree is exactly that ref.

`cran:<version>` installs the published tarball from the CRAN archive, cached under `bench/lib/`.
A git tag and the tarball CRAN actually shipped are not always identical, and for a release anchor the tarball is the more faithful answer to "what did users have".
Use `git:` for anchors inside the current development line and `cran:` for historical releases.

Dependencies are deliberately *not* per ref.
One shared library holds `data.table`, `ranger`, `mirai` and the rest, and each per-ref library holds only `arf`, with `R_LIBS` set to `lib/<ref>:lib/shared`.
Installing dependencies per ref would let versions diverge between refs and confound the comparison with somebody else's performance change.

## Tiers

quick: sequential only, `n` in {1e3, 1e4}, `p = 10`, `trees` in {10, 50}, minutes, runs on bertha or a laptop, and is what a routine PR claim cites.
`make bench` runs it, so it is cheap enough to fire off occasionally mid-work as a sanity check and catch a major regression the day it lands rather than at release.
Small and sequential by construction, so unlike the full tier it cannot take a machine down.
The five ops are `adversarial_rf`, `forde`, `lik`, `forge` and `expct`, matching the configs `submit-ops.sh` already defines.
`impute` is out of the grid until there is a reason to track it.

full: adds `backend` in {sequential, psock, foreach, mirai}, `n = 5e4`, `trees = 200`, `workers` in {1, 2, 4, 8, 16}, and the large `forge`/`expct` row-mode variants that `submit-ops.sh` already defines.
Hours, slurm, for releases and for any PR touching parallel code.

The split is deliberate: bertha for quick, the cluster for full.
Memory on bertha cannot be capped the way a slurm job's cgroup caps it, so a large-data run there risks taking the machine down with it, while the cluster gives each job its own cgroup and its own limit.
That is an argument about the full tier's data sizes, not about bertha, which is why quick stays welcome anywhere.

This also mostly dissolves the metric-comparability worry.
The quick tier is sequential and single-process, with no fork workers, no daemons and no mori regions, so there is no shared memory to account for and even plain RSS would answer correctly.
The `cgroup-anon+shmem` metric exists to compare copy-per-worker against shared-memory backends, which only the full tier does, and that tier runs under slurm where a per-job cgroup is guaranteed.
cgroup v2 is in fact available on bertha, so quick runs there get the preferred metric anyway; a laptop that lacks it still answers correctly for the same single-process reason.
Numbers still carry their `metric` column and are never compared across metrics; a delta between refs measured on one machine with one metric stays valid.

## Report

`report.md` is a delta table per op and cell, ready to paste into a PR, with the 3% and 10% thresholds applied and any `OUTPUT DIFFERS` lines directly beneath the affected rows.
`viz.qmd` gains a panel plotting `history.csv` across releases for runtime and peak memory.

## Non-goals

CI integration.
GitHub runners are noisy, shared, and give no cgroup isolation, so a memory number from one would be indefensible.
The suite runs on bertha or the cluster, by hand or by `sbatch`.

Replacing the backend comparison.
The existing sweep answers "which backend, at what size", which is a different question from "did this change make things worse".

## Anchor compatibility

With `cran:0.2.5` as the only anchor the problem is latent rather than live: 0.2.5 shares essentially its whole API with the development line, and the `mtry` default change that version introduced is shared by both, so no current comparison straddles a behaviour gap.
The mechanism still goes in, because it is what makes adding an older anchor later a safe decision instead of a silent one.

Before running, each ref's `names(formals())` is checked against the arguments the cells actually pass.
A missing argument aborts that cell with a message naming the ref.
This matters because R would otherwise be unhelpfully forgiving: an unmatched argument can be dropped or partially matched to something else, the call succeeds against that version's defaults, and the result is a plausible number measuring a different computation.
That is the one failure this suite must not have.

## Settled

Memory accounting: cgroup v2 is available on bertha, confirmed, so both tiers use the preferred metric.
Nothing depends on it in the quick tier regardless, per the tier split above, so a laptop without it is still usable.
To check a machine, read-only:

```sh
stat -fc %T /sys/fs/cgroup    # cgroup2fs = v2 unified, tmpfs = v1 or hybrid
cg=$(awk -F: '$1=="0"{print $3}' /proc/self/cgroup)
grep -E '^(anon|shmem) ' "/sys/fs/cgroup$cg/memory.stat"
```

Two printed lines mean the preferred metric works there.

Slurm template: the BIPS cluster template already in use, in `bips/bips-cluster`, wired into `makeClusterFunctionsSlurm()`.
The BIPS cluster configures batchtools globally in `/etc/xdg/batchtools/config.R`, which supplies `cluster.functions` with the site template plus `default.resources` (qos, clusters, partition) and `max.concurrent.jobs`.
So the slurm path reads that config and names no template; `ARF_BENCH_SLURM_TMPL` only overrides it.
The local path keeps `conf.file = NA` for determinism, which is the right choice there and the wrong one on the cluster: it would discard every site default.

Three units in that template are easy to get wrong, and all three were.
`ncpus` becomes `--cpus-per-task` and counts **hyperthreads**, since the site config notes one physical core is two threads, so a 16-worker cell needs 34 rather than 17 or it is oversubscribed, which the harness README warns degrades mirai catastrophically.
`memory` is **total** megabytes via `--mem`, and it is mutually exclusive with `mem_per_cpu`, which the template enforces with a `stop()`; since the site defaults already set `memory`, passing `mem_per_cpu` collides with them.
Walltime must fit the QoS: the default `medium` allows 1440 minutes, so the 8h default is inside it, and `short` at 60 minutes is not.
`ARF_BENCH_SLURM_CPUS`, `ARF_BENCH_SLURM_MEM`, `ARF_BENCH_SLURM_WALLTIME` and `ARF_BENCH_SLURM_PARTITION` override, and the submission prints the per-job totals it is asking for.

Cost, for the record: the full tier is 352 cells, and each cell runs `iters x mem_reps` = 25 executions of its op per ref, so 5280 cell children plus as many floor children across three refs.
`ARF_BENCH_OPS` and `ARF_BENCH_MAX_CELLS` cut that down for staging, and memory wants per-op sizing before the whole grid runs: `forde`, `lik` and `adversarial_rf` are modest, while `expct` and the large `forge` variants are what the 256 GB default exists for.

`history.csv`: appended from a run's anchor rows and committed by hand afterwards, for the runs worth keeping.

`make bench` runs the quick tier, alongside the existing `make test` and `make check` targets.

Anchors: `cran:0.2.5` only, for now.
Older releases are not worth anchoring.
`stepsize` is integral to how memory behaves and how `forge()` parallelizes, and `finite_bounds` is integral to the correctness checks and to downstream conditional sampling, so a comparison against a version where those had only just landed measures a different package rather than a slower one.
The ref set is therefore `HEAD`, `main`, `cran:0.2.5`.
