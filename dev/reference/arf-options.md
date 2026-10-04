# arf package options

Options controlling the parallel backend and its messaging, set via
[`options`](https://rdrr.io/r/base/options.html).

## Details

- `arf.backend`:

  Parallel backend used when `parallel = TRUE`: `"foreach"` or
  `"mirai"`. If unset, arf uses `"mirai"` when mirai daemons are running
  and `"foreach"` otherwise.

- `arf.verbose`:

  Report the selected backend once per backend configuration per
  session? Default `TRUE`; set `FALSE` to silence.

- `arf.block_rows`:

  Cap on rows materialized per block of conditions in
  [`expct`](https://bips-hb.github.io/arf/dev/reference/expct.md).
  Default `5e6`. Lower it to trade speed for memory on large forests
  with many conditions.

`arf.block_rows` does not affect results, only peak memory and speed.
For memory-constrained hardware, combine a small daemon count with a
lower `arf.block_rows` and the `batch`/`stepsize` arguments of
[`lik`](https://bips-hb.github.io/arf/dev/reference/lik.md),
[`forge`](https://bips-hb.github.io/arf/dev/reference/forge.md) and
[`expct`](https://bips-hb.github.io/arf/dev/reference/expct.md).

The `"foreach"` backend uses whatever adapter is registered (e.g.
`doParallel`, `doFuture`). The `"mirai"` backend uses `mirai` daemons
and shares large read-only inputs (training data, forest, learned
parameters) across workers via `mori`, so workers do not each copy them.
In our benchmarks (n = 20000, 200 trees, up to 16 workers) run time was
similar between the backends, though on a cold daemon pool `"mirai"` can
be noticeably slower. The memory benefit is largest for the
tree-parallel operations
([`forde`](https://bips-hb.github.io/arf/dev/reference/forde.md),
[`adversarial_rf`](https://bips-hb.github.io/arf/dev/reference/adversarial_rf.md))
on large forests with many workers: at 16 workers peak memory was about
30 percent lower for `forde` and about half for `adversarial_rf`. Expect
smaller gains on smaller problems, where the daemon pool's fixed
overhead can outweigh the sharing, so prefer `"mirai"` at scale and
either backend otherwise. Daemons started via `future.mirai` (e.g.
`plan(future.mirai::mirai_multisession)`) are detected like any other
mirai daemons, so futureverse users get the `"mirai"` backend
automatically. None of these packages are installed with arf: install
`mirai` and `mori` for the `"mirai"` backend, `doParallel` or `doFuture`
for `"foreach"`, and `doRNG` or `future.mirai` if you use them below.

Workers and `data.table` threads multiply. If their product exceeds the
core count, cap threads per worker with
[`data.table::setDTthreads()`](https://rdrr.io/pkg/data.table/man/openmp-utils.html)
and set `OMP_WAIT_POLICY=passive` before starting R so idle OpenMP
threads sleep instead of spinning.

Reproducibility of stochastic operations
([`forge`](https://bips-hb.github.io/arf/dev/reference/forge.md),
categorical
[`expct`](https://bips-hb.github.io/arf/dev/reference/expct.md)) under
parallel execution: [`set.seed`](https://rdrr.io/r/base/Random.html)
only governs the calling process, not the workers. With the `"mirai"`
backend, seed the daemons instead: `mirai::daemons(n, seed = 42)` gives
reproducible results provided the daemon count, the seed and the
sequence of calls on a fresh daemon pool are kept fixed (changing the
daemon count changes how work is chunked and therefore the random stream
assignment). This also applies to pools started via `future.mirai`:
`future`'s own seed machinery covers only work dispatched through its
API, which arf's backend bypasses, and `plan()` does not accept a `seed`
argument, so seed the pool with
[`mirai::daemons()`](https://mirai.r-lib.org/reference/daemons.html)
directly. For the `"foreach"` backend, register `doRNG` on top of the
adapter (`doRNG::registerDoRNG(42)` after `registerDoParallel()` or
`registerDoFuture()`): results are then reproducible, identical across
adapters, and independent of the worker count, provided `stepsize` is
set explicitly (its default depends on the worker count). Without
`doRNG`, `doFuture` flags the stochastic operations with "UNRELIABLE
VALUE" warnings. Sequential execution (`parallel = FALSE`) with
`set.seed` is exact as always.

## Examples

``` r
if (FALSE) { # \dontrun{
arf <- adversarial_rf(iris)

# foreach backend
doParallel::registerDoParallel(cores = 4)
psi <- forde(arf, iris)

# mirai backend: start daemons, then call as usual
mirai::daemons(4)
psi <- forde(arf, iris)
mirai::daemons(0)  # shut down when done

# futureverse: future.mirai daemons are detected automatically
future::plan(future.mirai::mirai_multisession, workers = 4)
psi <- forde(arf, iris)
future::plan("sequential")  # shuts the daemons down

# force a backend regardless of what is registered (NULL restores auto-detection)
options(arf.backend = "mirai")
options(arf.backend = NULL)

# silence the backend message
options(arf.verbose = FALSE)

# reproducible parallel sampling with mirai: seeded daemons
# (fixed daemon count, fresh pool)
evi <- data.frame(Species = sample(levels(iris$Species), 100, replace = TRUE))
mirai::daemons(4, seed = 42)
# stepsize = evidence rows per step; 25 = 100 conditions / 4 workers,
# i.e. the default sizing, made explicit
x <- forge(psi, n_synth = 1, evidence = evi, stepsize = 25)
mirai::daemons(0)

# reproducible parallel sampling with foreach: doRNG on top of the
# adapter; set stepsize explicitly (its default depends on worker count)
doParallel::registerDoParallel(cores = 4)
doRNG::registerDoRNG(42)
x <- forge(psi, n_synth = 1, evidence = evi, stepsize = 25)
} # }
```
