# arf (development version)
* Add mirai/mori parallel backend as an alternative to foreach/doParallel (#62)
  * Shares large read-only inputs (training data, forest, learned parameters) across workers via mori, lowering memory use in `adversarial_rf()`, `forde()`, `forge()`, `expct()`, and `lik()`
  * Enable with active mirai daemons or `options(arf.backend)`, see `?arf-options`
  * Parallel calls now report the backend in use once per session, including when `parallel = TRUE` finds no backend and runs sequentially; silence with `options(arf.verbose = FALSE)`
* Reduce peak memory of `expct()` by processing conditions in bounded blocks (16x lower on one large internal case, at 12-27% more run time; tune via `options(arf.block_rows)`, see `?arf-options`)
* Reduce memory and dispatch overhead in `forde()` (fused per-tree parameter pass, per-tree coverage) and `lik()` (per-batch reduction)
* Fix `forge()` and `lik()` failing under non-forking `foreach` adapters (`doParallel` PSOCK clusters, `doFuture` multisession) with "object not found" errors: parallel worker bodies now take all inputs as explicit arguments (#62)
* Fix `darf()`, `rarf()` and `earf()` ignoring a pre-trained `arf` or `params` passed via `...`: they retrained and recomputed anyway, then failed with "formal argument 'params' matched by multiple actual arguments"
* Fix a crash when arf runs with `parallel = TRUE` inside a forked worker (`future::plan(multicore)`, `mclapply()`, fork-based `doParallel`) while mirai daemons exist in the parent process: the child inherited the daemon connection and used it, which aborts the process since mirai is not fork-safe. Forked children now ignore the parent's daemons and use the `foreach` path; forcing `options(arf.backend = "mirai")` there errors with an explanation
* Fix `lik()` overestimating likelihoods for mixed (continuous and categorical) queries: leaves where one variable block had zero density contributed the other block's density instead of zero, so values were too high and batch-dependent, and the `arf` path disagreed with the slower path (#72)
  * Fix `lik()` assigning a repeated factor pattern the likelihood of the preceding row rather than of its match
  * Fix `lik()` erroring on purely categorical queries with duplicate rows and no `arf` ("column name 'obs' is not found")

# arf 0.2.5

* **Behavior change**: New `mtry` argument for `adversarial_rf()` with default `max(2, floor(sqrt(p)))` instead of ranger's `floor(sqrt(p))`, which gave `mtry = 1` for fewer than 4 features (#59)
* Export `sample_from_leaves()` for intra-leaf marginal sampling
* Avoid fractional recycling of factor column indices in `sample_from_leaves()` (#63)
* Fix `nomatch = "force"` fallback in `forge()` and `expct()` with `evidence_row_mode = "separate"` when evidence rows match no leaf (requires `finite_bounds != "no"` in `forde()`): errored with `data.table` input, and `expct()` silently filled impossible rows with misaligned values and returned `NA`-padded evidence columns for valid rows (#67)

# arf 0.2.4
* Let verbose=FALSE silence (some) warnings

# arf 0.2.3
* Add impute() function for direct missing data imputation with ARF
* Add one-line functions darf(), earf(), rarf()

# arf 0.2.2
* Faster and vectorized conditional sampling
* Use min.bucket argument from ranger to avoid pruning if possible
* Option to sample NAs in generated data if original data contains NAs
* Stepsize in forge() to reduce memory usage
* Option for local and global finite bounds

# arf 0.2.0
* Vectorized adversarial resampling
* Speed boost for compiling into a probabilistic circuit
* Conditional densities and sampling
* Bayesian solution for invariant continuous data within leaf nodes
* New function for computing (conditional) expectations
* Options for missing data

# arf 0.1.3
* Speed boost for the adversarial resampling step 
* Early stopping option for adversarial training
* alpha parameter for regularizing multinomial distributions in forde
* Unified treatment of colnames with internal semantics (y, obs, tree, leaf)
