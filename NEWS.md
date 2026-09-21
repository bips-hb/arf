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
