# Changelog

## arf 0.2.5

- Export sample_from_leaves() for intra-leaf marginal sampling

## arf 0.2.4

CRAN release: 2025-02-24

- Let verbose=FALSE silence (some) warnings

## arf 0.2.3

- Add impute() function for direct missing data imputation with ARF
- Add one-line functions darf(), earf(), rarf()

## arf 0.2.2

- Faster and vectorized conditional sampling
- Use min.bucket argument from ranger to avoid pruning if possible
- Option to sample NAs in generated data if original data contains NAs
- Stepsize in forge() to reduce memory usage
- Option for local and global finite bounds

## arf 0.2.0

CRAN release: 2024-01-24

- Vectorized adversarial resampling
- Speed boost for compiling into a probabilistic circuit
- Conditional densities and sampling
- Bayesian solution for invariant continuous data within leaf nodes
- New function for computing (conditional) expectations
- Options for missing data

## arf 0.1.3

CRAN release: 2023-02-06

- Speed boost for the adversarial resampling step
- Early stopping option for adversarial training
- alpha parameter for regularizing multinomial distributions in forde
- Unified treatment of colnames with internal semantics (y, obs, tree,
  leaf)
