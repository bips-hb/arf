# Package index

## Fitting an ARF

Adversarial random forest training, the starting point for everything
else.

- [`arf`](https://bips-hb.github.io/arf/dev/reference/arf-package.md)
  [`arf-package`](https://bips-hb.github.io/arf/dev/reference/arf-package.md)
  : arf: Adversarial Random Forests
- [`adversarial_rf()`](https://bips-hb.github.io/arf/dev/reference/adversarial_rf.md)
  : Adversarial Random Forests

## Density estimation and generation

Learn leaf parameters, then sample, impute or evaluate likelihoods and
expectations, optionally conditioned on evidence.

- [`forde()`](https://bips-hb.github.io/arf/dev/reference/forde.md) :
  Forests for Density Estimation
- [`forge()`](https://bips-hb.github.io/arf/dev/reference/forge.md) :
  Forests for Generative Modeling
- [`lik()`](https://bips-hb.github.io/arf/dev/reference/lik.md) :
  Likelihood Estimation
- [`expct()`](https://bips-hb.github.io/arf/dev/reference/expct.md) :
  Expected Value
- [`impute()`](https://bips-hb.github.io/arf/dev/reference/impute.md) :
  Missing value imputation with ARF
- [`sample_from_leaves()`](https://bips-hb.github.io/arf/dev/reference/sample_from_leaves.md)
  : Generate synthetic data by sampling from the leaves of a random
  forest

## Shortcut functions

One-call wrappers that train, learn parameters and generate in one step.

- [`darf()`](https://bips-hb.github.io/arf/dev/reference/darf.md) :
  Shortcut likelihood function
- [`earf()`](https://bips-hb.github.io/arf/dev/reference/earf.md) :
  Shortcut expectation function
- [`rarf()`](https://bips-hb.github.io/arf/dev/reference/rarf.md) :
  Shortcut sampling function

## Package options

- [`arf-options`](https://bips-hb.github.io/arf/dev/reference/arf-options.md)
  [`arf.backend`](https://bips-hb.github.io/arf/dev/reference/arf-options.md)
  [`arf.verbose`](https://bips-hb.github.io/arf/dev/reference/arf-options.md)
  [`arf.block_rows`](https://bips-hb.github.io/arf/dev/reference/arf-options.md)
  : arf package options
