# Likelihood Estimation

Compute the likelihood of input data, optionally conditioned on some
event(s).

## Usage

``` r
lik(
  params,
  query,
  evidence = NULL,
  arf = NULL,
  oob = FALSE,
  log = TRUE,
  batch = NULL,
  parallel = TRUE
)
```

## Arguments

- params:

  Circuit parameters learned via
  [`forde`](https://bips-hb.github.io/arf/reference/forde.md).

- query:

  Data frame of samples, optionally comprising just a subset of training
  features. Likelihoods will be computed for each sample. Missing
  features will be marginalized out. See Details.

- evidence:

  Optional set of conditioning events. This can take one of three
  forms: (1) a partial sample, i.e. a single row of data with some but
  not all columns; (2) a data frame of conditioning events, which allows
  for inequalities; or (3) a posterior distribution over leaves. See
  Details.

- arf:

  Pre-trained
  [`adversarial_rf`](https://bips-hb.github.io/arf/reference/adversarial_rf.md)
  or other object of class `ranger`. This is not required but speeds up
  computation considerably for total evidence queries. (Ignored for
  partial evidence queries.)

- oob:

  Only use out-of-bag leaves for likelihood estimation? If `TRUE`, `x`
  must be the same dataset used to train `arf`. Only applicable for
  total evidence queries.

- log:

  Return likelihoods on log scale? Recommended to prevent underflow.

- batch:

  Batch size. The default is to compute densities for all of queries in
  one round, which is always the fastest option if memory allows.
  However, with large samples or many trees, it can be more memory
  efficient to split the data into batches. This has no impact on
  results.

- parallel:

  Compute in parallel? Must register backend beforehand, e.g. via
  `doParallel` or `doFuture`; see examples.

## Value

A vector of likelihoods, optionally on the log scale.

## Details

This function computes the likelihood of input data, optionally
conditioned on some event(s). Queries may be partial, i.e. covering some
but not all features, in which case excluded variables will be
marginalized out.

There are three methods for (optionally) encoding conditioning events
via the `evidence` argument. The first is to provide a partial sample,
where some but not all columns from the training data are present. The
second is to provide a data frame with three columns: `variable`,
`relation`, and `value`. This supports inequalities via `relation`.
Alternatively, users may directly input a pre-calculated posterior
distribution over leaves, with columns `f_idx` and `wt`. This may be
preferable for complex constraints. See Examples.

## References

Watson, D., Blesch, K., Kapar, J., & Wright, M. (2023). Adversarial
random forests for density estimation and generative modeling. In
*Proceedings of the 26th International Conference on Artificial
Intelligence and Statistics*, pp. 5357-5375.

## See also

[`arf`](https://bips-hb.github.io/arf/reference/arf-package.md),
[`adversarial_rf`](https://bips-hb.github.io/arf/reference/adversarial_rf.md),
[`forde`](https://bips-hb.github.io/arf/reference/forde.md),
[`forge`](https://bips-hb.github.io/arf/reference/forge.md),
[`expct`](https://bips-hb.github.io/arf/reference/expct.md)

## Examples

``` r
# Train ARF and estimate leaf parameters
arf <- adversarial_rf(iris)
#> Iteration: 0, Accuracy: 74.16%
#> Iteration: 1, Accuracy: 43.05%
psi <- forde(arf, iris)

# Estimate average log-likelihood
ll <- lik(psi, iris, arf = arf, log = TRUE)
mean(ll)
#> [1] -0.4781229

# Identical but slower
ll <- lik(psi, iris, log = TRUE)
#> Warning: For total evidence queries, it is faster to include the pre-trained arf.
mean(ll)
#> [1] -0.4781229

# Partial evidence query
lik(psi, query = iris[1, 1:3])
#> [1] 0.1577479

# Condition on Species = "setosa"
evi <- data.frame(Species = "setosa")
lik(psi, query = iris[1, 1:3], evidence = evi)
#> [1] 1.255691

# Condition on Species = "setosa" and Petal.Width > 0.3
evi <- data.frame(Species = "setosa", 
                  Petal.Width = ">0.3")
lik(psi, query = iris[1, 1:3], evidence = evi)
#> [1] 1.216418

if (FALSE) { # \dontrun{
# Parallelization with doParallel
doParallel::registerDoParallel(cores = 4)

# ... or with doFuture
doFuture::registerDoFuture()
future::plan("multisession", workers = 4)
} # }
```
