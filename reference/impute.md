# Missing value imputation with ARF

Perform single or multiple imputation with ARFs. Calls `adversarial_rf`,
`forde` and `expct`/`forge`.

## Usage

``` r
impute(
  x,
  m = 1,
  expectation = ifelse(m == 1, TRUE, FALSE),
  num_trees = 100L,
  min_node_size = 10L,
  round = TRUE,
  finite_bounds = "local",
  epsilon = 1e-14,
  verbose = FALSE,
  ...
)
```

## Arguments

- x:

  Input data.

- m:

  Number of imputed datasets to generate. The default is single
  imputation (`m = 1`).

- expectation:

  Return expected value instead of multiple imputations. By default, for
  single imputation (`m = 1`), the expected value is returned.

- num_trees:

  Number of trees to grow in the ARF.

- min_node_size:

  Minimal number of real data samples in leaf nodes.

- round:

  Round continuous variables to their respective maximum precision in
  the real data set?

- finite_bounds:

  Impose finite bounds on all continuous variables? See
  [`forde`](https://bips-hb.github.io/arf/reference/forde.md).

- epsilon:

  Slack parameter on empirical bounds; see
  [`forde`](https://bips-hb.github.io/arf/reference/forde.md).

- verbose:

  Print progress for `adversarial_rf`?

- ...:

  Extra parameters to be passed to `adversarial_rf`, `forde` and
  `expct`/`forge`.

## Value

Imputed data. A single dataset is returned for `m = 1`, a list of
datasets for `m > 1`.

## See also

[`arf`](https://bips-hb.github.io/arf/reference/arf-package.md),
[`forde`](https://bips-hb.github.io/arf/reference/forde.md),
[`forge`](https://bips-hb.github.io/arf/reference/forge.md),
[`expct`](https://bips-hb.github.io/arf/reference/expct.md),
[`lik`](https://bips-hb.github.io/arf/reference/lik.md)

## Examples

``` r
# Generate some missings
iris_na <- iris
for (j in 1:ncol(iris)) {
  iris_na[sample(1:nrow(iris), 5), j] <- NA
}

# Single imputation
iris_imputed <- arf::impute(iris_na, num_trees = 10, m = 1)

# Multiple imputation
iris_imputed <- arf::impute(iris_na, num_trees = 10, m = 10)

if (FALSE) { # \dontrun{
# Parallelization with doParallel
doParallel::registerDoParallel(cores = 4)

# ... or with doFuture
doFuture::registerDoFuture()
future::plan("multisession", workers = 4)
} # }
```
