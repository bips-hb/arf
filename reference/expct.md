# Expected Value

Compute the expectation of some query variable(s), optionally
conditioned on some event(s).

## Usage

``` r
expct(
  params,
  query = NULL,
  evidence = NULL,
  evidence_row_mode = c("separate", "or"),
  round = FALSE,
  nomatch = c("force", "na"),
  verbose = TRUE,
  stepsize = 0,
  parallel = TRUE
)
```

## Arguments

- params:

  Circuit parameters learned via
  [`forde`](https://bips-hb.github.io/arf/reference/forde.md).

- query:

  Optional character vector of variable names. Estimates will be
  computed for each. If `NULL`, all variables other than those in
  `evidence` will be estimated. If `evidence` contains `NA`s, those
  values will be imputed and a full dataset is returned.

- evidence:

  Optional set of conditioning events. This can take one of three
  forms: (1) a partial sample, i.e. a single row of data with some but
  not all columns; (2) a data frame of conditioning events, which allows
  for inequalities and intervals; or (3) a posterior distribution over
  leaves. See Details and Examples.

- evidence_row_mode:

  Interpretation of rows in multi-row evidence. If `"separate"`, each
  row in `evidence` is a unique conditioning event for which `n_synth`
  synthetic samples are generated. If `"or"`, the rows are combined with
  a logical OR. See Examples.

- round:

  Round continuous variables to their respective maximum precision in
  the real data set?

- nomatch:

  What to do if no leaf matches a condition in `evidence`? Options are
  to force sampling from a random leaf (`"force"`) or return `NA`
  (`"na"`). The default is `"force"`.

- verbose:

  Show warnings, e.g. when no leaf matches a condition?

- stepsize:

  How many rows of evidence should be handled at each step? Defaults to
  `nrow(evidence) / num_registered_workers` for `parallel == TRUE`.

- parallel:

  Compute in parallel? Must register backend beforehand, e.g. via
  `doParallel` or `doFuture`; see Examples.

## Value

A one row data frame with values for all query variables.

## Details

This function computes expected values for any subset of features,
optionally conditioned on some event(s).

There are three methods for (optionally) encoding conditioning events
via the `evidence` argument. The first is to provide a partial sample,
where some columns from the training data are missing or set to `NA`.
The second is to provide a data frame with condition events. This
supports inequalities and intervals. Alternatively, users may directly
input a pre-calculated posterior distribution over leaves, with columns
`f_idx` and `wt`. This may be preferable for complex constraints. See
Examples.

Please note that results for continuous features which are both included
in `query` and in `evidence` with an interval condition are currently
inconsistent.

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
[`lik`](https://bips-hb.github.io/arf/reference/lik.md)

## Examples

``` r
# Train ARF and estimate leaf parameters
arf <- adversarial_rf(iris)
#> Iteration: 0, Accuracy: 77.03%
#> Iteration: 1, Accuracy: 33.22%
psi <- forde(arf, iris)

# What is the expected value of Sepal.Length?
expct(psi, query = "Sepal.Length")
#>   Sepal.Length
#> 1     5.843333

# What if we condition on Species = "setosa"?
evi <- data.frame(Species = "setosa")
expct(psi, query = "Sepal.Length", evidence = evi)
#>   Sepal.Length
#> 1     5.026998

# Compute expectations for all features other than Species
expct(psi, evidence = evi)
#>   Sepal.Length Sepal.Width Petal.Length Petal.Width
#> 1     5.026998    3.414044     1.542446   0.2771403

# Condition on Species = "setosa" and Petal.Width > 0.3
evi <- data.frame(Species = "setosa", 
                  Petal.Width = ">0.3")
expct(psi, evidence = evi)
#>   Sepal.Length Sepal.Width Petal.Length
#> 1     5.208278    3.515179     1.755113

# Condition on first two rows with some missing values
evi <- iris[1:2,]
evi[1, 1] <- NA_real_
evi[1, 5] <- NA_character_
evi[2, 2] <- NA_real_
x_synth <- expct(psi, evidence = evi)

if (FALSE) { # \dontrun{
# Parallelization with doParallel
doParallel::registerDoParallel(cores = 4)

# ... or with doFuture
doFuture::registerDoFuture()
future::plan("multisession", workers = 4)
} # }
```
