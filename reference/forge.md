# Forests for Generative Modeling

Uses pre-trained FORDE model to simulate synthetic data.

## Usage

``` r
forge(
  params,
  n_synth,
  evidence = NULL,
  evidence_row_mode = c("separate", "or"),
  round = TRUE,
  sample_NAs = FALSE,
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

- n_synth:

  Number of synthetic samples to generate.

- evidence:

  Optional set of conditioning events. This can take one of three
  forms: (1) a partial sample, i.e. a single row of data with some but
  not all columns; (2) a data frame of conditioning events, which allows
  for inequalities; or (3) a posterior distribution over leaves. See
  Details.

- evidence_row_mode:

  Interpretation of rows in multi-row evidence. If `"separate"`, each
  row in `evidence` is a unique conditioning event for which `n_synth`
  synthetic samples are generated. If `"or"`, the rows are combined with
  a logical OR. See Examples.

- round:

  Round continuous variables to their respective maximum precision in
  the real data set?

- sample_NAs:

  Sample `NA`s respecting the probability for missing values in the
  original data?

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
  `doParallel` or `doFuture`; see examples.

## Value

A dataset of `n_synth` synthetic samples.

## Details

`forge` simulates a synthetic dataset of `n_synth` samples. First,
leaves are sampled in proportion to either their coverage (if
`evidence = NULL`) or their posterior probability. Then, each feature is
sampled independently within each leaf according to the probability mass
or density function learned by
[`forde`](https://bips-hb.github.io/arf/reference/forde.md). This will
create realistic data so long as the adversarial RF used in the previous
step satisfies the local independence criterion. See Watson et al.
(2023).

There are three methods for (optionally) encoding conditioning events
via the `evidence` argument. The first is to provide a partial sample,
where some columns from the training data are missing or set to `NA`.
The second is to provide a data frame with condition events. This
supports inequalities and intervals. Alternatively, users may directly
input a pre-calculated posterior distribution over leaves, with columns
`f_idx` and `wt`. This may be preferable for complex constraints. See
Examples.

## References

Watson, D., Blesch, K., Kapar, J., & Wright, M. (2023). Adversarial
random forests for density estimation and generative modeling. In
*Proceedings of the 26th International Conference on Artificial
Intelligence and Statistics*, pp. 5357-5375.

## See also

[`arf`](https://bips-hb.github.io/arf/reference/arf-package.md),
[`adversarial_rf`](https://bips-hb.github.io/arf/reference/adversarial_rf.md),
[`forde`](https://bips-hb.github.io/arf/reference/forde.md),
[`expct`](https://bips-hb.github.io/arf/reference/expct.md),
[`lik`](https://bips-hb.github.io/arf/reference/lik.md)

## Examples

``` r
# Train ARF and estimate leaf parameters
arf <- adversarial_rf(iris)
#> Iteration: 0, Accuracy: 76.61%
#> Iteration: 1, Accuracy: 38.8%
psi <- forde(arf, iris)

# Generate 100 synthetic samples from the iris dataset
x_synth <- forge(psi, n_synth = 100)

# Condition on Species = "setosa"
evi <- data.frame(Species = "setosa")
x_synth <- forge(psi, n_synth = 100, evidence = evi)

# Condition on Species = "setosa" and Sepal.Length > 6
evi <- data.frame(Species = "setosa",
                  Sepal.Length = "(6, Inf)")
x_synth <- forge(psi, n_synth = 100, evidence = evi)

# Alternative syntax for </> conditions
evi <- data.frame(Sepal.Length = ">6")
x_synth <- forge(psi, n_synth = 100, evidence = evi)

# Negation condition, i.e. all classes except "setosa"
evi <- data.frame(Species = "!setosa")
x_synth <- forge(psi, n_synth = 100, evidence = evi)

# Condition on first two data rows with some missing values
evi <- iris[1:2,]
evi[1, 1] <- NA_real_
evi[1, 5] <- NA_character_
evi[2, 2] <- NA_real_
x_synth <- forge(psi, n_synth = 1, evidence = evi)

# Or just input some distribution on leaves
# (Weights that do not sum to unity are automatically scaled)
n_leaves <- nrow(psi$forest)
evi <- data.frame(f_idx = psi$forest$f_idx, wt = rexp(n_leaves))
x_synth <- forge(psi, n_synth = 100, evidence = evi)

if (FALSE) { # \dontrun{
# Parallelization with doParallel
doParallel::registerDoParallel(cores = 4)

# ... or with doFuture
doFuture::registerDoFuture()
future::plan("multisession", workers = 4)
} # }
```
