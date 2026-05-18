# Forests for Density Estimation

Uses a pre-trained ARF model to estimate leaf and distribution
parameters.

## Usage

``` r
forde(
  arf,
  x,
  oob = FALSE,
  family = "truncnorm",
  finite_bounds = c("no", "local", "global"),
  alpha = 0,
  epsilon = 0,
  parallel = TRUE
)
```

## Arguments

- arf:

  Pre-trained
  [`adversarial_rf`](https://bips-hb.github.io/arf/reference/adversarial_rf.md).
  Alternatively, any object of class `ranger`.

- x:

  Training data for estimating parameters.

- oob:

  Only use out-of-bag samples for parameter estimation? If `TRUE`, `x`
  must be the same dataset used to train `arf`. Set to `"inbag"` to only
  use in-bag samples. Default is `FALSE`, i.e. use all observations.

- family:

  Distribution to use for density estimation of continuous features.
  Current options include truncated normal (the default
  `family = "truncnorm"`) and uniform (`family = "unif"`). See Details.

- finite_bounds:

  Impose finite bounds on all continuous variables? If `"local"`,
  infinite bounds are set to empirical extrema within leaves. If
  `"global"`, infinite bounds are set to global empirical extrema. if
  `"no"` (the default), infinite bounds are left unchanged.

- alpha:

  Optional pseudocount for Laplace smoothing of categorical features.
  This avoids zero-mass points when test data fall outside the support
  of training data. Effectively parameterizes a flat Dirichlet prior on
  multinomial likelihoods.

- epsilon:

  Optional slack parameter on empirical bounds when
  `finite_bounds != "no"`. This avoids zero-density points when test
  data fall outside the support of training data. The gap between lower
  and upper bounds is expanded by a factor of `1 + epsilon`.

- parallel:

  Compute in parallel? Must register backend beforehand, e.g. via
  `doParallel` or `doFuture`; see examples.

## Value

A `list` with 5 elements: (1) parameters for continuous data; (2)
parameters for discrete data; (3) leaf indices and coverage; (4)
metadata on variables; and (5) the data input class. This list is used
for estimating likelihoods with
[`lik`](https://bips-hb.github.io/arf/reference/lik.md) and generating
data with [`forge`](https://bips-hb.github.io/arf/reference/forge.md).

## Details

`forde` extracts leaf parameters from a pretrained forest and learns
distribution parameters for data within each leaf. The former includes
coverage (proportion of data falling into the leaf) and split criteria.
The latter includes proportions for categorical features and
mean/variance for continuous features. The result is a probabilistic
circuit, stored as a `data.table`, which can be used for various
downstream inference tasks.

Currently, `forde` only provides support for a limited number of
distributional families: truncated normal or uniform for continuous
data, and multinomial for discrete data.

Though `forde` was designed to take an adversarial random forest as
input, the function's first argument can in principle be any object of
class `ranger`. This allows users to test performance with alternative
pipelines (e.g., with supervised forest input). There is also no
requirement that `x` be the data used to fit `arf`, unless `oob = TRUE`.
In fact, using another dataset here may protect against overfitting.
This connects with Wager & Athey's (2018) notion of "honest trees".

## References

Watson, D., Blesch, K., Kapar, J., & Wright, M. (2023). Adversarial
random forests for density estimation and generative modeling. In
*Proceedings of the 26th International Conference on Artificial
Intelligence and Statistics*, pp. 5357-5375.

Wager, S. & Athey, S. (2018). Estimation and inference of heterogeneous
treatment effects using random forests. *J. Am. Stat. Assoc.*,
*113*(523): 1228-1242.

## See also

[`arf`](https://bips-hb.github.io/arf/reference/arf-package.md),
[`adversarial_rf`](https://bips-hb.github.io/arf/reference/adversarial_rf.md),
[`forge`](https://bips-hb.github.io/arf/reference/forge.md),
[`expct`](https://bips-hb.github.io/arf/reference/expct.md),
[`lik`](https://bips-hb.github.io/arf/reference/lik.md)

## Examples

``` r
# Train ARF and estimate leaf parameters
arf <- adversarial_rf(iris)
#> Iteration: 0, Accuracy: 78.52%
#> Iteration: 1, Accuracy: 41.16%
psi <- forde(arf, iris)

# Generate 100 synthetic samples from the iris dataset
x_synth <- forge(psi, n_synth = 100)

# Condition on Species = "setosa" and Sepal.Length > 6
evi <- data.frame(Species = "setosa",
                  Sepal.Length = "(6, Inf)")
x_synth <- forge(psi, n_synth = 100, evidence = evi)

# Estimate average log-likelihood
ll <- lik(psi, iris, arf = arf, log = TRUE)
mean(ll)
#> [1] -0.4118479

# Expectation of Sepal.Length for class setosa
evi <- data.frame(Species = "setosa")
expct(psi, query = "Sepal.Length", evidence = evi)
#>   Sepal.Length
#> 1      5.01768

if (FALSE) { # \dontrun{
# Parallelization with doParallel
doParallel::registerDoParallel(cores = 4)

# ... or with doFuture
doFuture::registerDoFuture()
future::plan("multisession", workers = 4)
} # }

```
