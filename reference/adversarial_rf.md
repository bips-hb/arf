# Adversarial Random Forests

Implements an adversarial random forest to learn independence-inducing
splits.

## Usage

``` r
adversarial_rf(
  x,
  num_trees = 10L,
  min_node_size = 2L,
  delta = 0,
  max_iters = 10L,
  early_stop = TRUE,
  prune = TRUE,
  verbose = TRUE,
  parallel = TRUE,
  ...
)
```

## Arguments

- x:

  Input data. Integer variables are recoded as ordered factors with a
  warning. See Details.

- num_trees:

  Number of trees to grow in each forest. The default works well for
  most generative modeling tasks, but should be increased for likelihood
  estimation. See Details.

- min_node_size:

  Minimal number of real data samples in leaf nodes.

- delta:

  Tolerance parameter. Algorithm converges when OOB accuracy is \< 0.5 +
  `delta`.

- max_iters:

  Maximum iterations for the adversarial loop.

- early_stop:

  Terminate loop if performance fails to improve from one round to the
  next?

- prune:

  Impose `min_node_size` by pruning?

- verbose:

  Print discriminator accuracy after each round? Will also show
  additional warnings.

- parallel:

  Compute in parallel? Must register backend beforehand, e.g. via
  `doParallel` or `doFuture`; see examples.

- ...:

  Extra parameters to be passed to `ranger`.

## Value

A random forest object of class `ranger`.

## Details

The adversarial random forest (ARF) algorithm partitions data into fully
factorized leaves where features are jointly independent. ARFs are
trained iteratively, with alternating rounds of generation and
discrimination. In the first instance, synthetic data is generated via
independent bootstraps of each feature, and a RF classifier is trained
to distinguish between real and fake samples. In subsequent rounds,
synthetic data is generated separately in each leaf, using splits from
the previous forest. This creates increasingly realistic data that
satisfies local independence by construction. The algorithm converges
when a RF cannot reliably distinguish between the two classes, i.e. when
OOB accuracy falls below 0.5 + `delta`.

ARFs are useful for several unsupervised learning tasks, such as density
estimation (see
[`forde`](https://bips-hb.github.io/arf/reference/forde.md)) and data
synthesis (see
[`forge`](https://bips-hb.github.io/arf/reference/forge.md)). For the
former, we recommend increasing the number of trees for improved
performance (typically on the order of 100-1000 depending on sample
size).

Integer variables are recoded with a warning (set `verbose = FALSE` to
silence these). Default behavior is to convert integer variables with
six or more unique values to numeric, while those with up to five unique
values are treated as ordered factors. To override this behavior,
explicitly recode integer variables to the target type prior to
training.

Note: convergence is not guaranteed in finite samples. The `max_iters`
argument sets an upper bound on the number of training rounds. Similar
results may be attained by increasing `delta`. Even a single round can
often give good performance, but data with strong or complex
dependencies may require more iterations. With the default
`early_stop = TRUE`, the adversarial loop terminates if performance does
not improve from one round to the next, in which case further training
may be pointless.

## References

Watson, D., Blesch, K., Kapar, J., & Wright, M. (2023). Adversarial
random forests for density estimation and generative modeling. In
*Proceedings of the 26th International Conference on Artificial
Intelligence and Statistics*, pp. 5357-5375.

## See also

[`arf`](https://bips-hb.github.io/arf/reference/arf-package.md),
[`forde`](https://bips-hb.github.io/arf/reference/forde.md),
[`forge`](https://bips-hb.github.io/arf/reference/forge.md),
[`expct`](https://bips-hb.github.io/arf/reference/expct.md),
[`lik`](https://bips-hb.github.io/arf/reference/lik.md)

## Examples

``` r
# Train ARF and estimate leaf parameters
arf <- adversarial_rf(iris)
#> Iteration: 0, Accuracy: 70.71%
#> Iteration: 1, Accuracy: 35.47%
#> Warning: executing %dopar% sequentially: no parallel backend registered
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
#> [1] -0.4833021

# Expectation of Sepal.Length for class setosa
evi <- data.frame(Species = "setosa")
expct(psi, query = "Sepal.Length", evidence = evi)
#>   Sepal.Length
#> 1      5.01668

if (FALSE) { # \dontrun{
# Parallelization with doParallel
doParallel::registerDoParallel(cores = 4)

# ... or with doFuture
doFuture::registerDoFuture()
future::plan("multisession", workers = 4)
} # }
```
