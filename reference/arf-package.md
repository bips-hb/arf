# arf: Adversarial Random Forests

Adversarial random forests (ARFs) recursively partition data into fully
factorized leaves, where features are jointly independent. The procedure
is iterative, with alternating rounds of generation and discrimination.
Data becomes increasingly realistic at each round, until original and
synthetic samples can no longer be reliably distinguished. This is
useful for several unsupervised learning tasks, such as density
estimation and data synthesis. Methods for both are implemented in this
package. ARFs naturally handle unstructured data with mixed continuous
and categorical covariates. They inherit many of the benefits of random
forests, including speed, flexibility, and solid performance with
default parameters. For details, see Watson et al. (2023)
<https://proceedings.mlr.press/v206/watson23a.html>.

## See also

[`adversarial_rf`](https://bips-hb.github.io/arf/reference/adversarial_rf.md),
[`forde`](https://bips-hb.github.io/arf/reference/forde.md),
[`forge`](https://bips-hb.github.io/arf/reference/forge.md),
[`expct`](https://bips-hb.github.io/arf/reference/expct.md),
[`lik`](https://bips-hb.github.io/arf/reference/lik.md)

Useful links:

- <https://github.com/bips-hb/arf>

- <https://bips-hb.github.io/arf/>

- Report bugs at <https://github.com/bips-hb/arf/issues>

## Author

**Maintainer**: Marvin N. Wright <cran@wrig.de>
([ORCID](https://orcid.org/0000-0002-8542-6291))

Authors:

- David S. Watson <david.s.watson11@gmail.com>
  ([ORCID](https://orcid.org/0000-0001-9632-2159))

- Kristin Blesch ([ORCID](https://orcid.org/0000-0001-6241-3079))

- Jan Kapar ([ORCID](https://orcid.org/0009-0000-6408-2840))

## Examples

``` r
# Train ARF and estimate leaf parameters
arf <- adversarial_rf(iris)
#> Iteration: 0, Accuracy: 78.19%
#> Iteration: 1, Accuracy: 42.47%
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
#> [1] -0.4604258

# Expectation of Sepal.Length for class setosa
evi <- data.frame(Species = "setosa")
expct(psi, query = "Sepal.Length", evidence = evi)
#>   Sepal.Length
#> 1     5.022068

if (FALSE) { # \dontrun{
# Parallelization with doParallel
doParallel::registerDoParallel(cores = 4)

# ... or with doFuture
doFuture::registerDoFuture()
future::plan("multisession", workers = 4)
} # }
```
