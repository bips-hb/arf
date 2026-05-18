# Shortcut likelihood function

Calls `adversarial_rf`, `forde` and `lik`. For repeated application, it
is faster to save outputs of `adversarial_rf` and `forde` and pass them
via `...` or directly use `lik`.

## Usage

``` r
darf(x, query = NULL, ...)
```

## Arguments

- x:

  Input data. Integer variables are recoded as ordered factors with a
  warning. See Details.

- query:

  Data frame of samples, optionally comprising just a subset of training
  features. See Details of `lik`. Is set to `x` if `zero`.

- ...:

  Extra parameters to be passed to `adversarial_rf`, `forde` and `lik`.

## Value

A vector of likelihoods, optionally on the log scale. A dataset of
`n_synth` synthetic samples or of `nrow(x)` synthetic samples if
`n_synth` is undefined.

## References

Watson, D., Blesch, K., Kapar, J., & Wright, M. (2023). Adversarial
random forests for density estimation and generative modeling. In
*Proceedings of the 26th International Conference on Artificial
Intelligence and Statistics*, pp. 5357-5375.

## See also

[`arf`](https://bips-hb.github.io/arf/reference/arf-package.md),
[`adversarial_rf`](https://bips-hb.github.io/arf/reference/adversarial_rf.md),
[`forde`](https://bips-hb.github.io/arf/reference/forde.md),
[`forge`](https://bips-hb.github.io/arf/reference/forge.md)

## Examples

``` r
# Estimate log-likelihoods
ll <- darf(iris)

# Partial evidence query
ll <- darf(iris, query = iris[1, 1:3])

# Condition on Species = "setosa"
ll <- darf(iris, query = iris[1, 1:3], evidence = data.frame(Species = "setosa"))

```
