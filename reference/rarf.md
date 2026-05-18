# Shortcut sampling function

Calls `adversarial_rf`, `forde` and `forge`. For repeated application,
it is faster to save outputs of `adversarial_rf` and `forde` and pass
them via `...` or directly use `forge`.

## Usage

``` r
rarf(x, n_synth = NULL, ...)
```

## Arguments

- x:

  Input data. Integer variables are recoded as ordered factors with a
  warning. See Details.

- n_synth:

  Number of synthetic samples to generate for unconditional generation
  with no `evidence` given. Number of synthetic samples to generate per
  `evidence` row if `evidence` is provided. If `NULL`, defaults to
  `nrow(x)` if no `evidence` is provided and to `1` otherwise.

- ...:

  Extra parameters to be passed to `adversarial_rf`, `forde` and
  `forge`.

## Value

A dataset of `n_synth` synthetic samples or of `nrow(x)` synthetic
samples if `n_synth` is undefined.

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
# Generate 150 (size of original iris dataset) synthetic samples from the iris dataset
x_synth <- rarf(iris)

# Generate 100 synthetic samples from the iris dataset
x_synth <- rarf(iris, n_synth = 100)

# Condition on Species = "setosa"
x_synth <- rarf(iris, evidence = data.frame(Species = "setosa"))
```
