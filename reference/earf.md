# Shortcut expectation function

Calls `adversarial_rf`, `forde` and `expct`. For repeated application,
it is faster to save outputs of `adversarial_rf` and `forde` and pass
them via `...` or directly use `expct`.

## Usage

``` r
earf(x, ...)
```

## Arguments

- x:

  Input data. Integer variables are recoded as ordered factors with a
  warning. See Details.

- ...:

  Extra parameters to be passed to `adversarial_rf`, `forde` and
  `expct`.

## Value

A one row data frame with values for all query variables.

## References

Watson, D., Blesch, K., Kapar, J., & Wright, M. (2023). Adversarial
random forests for density estimation and generative modeling. In
*Proceedings of the 26th International Conference on Artificial
Intelligence and Statistics*, pp. 5357-5375.

## See also

[`arf`](https://bips-hb.github.io/arf/reference/arf-package.md),
[`adversarial_rf`](https://bips-hb.github.io/arf/reference/adversarial_rf.md),
[`forde`](https://bips-hb.github.io/arf/reference/forde.md),
[`expct`](https://bips-hb.github.io/arf/reference/expct.md)

## Examples

``` r
# What is the expected values of each feature?
earf(iris)
#>   Sepal.Length Sepal.Width Petal.Length Petal.Width   Species
#> 1     5.843333    3.057333        3.758    1.199333 virginica

#' # What is the expected values of Sepal.Length?
earf(iris, query = "Sepal.Length")
#>   Sepal.Length
#> 1     5.843333

# What if we condition on Species = "setosa"?
earf(iris, query = "Sepal.Length", evidence = data.frame(Species = "setosa"))
#>   Sepal.Length
#> 1     5.010467

```
