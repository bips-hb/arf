# Post-process data

This function prepares output data for forge.

## Usage

``` r
post_x(x, params, round = TRUE)
```

## Arguments

- x:

  Input data.frame.

- params:

  Circuit parameters learned via
  [`forde`](https://bips-hb.github.io/arf/reference/forde.md).

- round:

  Round continuous variables to their respective maximum precision in
  the real data set?
