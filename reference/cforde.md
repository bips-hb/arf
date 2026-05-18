# Compute conditional circuit parameters

Compute conditional circuit parameters

## Usage

``` r
cforde(
  params,
  evidence,
  row_mode = c("separate", "or"),
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

- evidence:

  Data frame of conditioning event(s).

- row_mode:

  Interpretation of rows in multi-row conditions.

- nomatch:

  What to do if no leaf matches a condition in `evidence`? Options are
  to force sampling from a random leaf (`"force"`) or return `NA`
  (`"na"`). The default is `"force"`.

- verbose:

  Show warnings, e.g. when no leaf matches a condition?

- stepsize:

  Stepsize defining number of condition rows handled in one for each
  step.

- parallel:

  Compute in parallel? Must register backend beforehand, e.g. via
  `doParallel` or `doFuture`; see examples.

## Value

List with conditions (`evidence_input`), prepared conditions
(`evidence_prepped`) and leaves that match the conditions in evidence
with continuous data (`cnt`) and categorical data (`cat`) as well as
leaf info (`forest`).
