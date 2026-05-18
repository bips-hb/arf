# Adaptive column renaming

This function renames columns in case the input colnames includes any
colnames required by internal functions (e.g., `"y"`).

## Usage

``` r
col_rename(cn, old_name)
```

## Arguments

- cn:

  Column names.

- old_name:

  Name of column to be renamed.
