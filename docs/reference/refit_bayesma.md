# Refit a bayesma model on new (subset) data

Uses the stored call_args from the original fit to re-run bayesma() with
different data.

## Usage

``` r
refit_bayesma(model, newdata)
```

## Arguments

- newdata:

  A data frame (subset of original)

- x:

  A bayesma object

## Value

A new bayesma object
