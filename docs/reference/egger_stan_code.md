# Generate Stan code for an Egger's test specification

Generate Stan code for an Egger's test specification

## Usage

``` r
egger_stan_code(spec, format = TRUE)
```

## Arguments

- spec:

  An `egger_spec` object.

- format:

  Logical. If `TRUE` (default), format via `stanc --auto-format`; falls
  back to the raw program if unavailable.

## Value

An object of class `"egger_stan_code"`.
