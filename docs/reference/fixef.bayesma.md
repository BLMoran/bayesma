# Extract fixed effect summary from a bayesma object

Returns a named numeric vector mimicking the structure of brms::fixef(),
with elements: Estimate, Est.Error, Q2.5, Q97.5

## Usage

``` r
fixef.bayesma(x)
```

## Arguments

- x:

  A bayesma object

## Value

A 1x4 matrix with columns Estimate, Est.Error, Q2.5, Q97.5
