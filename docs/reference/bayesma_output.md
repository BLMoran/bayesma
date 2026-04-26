# Assemble a `bayesma` object from pipeline outputs

Combines the outputs of the earlier stages into the standard `bayesma`
return object consumed by downstream tools
([`forest()`](https://blmoran.github.io/bayesma/reference/forest.md),
[`diagnostics()`](https://blmoran.github.io/bayesma/reference/diagnostics.md),
etc.).

## Usage

``` r
bayesma_output(spec, fit, effects)
```

## Arguments

- spec:

  A `bayesma_spec` object.

- fit:

  A `bayesma_fit` object.

- effects:

  A `bayesma_effects` object (from
  [`bayesma_extract()`](https://blmoran.github.io/bayesma/reference/bayesma_extract.md)).

## Value

A list of class `"bayesma"`.
