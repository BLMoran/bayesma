# Print a summary for a fitted model represented by a `bayesma` object

Print a summary for a fitted model represented by a `bayesma` object

## Usage

``` r
# S3 method for class 'bayesma'
print(x, digits = 2, ...)
```

## Arguments

- x:

  An object of class `bayesma`.

- digits:

  The number of significant digits for printing out the summary;
  defaults to 2. Bulk_ESS and Tail_ESS are always rounded to integers.

- ...:

  Additional arguments passed to
  [`summary.bayesma`](https://blmoran.github.io/bayesma/reference/summary.bayesma.md).

## See also

[`summary.bayesma`](https://blmoran.github.io/bayesma/reference/summary.bayesma.md)
