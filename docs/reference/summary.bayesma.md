# Summarise a fitted `bayesma` model

Collects all relevant posterior summaries into a `bayesma_summary`
object that `print.bayesma_summary` renders in brms-style layout.

## Usage

``` r
# S3 method for class 'bayesma'
summary(object, ...)
```

## Arguments

- object:

  An object of class `bayesma`.

- ...:

  Currently unused.

## Value

An object of class `bayesma_summary` (a named list).
