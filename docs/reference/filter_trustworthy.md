# Filter Studies by INSPECT-SR Trustworthiness

Filters a meta-analysis dataset based on INSPECT-SR results.

## Usage

``` r
filter_trustworthy(
  data,
  inspect_results,
  studyvar = study,
  exclude = c("serious", "some"),
  domain = "Overall"
)
```

## Arguments

- data:

  A data frame with one row per study.

- inspect_results:

  An object returned by
  [`inspect_sr()`](https://blmoran.github.io/bayesma/reference/inspect_sr.md).

- studyvar:

  Unquoted column name for the study identifier in `data` (tidyeval).
  Must match the Study column in `inspect_results`.

- exclude:

  `"serious"` (default) drops studies with "Serious concerns"; `"some"`
  drops both "Some concerns" and "Serious concerns".

- domain:

  Which domain: `"Overall"` (default), `"D1"`, `"D2"`, `"D3"`, or
  `"D4"`.

## Value

The filtered data frame, with attribute `excluded_studies`.

## Examples

``` r
if (FALSE) { # \dontrun{
res <- inspect_sr(inspect_sr_example)
clean <- filter_trustworthy(my_data, res, studyvar = study)
} # }
```
