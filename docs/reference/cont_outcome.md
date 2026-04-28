# Example Continuous Outcome Dataset

A simulated 12-study meta-analysis of a continuous outcome used to
demonstrate bayesma workflows for Gaussian likelihoods.

## Usage

``` r
cont_outcome
```

## Format

A tibble with 12 rows and 15 columns:

- `Author`:

  Character. Study author/identifier.

- `Year`:

  Numeric. Publication year.

- `Subgroup`:

  Character. Subgroup label.

- `Mean_Control`:

  Numeric. Control arm mean.

- `SD_Control`:

  Numeric. Control arm standard deviation.

- `N_Control`:

  Numeric. Control arm sample size.

- `Mean_Intervention`:

  Numeric. Intervention arm mean.

- `SD_Intervention`:

  Numeric. Intervention arm standard deviation.

- `N_Intervention`:

  Numeric. Intervention arm sample size.

- `D1`, `D2`, `D3`, `D4`, `D5`:

  Character. Risk-of-bias domain judgements: `"Low"`, `"Some concerns"`,
  or `"High"`.

- `Overall`:

  Character. Overall risk-of-bias judgement.

## Source

Simulated.
