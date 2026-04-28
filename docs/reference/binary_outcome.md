# Example Binary Outcome Dataset

A simulated 12-study meta-analysis of a binary outcome (postoperative
nausea and vomiting) used to demonstrate bayesma workflows for binomial
likelihoods.

## Usage

``` r
binary_outcome
```

## Format

A tibble with 12 rows and 18 columns:

- `Author`:

  Character. Study author/identifier.

- `Year`:

  Numeric. Publication year.

- `Subgroup`:

  Character. Surgical specialty subgroup.

- `Control`:

  Character. Control intervention label.

- `Intervention`:

  Character. Active intervention label.

- `N_Total`:

  Numeric. Total sample size.

- `N_Control`:

  Numeric. Control arm sample size.

- `N_Intervention`:

  Numeric. Intervention arm sample size.

- `Event_Control`:

  Numeric. Events in control arm.

- `Outcome_Control_No`:

  Numeric. Non-events in control arm.

- `Event_Intervention`:

  Numeric. Events in intervention arm.

- `Outcome_Intervention_No`:

  Numeric. Non-events in intervention arm.

- `D1`, `D2`, `D3`, `D4`, `D5`:

  Character. Risk-of-bias domain judgements: `"Low"`, `"Some concerns"`,
  or `"High"`.

- `Overall`:

  Character. Overall risk-of-bias judgement.

## Source

Simulated.
