# Example INSPECT-SR Dataset: EEG-Guided Anaesthesia and Delirium

A simulated eight-trial systematic review used to demonstrate the full
INSPECT-SR workflow. The eight trials are hand-crafted to span a range
of trustworthiness profiles and to exercise every automated check.

## Usage

``` r
inspect_sr_example
```

## Format

A tibble with one row per study. Every INSPECT-SR item has its own
column so that the expected layout is immediately obvious.
Variable-length inputs to the automated Domain 4 checks (Table 1
baselines and reported test statistics) are stored in list-columns.

- `study`:

  Character. Study identifier (must be unique).

- `n_randomised_int`, `n_randomised_ctrl`, `n_randomised_total`,
  `n_analysed_int`, `n_analysed_ctrl`, `n_lost_int`, `n_lost_ctrl`:

  Integer participant-flow columns feeding check 4.6.

- `d1_1`..`d1_3`, `d2_1`..`d2_5`, `d3_1`..`d3_2`, `d4_1`..`d4_11`:

  Character. Manual judgement for each INSPECT-SR item. Allowed values:
  `"No concerns"`, `"Some concerns"`, `"Serious concerns"`, `NA`. The
  four automated items (`d4_3`, `d4_6`, `d4_8`, `d4_9`) are left as `NA`
  and filled in by
  [`inspect_sr()`](https://blmoran.github.io/bayesma/reference/inspect_sr.md).

- `d1_1_comment`..`d4_11_comment`:

  Character. Optional free-text comment attached to each judgement.

- `baseline`:

  List-column. Each element is a data frame of Table 1 variables with
  columns `variable`, `mean_int`, `sd_int`, `mean_ctrl`, `sd_ctrl`,
  `p_value`, `integer_scale`. Feeds GRIM (4.8) and Carlisle (4.3).

- `statistics`:

  List-column. Each element is a data frame of reported tests with
  columns `test_type`, `statistic`, `df`, `df2`, `reported_p`,
  `context`. Feeds the p-value check (4.9).

- `outcome_estimate`, `outcome_ci_lower`, `outcome_ci_upper`,
  `outcome_se`, `outcome_log_scale`:

  Primary outcome (odds ratio for delirium incidence) — values are on
  the log scale.

## Source

Simulated.

## Details

The eight trials, their expected behaviour under the automated checks,
and the INSPECT-SR item they illustrate:

|  |  |  |
|----|----|----|
| **Study** | **Profile** | **Expected flag** |
| Callahan (1994) | Clean small trial | None |
| Peek (2005) | Clean large multi-centre | None |
| Clayton (2010) | Clean | None |
| Singer (2003) | Clean (large trial) | None |
| Underwood (2013) | Some concerns (late registration) | Manual D2 |
| Rosa (2008) | Too-perfect baseline balance | Carlisle (4.3) |
| Doe (1995) | Retracted; GRIM failures on integer scale | GRIM (4.8), manual D1/D2/D3 |
| Kint (1995) | CONSORT arithmetic off | N consistency (4.6) |

All values are fabricated for didactic purposes. The layout matches the
schema expected by
[`inspect_sr()`](https://blmoran.github.io/bayesma/reference/inspect_sr.md),
so the whole review can be assessed in one call.

## Examples

``` r
data(inspect_sr_example)

# Frequentist assessment
res <- inspect_sr(inspect_sr_example, verbose = FALSE)
res
#> 
#> INSPECT-SR Trustworthiness Assessment
#> ================================================== 
#> 
#> Study                     D1    D2    D3    D4       Overall           
#> ---------------------------------------------------------------------- 
#> Callahan (1994)           OK    OK    OK    OK       No concerns       
#> Peek (2005)               OK    OK    OK    OK       No concerns       
#> Clayton (2010)            OK    OK    OK    OK       No concerns       
#> Singer (2003)             OK    OK    OK    OK       No concerns       
#> Underwood (2013)          OK    SOME  --    OK       Some concerns     
#> Rosa (2008)               OK    SOME  --    SERIOUS  Serious concerns  
#> Doe (1995)                SERIOUS SERIOUS SOME  SERIOUS  Serious concerns  
#> Kint (1995)               --    --    --    SERIOUS  Serious concerns  
#> 
#> Domains: D1 post-publication, D2 conduct/governance,
#>          D3 text/figures, D4 results (auto-filled for 4.3/4.6/4.8/4.9)
#> OK = No concerns, SOME = Some concerns, SERIOUS = Serious concerns
#> -- = Not assessed
#> 
#> For a per-check table, call inspect_sr_table().
#> 

# Per-check table
inspect_sr_table(res)


  


INSPECT-SR Automated Check Details
```
