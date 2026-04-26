# Generating a Bayesian meta-analysis report

## Overview

[`bayesma_report()`](https://blmoran.github.io/bayesma/reference/bayesma_report.md)
renders a parameterised Quarto template into a publication-style
statistical report. The template structure follows the [workflow
guide](https://blmoran.github.io/bayesma/articles/workflow-guide.qmd)
(data prep → fitting → comparison → bias → sensitivity → interpretation)
and the prose style of a full systematic review (Introduction → Methods
→ Results → Overall Interpretation).

You supply the data, a few measure-specific column names, and a handful
of narrative strings; the function fills the template, fits the
requested models, and writes a `.qmd` (and optionally a rendered
`.html`) to disk.

The split is:

- **Inputs** — data, columns, priors, null range. Always required.
- **Narrative slots** — Background, Hypothesis, Inclusion/Exclusion,
  etc. Optional; sensible defaults.
- **Section toggles** — turn whole analyses (RoBMA, meta-regression,
  INSPECT-SR, …) on or off.

Once rendered, the output is a normal `.qmd` you can edit by hand.

------------------------------------------------------------------------

## Quick start

``` r
library(bayesma)

bayesma_report(
  data         = eeg_data,
  studyvar     = "Author",
  yearvar      = "Year",
  estimand     = "OR",
  event_ctrl   = "c_event",
  event_int    = "i_event",
  n_ctrl       = "c_n",
  n_int        = "i_n",
  outcome_name = "Postoperative Delirium",
  null_range   = c(0.9, 1.1),
  title        = "EEG-guided anaesthesia and POD",
  author       = "Dr Benjamin Moran",
  output_file  = "eeg_pod_report.qmd",
  render       = TRUE
)
```

The call returns the path to the generated `.qmd`. Alongside it, a
companion `*_bundle.rds` is written containing `data`, `rob_data`, and
`priors` — the template reads from this at render time.

With `render = TRUE`, the `.qmd` is rendered to HTML next to itself.

------------------------------------------------------------------------

## Required arguments

| Argument | Purpose |
|----|----|
| `data` | A data frame with one row per study (or per arm). |
| `studyvar` | Column with study identifiers. |
| `yearvar` | Column with publication year. |
| `estimand` | Target effect: `"OR"`, `"RR"`, `"HR"`, `"IRR"`, `"MD"`, `"SMD"` (relative / mean-difference) or `"RD"` / `"ARR"`, `"ATE"`, `"ATT"`, `"CATE"` (marginal — computed via [`bayesma_marginal()`](https://blmoran.github.io/bayesma/reference/bayesma_marginal.md)). |
| `outcome_name` | Human-readable outcome label (used in titles and prose). |
| `null_range` | Numeric length-2 vector. ROPE / null range on the natural scale (e.g. `c(0.9, 1.1)` for an OR). |
| `output_file` | Path for the generated `.qmd`. |

**Estimand-specific columns** — pass whichever set matches the data
shape implied by `estimand`. Names match
[`bayesma()`](https://blmoran.github.io/bayesma/reference/bayesma.md)
itself:

``` r
# Binary / count estimands (OR, RR, HR, IRR, RD, ARR, ATE, ATT, CATE)
event_ctrl = "c_event", event_int = "i_event",
n_ctrl     = "c_n",     n_int     = "i_n"

# Continuous estimands (MD, SMD)
mean_ctrl = "c_mean", mean_int = "i_mean",
sd_ctrl   = "c_sd",   sd_int   = "i_sd",
n_ctrl    = "c_n",    n_int    = "i_n"
```

------------------------------------------------------------------------

## Narrative slots

These are character strings (or character vectors) that fill prose
blocks. All are optional; defaults are generic.

| Slot | Where it appears |
|----|----|
| `background` | Introduction → Background. One or more paragraphs. |
| `hypothesis` | Methods → Hypothesis and Aims. |
| `inclusion_criteria` | Methods → Trial Design → Inclusion. Bulleted list (character vector). |
| `exclusion_criteria` | Methods → Trial Design → Exclusion. |
| `quality_tool` | Methods → Quality Assessment (default: `"Cochrane RoB 2"`). |
| `prospero_id` | Methods → Trial Design (registration paragraph). |
| `software_note` | Methods → Software (default cites `bayesma`, `cmdstanr`, R version). |

Example:

``` r
bayesma_report(
  ...,
  background = c(
    "Postoperative delirium (POD) is a common complication ...",
    "Processed EEG monitors have been proposed to ..."
  ),
  hypothesis = "pEEG-guided anaesthesia reduces the incidence of POD compared to usual care.",
  inclusion_criteria = c(
    "Randomised controlled trials (parallel-group or crossover)",
    "Adults undergoing general anaesthesia for any surgery",
    "EEG-guided anaesthesia versus standard care"
  ),
  exclusion_criteria = c(
    "Non-RCT or no full text",
    "Sedation other than general anaesthesia",
    "Non-EEG depth-of-anaesthesia monitors"
  ),
  prospero_id = "CRD420251081478"
)
```

------------------------------------------------------------------------

## Section toggles

Each toggle defaults to a sensible value (`TRUE` for primary analyses,
`FALSE` for the heavier optional ones). Setting a toggle to `FALSE`
removes both the Methods description and the Results section for that
analysis — there are no orphan headings.

| Toggle | Default | Adds |
|----|----|----|
| `run_inspect_sr` | `FALSE` | INSPECT-SR trustworthiness assessment |
| `run_rob` | `FALSE` | RoB 2 traffic-light + summary |
| `run_egger` | `TRUE` | Bayesian Egger regression + funnel plot |
| `run_robma` | `TRUE` | Single-component bias-adjusted models + RoBMA |
| `run_metareg` | `FALSE` | Meta-regression on `moderators` |
| `run_subgroup` | `FALSE` | Subgroup analyses on `subgroup_var` |
| `run_sensitivity_priors` | `TRUE` | Prior sensitivity (vague / weakly-regularising / informative) |

Toggles that take companion arguments:

``` r
bayesma_report(
  ...,
  run_metareg  = TRUE,  moderators   = c("MeanAge", "PercentFemale"),
  run_subgroup = TRUE,  subgroup_var = c("SurgeryType", "Region"),
  run_rob      = TRUE,  rob_data     = my_rob_df
)
```

`subgroup_var` accepts a character vector — one heading is generated per
variable, with sub-headings per level.

------------------------------------------------------------------------

## Models fitted automatically

Regardless of toggles, the primary analysis section fits this candidate
suite:

- One-stage common-effect (`one_stage_ce`)
- Two-stage common-effect (`two_stage_ce`)
- One-stage random-effects (`one_stage_re`)
- Two-stage random-effects (`two_stage_re`)
- Two-stage random-effects with HKSJ (`two_stage_re_hksj`)
- Two-stage random-effects with mixture RE distribution
  (`two_stage_re_mixture`)

The **primary model is selected via LOSO-CV**
([`compare_models()`](https://blmoran.github.io/bayesma/reference/compare_models.md)),
with diagnostics shown alongside. The selected model drives the forest
plot, posterior, and any downstream sensitivity / RoBMA work.

To override automatic selection:

``` r
bayesma_report(..., primary_model = "two_stage_re_hksj")
```

When `run_robma = TRUE`, an additional set of single-component
bias-adjusted models is fitted alongside RoBMA itself:

- `bias_corrected` — parametric bias-correction (Verde 2021)
- `bc_bnp` — Bias-Corrected Bayesian Non-Parametric (Verde 2025)
- `selection_copas` — Robust Bayesian Copas selection model (Bai et
  al. 2020)
- `selection_weight` — Vevea-Hedges weight-function selection model
- `pet_peese` — small-study adjustment via funnel-plot regression
- `mixture_model` — Bayesian meta-analytic mixture (Maier 2024)

------------------------------------------------------------------------

## Priors

By default the report fits with `bayesma`’s weakly-regularising priors
and reports a prior sensitivity panel covering `vague`, `weak_reg`, and
`informative` specifications.

To pass custom priors:

``` r
my_priors <- list(
  weak_reg = list(
    mu_prior  = normal(0, 1),
    tau_prior = half_student_t(3, 0, 0.5)
  ),
  informative = list(
    mu_prior  = normal(-0.3, 0.3),
    tau_prior = half_normal(0, 0.3)
  )
)

bayesma_report(..., priors = my_priors)
```

Priors are stashed in the companion `*_bundle.rds` and read by the
template — no need to inline them into YAML.

------------------------------------------------------------------------

## Marginal estimands (RD / ATE / ATT / CATE)

When `estimand` is one of `"RD"`/`"ARR"`, `"ATE"`, `"ATT"`, or `"CATE"`,
each fitted model also carries a `$marginal` element produced by
[`bayesma_marginal()`](https://blmoran.github.io/bayesma/reference/bayesma_marginal.md):

- **RD / ARR** and **ATE** for one-stage binomial fits are computed from
  posterior draws of the per-arm linear predictor.
- **RD / ARR** and **ATE** for two-stage binomial fits require a
  reference baseline. Pass `baseline_risk = 0.15` (numeric) or
  `baseline_risk = "study_mean"` (the unweighted mean of observed
  control rates).
- **ATT** is computed as the ATE weighted by treated-arm sample size;
  without IPD this is an arm-size-weighted ATE rather than a true causal
  ATT.
- **CATE** routes to a meta-regression on `cate_covariate`. Set
  `cate_covariate = "MeanAge"` (or similar). The user evaluates the
  posterior at the covariate value of interest downstream.

The relative-effect estimands (`OR`/`RR`/`HR`/`IRR`/`MD`/`SMD`) behave
exactly as before.

## Risk of bias

Pass a pre-scored RoB data frame via `rob_data` when `run_rob = TRUE`:

``` r
bayesma_report(..., run_rob = TRUE, rob_data = rob_2_scored)
```

If `rob_data` is `NULL` but `data` itself contains `rob_*` columns,
those are used instead. Otherwise the RoB section is skipped with a
message.

------------------------------------------------------------------------

## Worked example

Reproducing an EEG-guided anaesthesia report:

``` r
library(bayesma)
library(readxl)

eeg <- read_excel("EEG_POD_data.xlsx", sheet = "delirium")

bayesma_report(
  data         = eeg,
  studyvar     = "Author",
  yearvar      = "Year",
  estimand     = "OR",
  event_ctrl   = "c_event", event_int = "i_event",
  n_ctrl       = "c_n",     n_int     = "i_n",
  outcome_name = "Postoperative Delirium",
  null_range   = c(0.9, 1.1),

  title       = "EEG-guided anaesthesia for reducing POD",
  subtitle    = "A Bayesian meta-analysis",
  author      = "Dr Benjamin Moran",
  prospero_id = "CRD420251081478",

  background = c(
    "Postoperative delirium is a common complication of surgery in older adults ...",
    "Processed EEG monitors have been proposed to guide anaesthetic depth ..."
  ),
  hypothesis = "pEEG-guided anaesthesia reduces the incidence of POD compared to usual care.",
  inclusion_criteria = c(
    "Randomised controlled trials",
    "Adults undergoing general anaesthesia",
    "pEEG-guided versus standard care"
  ),
  exclusion_criteria = c(
    "Non-RCT or no full text",
    "Sedation other than general anaesthesia",
    "Non-EEG depth monitors"
  ),

  run_inspect_sr         = FALSE,
  run_rob                = TRUE,
  run_egger              = TRUE,
  run_robma              = TRUE,
  run_subgroup           = TRUE,
  subgroup_var           = c("SurgeryType", "Region"),
  run_sensitivity_priors = TRUE,

  output_file = "eeg_pod_report.qmd",
  render      = TRUE
)
```

------------------------------------------------------------------------

## Editing the output

The generated `.qmd` is a regular Quarto document. Common post-render
edits:

- Tighten the prose slots (Background, Limitations, Conclusion) — these
  are the parts a templated report cannot get exactly right.
- Reorder Results sections, drop a model that didn’t converge, add
  domain-specific commentary.
- Swap the `format:` block for `pdf` or `docx` when submitting.

Re-rendering after edits is just `quarto render eeg_pod_report.qmd` —
[`bayesma_report()`](https://blmoran.github.io/bayesma/reference/bayesma_report.md)
does not need to be called again. (The companion `*_bundle.rds` must
remain alongside the `.qmd`.)

------------------------------------------------------------------------

## What `bayesma_report()` does *not* do

- It does not write your Background / Discussion / Limitations for you.
  The slots accept your prose verbatim; the function does not generate
  clinical narrative.
- It does not assess study quality. Pass a pre-scored RoB data frame via
  `rob_data =` if `run_rob = TRUE`.
- It does not search literature, screen studies, or extract data — it
  starts from a clean tibble.

------------------------------------------------------------------------

## See also

- [`workflow-guide`](https://blmoran.github.io/bayesma/articles/workflow-guide.qmd)
  — full walk-through of the underlying functions.
- [`robma`](https://blmoran.github.io/bayesma/articles/robma.qmd) —
  RoBMA-specific details.
- [`primed`](https://blmoran.github.io/bayesma/articles/primed.qmd) —
  preliminary data exploration.
