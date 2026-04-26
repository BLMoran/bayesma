# Models overview

## Bayesian meta-analysis models in **bayesma**

The **bayesma** package implements a range of Bayesian meta-analysis
models for synthesizing evidence across studies.

These models differ in how they handle:

- between-study heterogeneity
- study-level bias
- publication bias
- small-study effects
- model uncertainty

### Available model families

| Model                | Purpose                                      |
|----------------------|----------------------------------------------|
| Common-effect model  | Assumes a single true effect across studies  |
| Random-effects model | Allows study effects to vary                 |
| One-stage models     | Uses full likelihood of the original data    |
| Two-stage models     | Uses study-level effect estimates            |
| Bias-corrected model | Adjusts for systematic study bias            |
| Selection models     | Adjusts for publication bias                 |
| PET-PEESE            | Regression-based small-study bias correction |
| Robust mixtures      | Down-weights outliers                        |
| RoBMA                | Model-averaging across bias models           |

Each article below explains the statistical model, its assumptions, and
how it works.
