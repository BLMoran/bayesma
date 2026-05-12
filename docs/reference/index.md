# Package index

## Core meta-analysis

Functions for fitting common-effect and random-effects meta-analysis
models.

- [`bayesma()`](https://blmoran.github.io/bayesma/reference/bayesma.md)
  : Run a Bayesian Meta-Analysis in Stan
- [`bayesma_marginal()`](https://blmoran.github.io/bayesma/reference/bayesma_marginal.md)
  : Compute a marginal estimand from a bayesma fit
- [`bayesma_report()`](https://blmoran.github.io/bayesma/reference/bayesma_report.md)
  : Generate a Bayesian meta-analysis report

## Multivariate meta-analysis

- [`bayesma_mv()`](https://blmoran.github.io/bayesma/reference/bayesma_mv.md)
  : Run a Multivariate Bayesian Meta-Analysis in Stan

- [`bayesma_mv_spec()`](https://blmoran.github.io/bayesma/reference/bayesma_mv_spec.md)
  : Build a bivariate meta-analysis specification object

- [`bayesma_mv_stan_code()`](https://blmoran.github.io/bayesma/reference/bayesma_mv_stan_code.md)
  : Generate Stan code for a bivariate meta-analysis specification

- [`bayesma_mv_stan_data()`](https://blmoran.github.io/bayesma/reference/bayesma_mv_stan_data.md)
  : Build the Stan data list for a bivariate meta-analysis specification

- [`bayesma_mv_fit()`](https://blmoran.github.io/bayesma/reference/bayesma_mv_fit.md)
  : Compile and sample a bivariate meta-analysis model

- [`bayesma_mv_extract()`](https://blmoran.github.io/bayesma/reference/bayesma_mv_extract.md)
  : Extract tidy effect components from a bivariate meta-analysis fit

- [`bayesma_mv_output()`](https://blmoran.github.io/bayesma/reference/bayesma_mv_output.md)
  :

  Assemble a `bayesma_mv` object from pipeline outputs

## Egger’s test

- [`egger()`](https://blmoran.github.io/bayesma/reference/egger.md) :
  Egger's Regression Test for Small-Study Effects (Bayesian)
- [`egger_plot()`](https://blmoran.github.io/bayesma/reference/egger_plot.md)
  : Plot method for bayesma_egger

## Meta-regression

- [`meta_reg()`](https://blmoran.github.io/bayesma/reference/meta_reg.md)
  : Bayesian Meta-Regression
- [`coefficient_evidence()`](https://blmoran.github.io/bayesma/reference/coefficient_evidence.md)
  : Bayesian Evidence for Meta-Regression Coefficients
- [`metareg_mod_plot()`](https://blmoran.github.io/bayesma/reference/metareg_mod_plot.md)
  : Plot Method for bayesma_coef_evidence Objects
- [`bubble_plot()`](https://blmoran.github.io/bayesma/reference/bubble_plot.md)
  : Bubble Plot for Meta-Regression
- [`multi_bubble_plots()`](https://blmoran.github.io/bayesma/reference/multi_bubble_plots.md)
  : Multi-panel Bubble Plots

## Model comparison

- [`compare_models()`](https://blmoran.github.io/bayesma/reference/compare_models.md)
  : Compare Multiple Bayesian Meta-Analysis Models
- [`compare_table()`](https://blmoran.github.io/bayesma/reference/compare_table.md)
  : Create a Table for Model Comparison Results
- [`compare_plot()`](https://blmoran.github.io/bayesma/reference/compare_plot.md)
  : Plot Model Comparison Results

## RoBMA

- [`robma()`](https://blmoran.github.io/bayesma/reference/robma.md) :
  Robust Bayesian Model Averaging for Meta-Analysis
- [`robma_table()`](https://blmoran.github.io/bayesma/reference/robma_table.md)
  : Create a gt Table for RoBMA Results
- [`robma_sensitivity()`](https://blmoran.github.io/bayesma/reference/robma_sensitivity.md)
  : Run RoBMA Models Across Multiple Prior Specifications
- [`robma_default_priors()`](https://blmoran.github.io/bayesma/reference/robma_default_priors.md)
  : Default RoBMA prior set
- [`attach_robma_sensitivity()`](https://blmoran.github.io/bayesma/reference/attach_robma_sensitivity.md)
  : Attach RoBMA Sensitivity Fits to a bayesma Object
- [`has_robma_sensitivity()`](https://blmoran.github.io/bayesma/reference/has_robma_sensitivity.md)
  : Check if RoBMA Sensitivity Fits are Available

## Prior constructors

Functions for specifying prior distributions.

- [`normal()`](https://blmoran.github.io/bayesma/reference/priors.md)
  [`half_normal()`](https://blmoran.github.io/bayesma/reference/priors.md)
  [`cauchy()`](https://blmoran.github.io/bayesma/reference/priors.md)
  [`half_cauchy()`](https://blmoran.github.io/bayesma/reference/priors.md)
  [`half_student_t()`](https://blmoran.github.io/bayesma/reference/priors.md)
  [`exponential()`](https://blmoran.github.io/bayesma/reference/priors.md)
  [`uniform()`](https://blmoran.github.io/bayesma/reference/priors.md)
  [`dirichlet()`](https://blmoran.github.io/bayesma/reference/priors.md)
  [`beta()`](https://blmoran.github.io/bayesma/reference/priors.md)
  [`scaled_inv_chi_sq()`](https://blmoran.github.io/bayesma/reference/priors.md)
  [`lkj()`](https://blmoran.github.io/bayesma/reference/priors.md) :
  Prior distribution constructors
- [`prior_bias()`](https://blmoran.github.io/bayesma/reference/prior_bias.md)
  [`prior_weight_function()`](https://blmoran.github.io/bayesma/reference/prior_bias.md)
  [`prior_pet()`](https://blmoran.github.io/bayesma/reference/prior_bias.md)
  [`prior_peese()`](https://blmoran.github.io/bayesma/reference/prior_bias.md)
  [`prior_copas()`](https://blmoran.github.io/bayesma/reference/prior_bias.md)
  [`prior_jung()`](https://blmoran.github.io/bayesma/reference/prior_bias.md)
  [`prior_no_bias()`](https://blmoran.github.io/bayesma/reference/prior_bias.md)
  : RoBMA bias-prior constructors

## INSPECT-SR

Tools for assessing the trustworthiness of included trials.

- [`inspect_sr()`](https://blmoran.github.io/bayesma/reference/inspect_sr.md)
  : Run the INSPECT-SR Trustworthiness Assessment
- [`inspect_sr_table()`](https://blmoran.github.io/bayesma/reference/inspect_sr_table.md)
  : Per-Check INSPECT-SR Results Table
- [`inspect_plot()`](https://blmoran.github.io/bayesma/reference/inspect_plot.md)
  : Create INSPECT-SR Trustworthiness Plot
- [`inspect_summary_plot()`](https://blmoran.github.io/bayesma/reference/inspect_summary_plot.md)
  : Create INSPECT-SR Summary Bar Plot
- [`filter_trustworthy()`](https://blmoran.github.io/bayesma/reference/filter_trustworthy.md)
  : Filter Studies by INSPECT-SR Trustworthiness
- [`bayes_carlisle_test()`](https://blmoran.github.io/bayesma/reference/bayes_carlisle_test.md)
  : Bayesian Carlisle's Test for Baseline Balance
- [`bayes_grim_test()`](https://blmoran.github.io/bayesma/reference/bayes_grim_test.md)
  : Bayesian GRIM Test
- [`bayes_verify_pvalue()`](https://blmoran.github.io/bayesma/reference/bayes_verify_pvalue.md)
  : Bayesian P-Value Verification
- [`carlisle_test()`](https://blmoran.github.io/bayesma/reference/carlisle_test.md)
  : Carlisle's Test for Baseline Balance
- [`grim_test()`](https://blmoran.github.io/bayesma/reference/grim_test.md)
  : GRIM Test (Granularity-Related Inconsistency of Means)
- [`verify_pvalue()`](https://blmoran.github.io/bayesma/reference/verify_pvalue.md)
  : Verify a Reported P-Value
- [`check_n_consistency()`](https://blmoran.github.io/bayesma/reference/check_n_consistency.md)
  : Check Participant Number Consistency
- [`check_statistics_consistency()`](https://blmoran.github.io/bayesma/reference/check_statistics_consistency.md)
  : Check Internal Consistency of Summary Statistics
- [`inspect_sr_example`](https://blmoran.github.io/bayesma/reference/inspect_sr_example.md)
  : Example INSPECT-SR Dataset: EEG-Guided Anaesthesia and Delirium

## PRIMED

Preliminary investigation of meta-analytic data structure.

- [`primed()`](https://blmoran.github.io/bayesma/reference/primed.md) :
  PRIMED: Preliminary Investigation of Meta-analytic Databases

## Sensitivity analysis

- [`sensitivity_plot()`](https://blmoran.github.io/bayesma/reference/sensitivity_plot.md)
  : Generate Sensitivity Analysis Plot for Bayesian Meta-Analysis
- [`render_sensitivity_patchwork()`](https://blmoran.github.io/bayesma/reference/render_sensitivity_patchwork.md)
  : Render the stored components into a patchwork object
- [`sens_add_probs()`](https://blmoran.github.io/bayesma/reference/sens_add.md)
  [`sens_add_null()`](https://blmoran.github.io/bayesma/reference/sens_add.md)
  [`sens_add_titles()`](https://blmoran.github.io/bayesma/reference/sens_add.md)
  [`sens_add_x_lim()`](https://blmoran.github.io/bayesma/reference/sens_add.md)
  [`sens_add_plot_width()`](https://blmoran.github.io/bayesma/reference/sens_add.md)
  : Post-render modifications to a sensitivity plot

## Diagnostics

- [`diagnostics()`](https://blmoran.github.io/bayesma/reference/diagnostics.md)
  : Single-Page Model Diagnostics for bayesma Objects
- [`pp_check()`](https://blmoran.github.io/bayesma/reference/pp_check.md)
  : Posterior and Prior Predictive Checks for bayesma Objects

## Visualisations

- [`forest()`](https://blmoran.github.io/bayesma/reference/forest.md) :
  Create Bayesian Forest Plot for Meta-Analysis
- [`overall_plot()`](https://blmoran.github.io/bayesma/reference/overall_plot.md)
  : Create posterior plots for Bayesian meta-analysis
- [`funnel_plot()`](https://blmoran.github.io/bayesma/reference/funnel_plot.md)
  : Create a Funnel Plot for Bayesian Meta-Analysis
- [`rob_plot()`](https://blmoran.github.io/bayesma/reference/rob_plot.md)
  : Create Risk of Bias Plot
- [`ecdf_model_plot()`](https://blmoran.github.io/bayesma/reference/ecdf_model_plot.md)
  : ECDF Plot Comparing Model Strategies
- [`ecdf_prior_plot()`](https://blmoran.github.io/bayesma/reference/ecdf_prior_plot.md)
  : ECDF Plot Comparing Prior Sensitivity

## Datasets

- [`binary_outcome`](https://blmoran.github.io/bayesma/reference/binary_outcome.md)
  : Example Binary Outcome Dataset
- [`cont_outcome`](https://blmoran.github.io/bayesma/reference/cont_outcome.md)
  : Example Continuous Outcome Dataset

## Data entry

Tools for entering, formatting, and downloading study data.

- [`study_table()`](https://blmoran.github.io/bayesma/reference/study_table.md)
  : Create a Table of Included Studies
- [`download_template()`](https://blmoran.github.io/bayesma/reference/download_template.md)
  : Download a Data Entry Template
- [`launch_data_entry()`](https://blmoran.github.io/bayesma/reference/launch_data_entry.md)
  : Launch the bayesma Data Entry App

## Utilities

- [`stan_code()`](https://blmoran.github.io/bayesma/reference/stan_code.md)
  : Print Stan code from a fitted model
- [`interpret()`](https://blmoran.github.io/bayesma/reference/interpret.md)
  : Interpret a Bayesian meta-analysis workflow
