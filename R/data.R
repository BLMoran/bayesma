#' Example INSPECT-SR Dataset: EEG-Guided Anaesthesia and Delirium
#'
#' A simulated eight-trial systematic review used to demonstrate the full
#' INSPECT-SR workflow. The eight trials are hand-crafted to span a range of
#' trustworthiness profiles and to exercise every automated check.
#'
#' @format A tibble with one row per study. Every INSPECT-SR item has its own
#'   column so that the expected layout is immediately obvious. Variable-length
#'   inputs to the automated Domain 4 checks (Table 1 baselines and reported
#'   test statistics) are stored in list-columns.
#' \describe{
#'   \item{`study`}{Character. Study identifier (must be unique).}
#'   \item{`n_randomised_int`, `n_randomised_ctrl`, `n_randomised_total`,
#'         `n_analysed_int`, `n_analysed_ctrl`,
#'         `n_lost_int`, `n_lost_ctrl`}{Integer participant-flow columns
#'     feeding check 4.6.}
#'   \item{`d1_1`..`d1_3`, `d2_1`..`d2_5`, `d3_1`..`d3_2`,
#'         `d4_1`..`d4_11`}{Character. Manual judgement for each INSPECT-SR
#'     item. Allowed values: `"No concerns"`, `"Some concerns"`,
#'     `"Serious concerns"`, `NA`. The four automated items (`d4_3`, `d4_6`,
#'     `d4_8`, `d4_9`) are left as `NA` and filled in by [inspect_sr()].}
#'   \item{`d1_1_comment`..`d4_11_comment`}{Character. Optional free-text
#'     comment attached to each judgement.}
#'   \item{`baseline`}{List-column. Each element is a data frame of Table 1
#'     variables with columns `variable`, `mean_int`, `sd_int`, `mean_ctrl`,
#'     `sd_ctrl`, `p_value`, `integer_scale`. Feeds GRIM (4.8) and Carlisle
#'     (4.3).}
#'   \item{`statistics`}{List-column. Each element is a data frame of reported
#'     tests with columns `test_type`, `statistic`, `df`, `df2`, `reported_p`,
#'     `context`. Feeds the p-value check (4.9).}
#'   \item{`outcome_estimate`, `outcome_ci_lower`, `outcome_ci_upper`,
#'         `outcome_se`, `outcome_log_scale`}{Primary outcome (odds ratio for
#'     delirium incidence) — values are on the log scale.}
#' }
#'
#' @details
#' The eight trials, their expected behaviour under the automated checks, and
#' the INSPECT-SR item they illustrate:
#'
#' \tabular{lll}{
#'   **Study**         \tab **Profile**                               \tab **Expected flag** \cr
#'   Callahan (1994)   \tab Clean small trial                         \tab None \cr
#'   Peek (2005)       \tab Clean large multi-centre                  \tab None \cr
#'   Clayton (2010)    \tab Clean                                     \tab None \cr
#'   Singer (2003)     \tab Clean (large trial)                       \tab None \cr
#'   Underwood (2013)  \tab Some concerns (late registration)         \tab Manual D2 \cr
#'   Rosa (2008)       \tab Too-perfect baseline balance              \tab Carlisle (4.3) \cr
#'   Doe (1995)        \tab Retracted; GRIM failures on integer scale \tab GRIM (4.8), manual D1/D2/D3 \cr
#'   Kint (1995)       \tab CONSORT arithmetic off                    \tab N consistency (4.6) \cr
#' }
#'
#' All values are fabricated for didactic purposes. The layout matches the
#' schema expected by [inspect_sr()], so the whole review can be assessed
#' in one call.
#'
#' @source Simulated.
#'
#' @examples
#' data(inspect_sr_example)
#'
#' # Frequentist assessment
#' res <- inspect_sr(inspect_sr_example, verbose = FALSE)
#' res
#'
#' # Per-check table
#' inspect_sr_table(res)
#'
#' \dontrun{
#' # Bayesian assessment
#' bres <- inspect_sr(inspect_sr_example, bayes = TRUE, verbose = FALSE)
#' bres
#' }
"inspect_sr_example"

#' Example Binary Outcome Dataset
#'
#' A simulated 12-study meta-analysis of a binary outcome (postoperative nausea
#' and vomiting) used to demonstrate bayesma workflows for binomial likelihoods.
#'
#' @format A tibble with 12 rows and 18 columns:
#' \describe{
#'   \item{`Author`}{Character. Study author/identifier.}
#'   \item{`Year`}{Numeric. Publication year.}
#'   \item{`Subgroup`}{Character. Surgical specialty subgroup.}
#'   \item{`Control`}{Character. Control intervention label.}
#'   \item{`Intervention`}{Character. Active intervention label.}
#'   \item{`N_Total`}{Numeric. Total sample size.}
#'   \item{`N_Control`}{Numeric. Control arm sample size.}
#'   \item{`N_Intervention`}{Numeric. Intervention arm sample size.}
#'   \item{`Event_Control`}{Numeric. Events in control arm.}
#'   \item{`Outcome_Control_No`}{Numeric. Non-events in control arm.}
#'   \item{`Event_Intervention`}{Numeric. Events in intervention arm.}
#'   \item{`Outcome_Intervention_No`}{Numeric. Non-events in intervention arm.}
#'   \item{`D1`, `D2`, `D3`, `D4`, `D5`}{Character. Risk-of-bias domain
#'     judgements: `"Low"`, `"Some concerns"`, or `"High"`.}
#'   \item{`Overall`}{Character. Overall risk-of-bias judgement.}
#' }
#'
#' @source Simulated.
"binary_outcome"

#' Example Continuous Outcome Dataset
#'
#' A simulated 12-study meta-analysis of a continuous outcome used to
#' demonstrate bayesma workflows for Gaussian likelihoods.
#'
#' @format A tibble with 12 rows and 15 columns:
#' \describe{
#'   \item{`Author`}{Character. Study author/identifier.}
#'   \item{`Year`}{Numeric. Publication year.}
#'   \item{`Subgroup`}{Character. Subgroup label.}
#'   \item{`Mean_Control`}{Numeric. Control arm mean.}
#'   \item{`SD_Control`}{Numeric. Control arm standard deviation.}
#'   \item{`N_Control`}{Numeric. Control arm sample size.}
#'   \item{`Mean_Intervention`}{Numeric. Intervention arm mean.}
#'   \item{`SD_Intervention`}{Numeric. Intervention arm standard deviation.}
#'   \item{`N_Intervention`}{Numeric. Intervention arm sample size.}
#'   \item{`D1`, `D2`, `D3`, `D4`, `D5`}{Character. Risk-of-bias domain
#'     judgements: `"Low"`, `"Some concerns"`, or `"High"`.}
#'   \item{`Overall`}{Character. Overall risk-of-bias judgement.}
#' }
#'
#' @source Simulated.
"cont_outcome"
