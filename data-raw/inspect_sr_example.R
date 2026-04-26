# ============================================================================
# Build `inspect_sr_example` — a single tibble with one row per study.
#
# Every INSPECT-SR item (1.1–4.11) has its own column so users can see
# exactly what needs to be filled in.  Variable-length data that feeds the
# automated Domain 4 checks (baseline characteristics, test statistics,
# outcomes) is stored in list-columns.
#
# Re-run with:
#   source("data-raw/inspect_sr_example.R")
# ============================================================================

library(usethis)
library(tibble)

# --- Helper: build baseline data frames -----------------------------------

bl_callahan <- data.frame(
  variable      = c("Age", "BMI", "ASA_score"),
  mean_int      = c(71.2, 24.8, 2.5),
  sd_int        = c(8.3, 4.1, 0.5),
  mean_ctrl     = c(69.8, 25.1, 2.5),
  sd_ctrl       = c(9.1, 3.9, 0.5),
  p_value       = c(0.45, 0.72, 0.34),
  integer_scale = c(FALSE, FALSE, TRUE),
  stringsAsFactors = FALSE
)

bl_peek <- data.frame(
  variable      = c("Age", "BMI", "Duration_surgery_min", "Male_pct"),
  mean_int      = c(68.4, 26.7, 185, 62),
  sd_int        = c(10.2, 5.1, 55, NA),
  mean_ctrl     = c(67.9, 27.1, 190, 59),
  sd_ctrl       = c(9.8, 4.8, 60, NA),
  p_value       = c(0.58, 0.39, 0.37, 0.51),
  integer_scale = c(FALSE, FALSE, TRUE, FALSE),
  stringsAsFactors = FALSE
)

bl_clayton <- data.frame(
  variable      = c("Age", "MMSE_score", "Charlson_index"),
  mean_int      = c(74.1, 27.8, 3.2),
  sd_int        = c(7.2, 2.1, 1.8),
  mean_ctrl     = c(73.5, 27.5, 3.4),
  sd_ctrl       = c(7.8, 2.3, 1.9),
  p_value       = c(0.57, 0.33, 0.43),
  integer_scale = c(FALSE, FALSE, FALSE),
  stringsAsFactors = FALSE
)

bl_singer <- data.frame(
  variable      = c("Age", "BMI", "ASA_3_4_pct"),
  mean_int      = c(69.2, 28.4, 78),
  sd_int        = c(6.1, 5.7, NA),
  mean_ctrl     = c(69.5, 28.1, 76),
  sd_ctrl       = c(6.3, 5.5, NA),
  p_value       = c(0.38, 0.32, 0.42),
  integer_scale = c(FALSE, FALSE, FALSE),
  stringsAsFactors = FALSE
)

bl_underwood <- data.frame(
  variable      = c("Age", "BMI", "Duration_anaesthesia_min"),
  mean_int      = c(66.5, 23.1, 210),
  sd_int        = c(9.2, 3.4, 65),
  mean_ctrl     = c(67.1, 23.5, 198),
  sd_ctrl       = c(8.8, 3.6, 58),
  p_value       = c(0.68, 0.44, 0.21),
  integer_scale = c(FALSE, FALSE, TRUE),
  stringsAsFactors = FALSE
)

bl_rosa <- data.frame(
  variable      = c("Age", "BMI", "SBP", "HR", "Albumin", "Haemoglobin"),
  mean_int      = c(68.2, 24.1, 135.6, 74.2, 38.1, 128.3),
  sd_int        = c(9.1, 3.2, 15.8, 10.3, 4.2, 15.1),
  mean_ctrl     = c(68.1, 24.0, 135.5, 74.1, 38.0, 128.2),
  sd_ctrl       = c(9.0, 3.1, 15.7, 10.2, 4.1, 15.0),
  p_value       = c(0.93, 0.81, 0.95, 0.94, 0.85, 0.95),
  integer_scale = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
  stringsAsFactors = FALSE
)

bl_doe <- data.frame(
  variable      = c("Age", "ASA_score", "Pain_VAS", "Barthel_index"),
  mean_int      = c(72.3, 2.45, 4.75, 85.7),
  sd_int        = c(8.1, 0.5, 1.8, 12.3),
  mean_ctrl     = c(71.8, 2.55, 4.65, 86.2),
  sd_ctrl       = c(7.9, 0.5, 1.9, 11.8),
  p_value       = c(0.67, 0.58, 0.78, 0.78),
  integer_scale = c(FALSE, TRUE, TRUE, FALSE),
  stringsAsFactors = FALSE
)

bl_kint <- data.frame(
  variable      = c("Age", "BMI", "Duration_surgery_min"),
  mean_int      = c(70.5, 25.3, 175),
  sd_int        = c(8.6, 4.0, 50),
  mean_ctrl     = c(69.8, 25.7, 180),
  sd_ctrl       = c(8.2, 4.2, 55),
  p_value       = c(0.49, 0.42, 0.40),
  integer_scale = c(FALSE, FALSE, TRUE),
  stringsAsFactors = FALSE
)

# --- Helper: build statistics data frames ---------------------------------

st_callahan <- data.frame(
  test_type = "chi_sq", statistic = 1.23, df = 1, df2 = NA_real_,
  reported_p = 0.27, context = "Delirium incidence",
  stringsAsFactors = FALSE
)

st_peek <- data.frame(
  test_type  = c("chi_sq", "t"),
  statistic  = c(4.12, -0.85),
  df         = c(1, 437),
  df2        = c(NA_real_, NA_real_),
  reported_p = c(0.04, 0.40),
  context    = c("Delirium incidence", "Duration of delirium"),
  stringsAsFactors = FALSE
)

st_clayton <- data.frame(
  test_type = "chi_sq", statistic = 2.78, df = 1, df2 = NA_real_,
  reported_p = 0.10, context = "Delirium incidence",
  stringsAsFactors = FALSE
)

st_singer <- data.frame(
  test_type = "chi_sq", statistic = 0.18, df = 1, df2 = NA_real_,
  reported_p = 0.67, context = "Delirium incidence",
  stringsAsFactors = FALSE
)

st_underwood <- data.frame(
  test_type = "chi_sq", statistic = 5.02, df = 1, df2 = NA_real_,
  reported_p = 0.02, context = "Delirium incidence",
  stringsAsFactors = FALSE
)

st_rosa <- data.frame(
  test_type = "chi_sq", statistic = 7.84, df = 1, df2 = NA_real_,
  reported_p = 0.005, context = "Delirium incidence",
  stringsAsFactors = FALSE
)

st_doe <- data.frame(
  test_type  = c("chi_sq", "t"),
  statistic  = c(12.5, -3.21),
  df         = c(1, 173),
  df2        = c(NA_real_, NA_real_),
  reported_p = c(0.001, 0.002),
  context    = c("Delirium incidence", "Delirium duration"),
  stringsAsFactors = FALSE
)

st_kint <- data.frame(
  test_type = "chi_sq", statistic = 3.12, df = 1, df2 = NA_real_,
  reported_p = 0.08, context = "Delirium incidence",
  stringsAsFactors = FALSE
)

# ============================================================================
# Build the single tibble
# ============================================================================
#
# Column naming convention:
#   d{domain}_{item}  e.g. d1_1, d2_3, d4_8
#
# Judgement values:
#   "No concerns"      — no issue identified
#   "Some concerns"    — potential issue, warrants caution
#   "Serious concerns" — clear problem identified
#   NA                 — not assessed / not applicable
#
# INSPECT-SR signalling questions (Wilkinson et al. 2025):
#
# D1  Post-publication notices
#   1.1  Retraction, expression of concern, or correction issued?
#   1.2  Concerns raised on PubPeer, social media, or post-pub commentary?
#   1.3  Concerns raised about other trials by the same research group?
#
# D2  Conduct, governance, and transparency
#   2.1  Was the trial prospectively registered?
#   2.2  Was ethics approval reported?
#   2.3  Evidence of regulatory oversight (for drug/device trials)?
#   2.4  Is a protocol or statistical analysis plan available?
#   2.5  Has the corresponding author been responsive to queries?
#
# D3  Text and publication details
#   3.1  Concerns about text, tables, or figures (plagiarism, image reuse)?
#   3.2  Concerns about the publication history or provenance?
#
# D4  Results in the study
#   4.1  Are baseline data distributions plausible?
#   4.2  Are baseline data consistent across publications of the same trial?
#   4.3  Are any baseline data implausible (Carlisle's test)?        [AUTOMATED]
#   4.4  Are effect sizes plausible for the context?
#   4.5  Are the results "too good to be true"?
#   4.6  Unexplained inconsistencies in participant numbers?          [AUTOMATED]
#   4.7  Unexplained inconsistencies in results across publications?
#   4.8  Are means/variances of integer data impossible (GRIM)?       [AUTOMATED]
#   4.9  Are there errors in statistical results (p-value check)?     [AUTOMATED]
#   4.10 Are individual participant data available and consistent?
#   4.11 Any other concerns about the reported results?
# ============================================================================

inspect_sr_example <- tibble::tibble(

  study = c(
    "Callahan (1994)", "Peek (2005)",      "Clayton (2010)",  "Singer (2003)",
    "Underwood (2013)", "Rosa (2008)",      "Doe (1995)",      "Kint (1995)"
  ),

  # ── Participant flow numbers ──────────────────────────────────────────────
  n_randomised_int   = c( 46L, 230L, 101L, 614L,  80L, 120L,  90L, 150L),
  n_randomised_ctrl  = c( 46L, 230L, 101L, 618L,  80L, 120L,  90L, 150L),
  n_randomised_total = c( 92L, 460L, 202L, 1232L, 160L, 240L, 180L, 300L),
  n_analysed_int     = c( 43L, 221L,  98L, 598L,  76L, 118L,  88L, 142L),
  n_analysed_ctrl    = c( 44L, 218L,  96L, 605L,  74L, 117L,  87L, 138L),
  n_lost_int         = c(  3L,   9L,   3L,  16L,   4L,   2L,   2L,   5L),
  n_lost_ctrl        = c(  2L,  12L,   5L,  13L,   6L,   3L,   3L,   8L),

  # ── Domain 1: Post-publication notices ───────────────────────────────────
  d1_1 = c(
    "No concerns",      # Callahan: no retraction
    "No concerns",      # Peek
    "No concerns",      # Clayton
    "No concerns",      # Singer
    "No concerns",      # Underwood
    "No concerns",      # Rosa
    "Serious concerns", # Doe: retracted 2024
    NA_character_       # Kint: not assessed
  ),
  d1_1_comment = c(
    NA, NA, NA, NA, NA, NA,
    "Retracted 2024 for data integrity concerns",
    NA
  ),

  d1_2 = c(
    "No concerns", "No concerns", "No concerns", "No concerns",
    "No concerns", "No concerns",
    "Serious concerns", # Doe: PubPeer threads
    NA_character_
  ),
  d1_2_comment = c(
    NA, NA, NA, NA, NA, NA,
    "Multiple PubPeer comments flagging statistical anomalies",
    NA
  ),

  d1_3 = c(
    "No concerns", "No concerns", "No concerns", "No concerns",
    "No concerns", "No concerns",
    "Some concerns", # Doe: other trials by same group questioned
    NA_character_
  ),
  d1_3_comment = c(
    NA, NA, NA, NA, NA, NA,
    "Two other trials by senior author under investigation",
    NA
  ),

  # ── Domain 2: Conduct, governance, and transparency ─────────────────────
  d2_1 = c(
    "No concerns",      # Callahan: registered
    "No concerns",      # Peek: pre-registered, protocol published
    "No concerns",      # Clayton
    "No concerns",      # Singer: ENGAGES, well-documented
    "Some concerns",    # Underwood: registered after enrolment began
    "Some concerns",    # Rosa: no registration found
    "Serious concerns", # Doe: not registered
    NA_character_       # Kint: not assessed
  ),
  d2_1_comment = c(
    NA,
    "Pre-registered, protocol published",
    NA,
    "ENGAGES trial, registered and well-documented",
    "Registered after enrolment began",
    "No registration found",
    "No registration",
    NA
  ),

  d2_2 = c(
    "No concerns", "No concerns", "No concerns", "No concerns",
    "No concerns",
    "No concerns",
    "Serious concerns", # Doe: no ethics approval documented
    NA_character_
  ),
  d2_2_comment = c(
    NA, NA, NA, NA, NA, NA,
    "No ethics approval documented",
    NA
  ),

  d2_3 = c(
    NA_character_,      # Callahan: not a drug/device trial
    NA_character_,      # Peek
    NA_character_,      # Clayton
    NA_character_,      # Singer
    NA_character_,      # Underwood
    NA_character_,      # Rosa
    NA_character_,      # Doe
    NA_character_       # Kint
  ),
  d2_3_comment = c(NA, NA, NA, NA, NA, NA, NA, NA),

  d2_4 = c(
    "No concerns",      # Callahan
    "No concerns",      # Peek: protocol published
    "No concerns",      # Clayton
    "No concerns",      # Singer: protocol available
    "No concerns",      # Underwood
    "Some concerns",    # Rosa: no protocol
    "Serious concerns", # Doe: no protocol
    NA_character_       # Kint
  ),
  d2_4_comment = c(NA, "Protocol published", NA, "Protocol available", NA,
                   "No protocol located", "No protocol or SAP available", NA),

  d2_5 = c(
    "No concerns", "No concerns", "No concerns", "No concerns",
    "No concerns",
    "Some concerns",    # Rosa: author unresponsive
    "Serious concerns", # Doe: author unresponsive
    NA_character_
  ),
  d2_5_comment = c(NA, NA, NA, NA, NA,
                   "Author did not respond to two emails",
                   "Author uncontactable; affiliation could not be verified",
                   NA),

  # ── Domain 3: Text and publication details ───────────────────────────────
  d3_1 = c(
    "No concerns", "No concerns", "No concerns", "No concerns",
    NA_character_,      # Underwood: not assessed
    NA_character_,      # Rosa: not assessed
    "Some concerns",    # Doe: duplicated text
    NA_character_       # Kint
  ),
  d3_1_comment = c(
    NA, NA, NA, NA, NA, NA,
    "Duplicated text from earlier publication by same group",
    NA
  ),

  d3_2 = c(
    "No concerns", "No concerns", "No concerns", "No concerns",
    NA_character_, NA_character_,
    "Some concerns",    # Doe: publication history unclear
    NA_character_
  ),
  d3_2_comment = c(
    NA, NA, NA, NA, NA, NA,
    "Submitted simultaneously to two journals",
    NA
  ),

  # ── Domain 4: Results (manual items) ─────────────────────────────────────
  # Items 4.3, 4.6, 4.8, 4.9 are automated — leave NA here; inspect_sr()

  # fills them in from the baseline/statistics list-columns.

  d4_1 = c(
    "No concerns", "No concerns", "No concerns", "No concerns",
    "No concerns", "No concerns", "No concerns", "No concerns"
  ),
  d4_1_comment = c(NA, NA, NA, NA, NA, NA, NA, NA),

  d4_2 = c(
    NA_character_, NA_character_, NA_character_, NA_character_,
    NA_character_, NA_character_,
    "Some concerns",    # Doe: conference abstract differs from paper
    NA_character_
  ),
  d4_2_comment = c(NA, NA, NA, NA, NA, NA,
                   "Sample sizes in conference abstract differ from published paper",
                   NA),

  # d4_3: Carlisle's test — AUTOMATED (leave NA)
  d4_3 = rep(NA_character_, 8),
  d4_3_comment = rep(NA_character_, 8),

  d4_4 = c(
    "No concerns", "No concerns", "No concerns", "No concerns",
    "No concerns", "No concerns",
    "Some concerns",    # Doe: implausibly large effect
    "No concerns"
  ),
  d4_4_comment = c(NA, NA, NA, NA, NA, NA,
                   "Effect size much larger than all other trials",
                   NA),

  d4_5 = c(
    "No concerns", "No concerns", "No concerns", "No concerns",
    "No concerns", "No concerns",
    "Some concerns",    # Doe
    "No concerns"
  ),
  d4_5_comment = c(NA, NA, NA, NA, NA, NA,
                   "Very low dropout, very large treatment effect, very narrow CIs",
                   NA),

  # d4_6: Participant number consistency — AUTOMATED (leave NA)
  d4_6 = rep(NA_character_, 8),
  d4_6_comment = rep(NA_character_, 8),

  d4_7 = rep(NA_character_, 8),
  d4_7_comment = rep(NA_character_, 8),

  # d4_8: GRIM test — AUTOMATED (leave NA)
  d4_8 = rep(NA_character_, 8),
  d4_8_comment = rep(NA_character_, 8),

  # d4_9: P-value verification — AUTOMATED (leave NA)
  d4_9 = rep(NA_character_, 8),
  d4_9_comment = rep(NA_character_, 8),

  d4_10 = rep(NA_character_, 8),
  d4_10_comment = rep(NA_character_, 8),

  d4_11 = rep(NA_character_, 8),
  d4_11_comment = rep(NA_character_, 8),

  # ── Data for automated Domain 4 checks (list-columns) ───────────────────
  baseline = list(
    bl_callahan, bl_peek, bl_clayton, bl_singer,
    bl_underwood, bl_rosa, bl_doe, bl_kint
  ),

  statistics = list(
    st_callahan, st_peek, st_clayton, st_singer,
    st_underwood, st_rosa, st_doe, st_kint
  ),

  # ── Primary outcome (OR for delirium incidence) ─────────────────────────
  outcome_estimate = c(0.65, 0.82, 0.55, 0.95, 0.48, 0.30, 0.21, 0.58),
  outcome_ci_lower = c(0.38, 0.68, 0.31, 0.77, 0.26, 0.15, 0.11, 0.32),
  outcome_ci_upper = c(1.11, 0.99, 0.98, 1.17, 0.89, 0.60, 0.40, 1.05),
  outcome_se       = c(0.274, 0.095, 0.295, 0.108, 0.314, 0.355, 0.330, 0.304),
  outcome_log_scale = rep(TRUE, 8)
)

use_data(inspect_sr_example, overwrite = TRUE)
