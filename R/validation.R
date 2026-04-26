# ---- Core validation helpers ----

#' Check that a value is a bayesma_prior object (or NULL)
#'
#' Used by validate_prior_args(), validate_multi_arm_args(),
#' and validate_robust_args().
#'
#' @param x Object to check
#' @param name Argument name for error messages
#' @noRd
check_is_prior <- function(x, name, call = rlang::caller_env()) {
  if (!is.null(x) && !inherits(x, "bayesma_prior")) {
    cli::cli_abort(c(
      "{.arg {name}} must be a prior object (e.g. {.code normal(0, 10)}).",
      "i" = "Use {.fn normal}, {.fn half_normal}, {.fn half_cauchy}, etc."
    ), call = call)
  }
}


# ---- bayesma() input validation ----

#' Validate argument combinations for bayesma()
#'
#' Checks that model_type, stage, re_dist, and use_known_bias form a
#' valid combination. Aborts on invalid combinations.
#'
#' @noRd
validate_bayesma_args <- function(model_type, stage, re_dist,
                                  use_known_bias, data,
                                  call = rlang::caller_env()) {
  bias_models <- c("bias_corrected", "bc_bnp", "selection_copas",
                    "selection_weight", "pet_peese", "mixture_model")

  if (model_type == "common_effect" && re_dist != "normal") {
    cli::cli_abort(c(
      "Non-normal RE distributions require {.val random_effect} model.",
      "i" = "Set {.arg model_type = 'random_effect'} or {.arg re_dist = 'normal'}."
    ), call = call)
  }

  if (model_type %in% bias_models && stage != "two_stage") {
    cli::cli_abort(c(
      "{.val {model_type}} model is only available for two-stage models.",
      "i" = "Set {.arg stage = 'two_stage'}."
    ), call = call)
  }

  if (use_known_bias && !("biased" %in% names(data))) {
    cli::cli_abort(c(
      "{.arg use_known_bias = TRUE} requires a column named {.val biased} in {.arg data}.",
      "i" = "This column should contain 0/1 indicating known bias status."
    ), call = call)
  }

  invisible(NULL)
}


#' Reconcile conflicting bayesma() options with warnings
#'
#' Handles cases where user-supplied options are incompatible: emits a
#' warning and returns corrected values.
#'
#' @return Named list with corrected \code{re_dist} and \code{small_sample}.
#' @noRd
reconcile_bayesma_options <- function(model_type, stage, re_dist,
                                      small_sample) {
  bias_models <- c("bias_corrected", "bc_bnp", "selection_copas",
                    "selection_weight", "pet_peese", "mixture_model")

  rbc_compatible <- c("normal", "t")

  if (model_type == "selection_copas" && !(re_dist %in% rbc_compatible)) {
    cli::cli_warn(c(
      "{.val selection_copas} (RBC) supports {.val normal} or {.val t} \\
       random-effects distributions only.",
      "i" = "Ignoring {.arg re_dist = '{re_dist}'}; using {.val t}."
    ))
    re_dist <- "t"
  } else if (model_type %in% setdiff(bias_models,
             c("selection_copas", "bc_bnp", "mixture_model")) &&
             re_dist != "normal") {
    cli::cli_warn(c(
      "{.val {model_type}} model uses its own RE structure.",
      "i" = "Ignoring {.arg re_dist = '{re_dist}'}."
    ))
    re_dist <- "normal"
  }

  if (stage == "one_stage" && small_sample != "none") {
    cli::cli_warn("Small-sample adjustments are two-stage only. Ignoring.")
    small_sample <- "none"
  }

  list(re_dist = re_dist, small_sample = small_sample)
}


#' Validate multi-arm arguments for bayesma()
#'
#' @noRd
validate_multi_arm_args <- function(multi_arm, stage, re_dist,
                                    rho_prior, data,
                                    call = rlang::caller_env()) {
  if (is.null(multi_arm)) return(invisible(NULL))

  if (stage != "one_stage") {
    cli::cli_abort(c(
      "{.arg multi_arm} is only supported for {.val one_stage} models.",
      "i" = "Set {.arg stage = 'one_stage'} to use multi-arm study modelling."
    ), call = call)
  }

  if (!(multi_arm %in% names(data))) {
    cli::cli_abort(c(
      "Column {.val {multi_arm}} not found in {.arg data}.",
      "i" = "Specify the column name that identifies multi-arm study grouping."
    ), call = call)
  }

  check_is_prior(rho_prior, "rho_prior", call = call)

  if (re_dist == "mixture") {
    cli::cli_abort(c(
      "Multi-arm studies are not yet supported with {.arg re_dist = 'mixture'}.",
      "i" = "Use {.arg re_dist = 'normal'}, {.arg re_dist = 't'}, or {.arg re_dist = 'skew_normal'}."
    ), call = call)
  }

  invisible(NULL)
}


#' Validate robust modelling arguments for bayesma()
#'
#' @noRd
validate_robust_args <- function(robust, robust_prior, robust_weight,
                                 model_type,
                                 call = rlang::caller_env()) {
  if (!robust) return(invisible(NULL))

  bias_models <- c("bias_corrected", "bc_bnp", "selection_copas",
                    "selection_weight", "pet_peese", "mixture_model")

  if (model_type %in% bias_models) {
    cli::cli_warn(c(
      "{.arg robust = TRUE} with bias/selection models may lead to identifiability issues.",
      "i" = "Both mechanisms try to explain outlying studies."
    ))
  }

  check_is_prior(robust_prior, "robust_prior", call = call)
  check_is_prior(robust_weight, "robust_weight", call = call)

  invisible(NULL)
}


#' Validate that required data columns are supplied for the likelihood
#'
#' @noRd
validate_required_columns <- function(likelihood, event_ctrl, event_int,
                                      mean_ctrl, mean_int, sd_ctrl, sd_int,
                                      call = rlang::caller_env()) {
  if (likelihood %in% c("binomial", "poisson")) {
    if (is.null(event_ctrl) || is.null(event_int)) {
      cli::cli_abort(
        "{.val {likelihood}} likelihood requires {.arg event_ctrl} and {.arg event_int}.",
        call = call
      )
    }
  }

  if (likelihood == "gaussian") {
    missing_args <- c(
      if (is.null(mean_ctrl)) "mean_ctrl", if (is.null(mean_int)) "mean_int",
      if (is.null(sd_ctrl)) "sd_ctrl", if (is.null(sd_int)) "sd_int"
    )
    if (length(missing_args) > 0) {
      cli::cli_abort(c(
        "Gaussian likelihood requires all of {.arg mean_ctrl}, {.arg mean_int}, {.arg sd_ctrl}, {.arg sd_int}.",
        "x" = "Missing: {.arg {missing_args}}."
      ), call = call)
    }
  }

  invisible(NULL)
}


#' Validate prior arguments for bayesma()
#'
#' Checks that user-supplied priors are bayesma_prior objects and are
#' appropriate for the chosen model_type/re_dist/stage combination.
#'
#' @noRd
validate_prior_args <- function(stage, model_type, re_dist,
                                mu_prior, tau_prior, gamma_prior,
                                nu_prior, alpha_prior, mixture_priors,
                                b_prior, p_bias_prior, w_bias_prior,
                                use_known_bias, selection_priors,
                                call = rlang::caller_env()) {

  is_re <- model_type == "random_effect"
  is_bias <- model_type == "bias_corrected"
  is_bc_bnp <- model_type == "bc_bnp"
  is_selection <- model_type %in% c("selection_copas", "selection_weight", "pet_peese")

  check_is_prior(mu_prior, "mu_prior", call = call)
  check_is_prior(tau_prior, "tau_prior", call = call)
  check_is_prior(gamma_prior, "gamma_prior", call = call)
  check_is_prior(nu_prior, "nu_prior", call = call)
  check_is_prior(alpha_prior, "alpha_prior", call = call)
  check_is_prior(b_prior, "b_prior", call = call)
  check_is_prior(p_bias_prior, "p_bias_prior", call = call)
  check_is_prior(w_bias_prior, "w_bias_prior", call = call)

  if (!is.null(mixture_priors)) {
    if (!is.list(mixture_priors)) {
      cli::cli_abort("{.arg mixture_priors} must be a named list.", call = call)
    }
    purrr::walk2(mixture_priors, names(mixture_priors), function(p, nm) {
      check_is_prior(p, paste0("mixture_priors$", nm), call = call)
    })
  }

  # --- Context checks: standard model priors ---
  if (!is.null(gamma_prior) && stage == "two_stage") {
    cli::cli_warn("{.arg gamma_prior} is only used in one-stage models. Ignoring.")
  }

  if (!is.null(tau_prior) && model_type == "common_effect") {
    cli::cli_warn("{.arg tau_prior} is only used in random-effects or bias-corrected models. Ignoring.")
  }

  if (!is.null(nu_prior) && re_dist != "t") {
    cli::cli_abort(c(
      "{.arg nu_prior} is only used when {.arg re_dist = 't'}.",
      "i" = "Current {.arg re_dist} is {.val {re_dist}}."
    ), call = call)
  }

  if (!is.null(alpha_prior) && re_dist != "skew_normal") {
    cli::cli_abort(c(
      "{.arg alpha_prior} is only used when {.arg re_dist = 'skew_normal'}.",
      "i" = "Current {.arg re_dist} is {.val {re_dist}}."
    ), call = call)
  }

  uses_mixture <- (re_dist == "mixture") || (model_type == "mixture_model")

  if (!is.null(mixture_priors) && !uses_mixture) {
    cli::cli_abort(c(
      "{.arg mixture_priors} is only used when {.arg re_dist = 'mixture'} \\
       or {.arg model_type = 'mixture_model'}.",
      "i" = "Current: {.arg re_dist} = {.val {re_dist}}, \\
             {.arg model_type} = {.val {model_type}}."
    ), call = call)
  }

  if (uses_mixture && !is.null(mixture_priors)) {
    allowed <- if (model_type == "mixture_model") c("w", "mu_k", "tau_k")
    else if (stage == "two_stage") c("w", "mu_k", "tau_k")
    else c("w", "delta_k", "tau_k")
    unknown <- setdiff(names(mixture_priors), allowed)
    if (length(unknown) > 0) {
      cli::cli_abort(c(
        "Unknown keys in {.arg mixture_priors}: {.val {unknown}}.",
        "i" = "Allowed keys: {.val {allowed}}."
      ), call = call)
    }
  }

  # --- Context checks: bias-corrected priors ---
  bc_priors_supplied <- !is.null(b_prior) || !is.null(w_bias_prior)
  if (bc_priors_supplied && !is_bias) {
    cli::cli_abort(c(
      "{.arg b_prior} and {.arg w_bias_prior} are only used when \\
       {.arg model_type = 'bias_corrected'}.",
      "i" = "Current {.arg model_type} is {.val {model_type}}."
    ), call = call)
  }
  if (!is.null(p_bias_prior) && !is_bias && !is_bc_bnp) {
    cli::cli_abort(c(
      "{.arg p_bias_prior} is only used when {.arg model_type} is \\
       {.val bias_corrected} or {.val bc_bnp}.",
      "i" = "Current {.arg model_type} is {.val {model_type}}."
    ), call = call)
  }

  if (use_known_bias && !is_bias) {
    cli::cli_abort(c(
      "{.arg use_known_bias} is only used when {.arg model_type = 'bias_corrected'}.",
      "i" = "Current {.arg model_type} is {.val {model_type}}."
    ), call = call)
  }

  # --- Context checks: selection priors ---
  if (!is.null(selection_priors) && !is_selection && !is_bias) {
    cli::cli_abort(c(
      "{.arg selection_priors} is only used with selection/bias model types.",
      "i" = "Current {.arg model_type} is {.val {model_type}}."
    ), call = call)
  }

  if (!is.null(selection_priors)) {
    if (!is.list(selection_priors)) {
      cli::cli_abort("{.arg selection_priors} must be a named list.", call = call)
    }
    purrr::walk2(selection_priors, names(selection_priors), function(p, nm) {
      check_is_prior(p, paste0("selection_priors$", nm), call = call)
    })
    allowed_keys <- switch(model_type,
                           selection_copas  = c("gamma0", "gamma1", "rho"),
                           selection_weight = c("omega"),
                           pet_peese        = c("beta_bias"),
                           NULL
    )
    if (!is.null(allowed_keys)) {
      unknown <- setdiff(names(selection_priors), allowed_keys)
      if (length(unknown) > 0) {
        cli::cli_abort(c(
          "Unknown keys in {.arg selection_priors}: {.val {unknown}}.",
          "i" = "Allowed keys for {.val {model_type}}: {.val {allowed_keys}}."
        ), call = call)
      }
    }
  }

  # --- Bounded parameter warnings ---
  bounded_check <- function(prior, name) {
    if (!is.null(prior) && prior$family == "normal") {
      cli::cli_warn(c(
        "{.arg {name}} is lower-bounded at 0 but prior is {.fn normal}.",
        "i" = "Consider {.fn half_normal}, {.fn half_cauchy}, or {.fn half_student_t}."
      ))
    }
  }
  if (is_re || is_bias || is_selection) bounded_check(tau_prior, "tau_prior")

  invisible(NULL)
}


# ---- Plot/post-hoc validation ----

#' Validate inputs for bayes_forest (bayesma pathway)
#'
#' @noRd
validate_inputs_bayesma <- function(
    model,
    data,
    measure,
    studyvar,
    year,
    subgroup,
    subgroup_var,
    sort_studies_by,
    shrinkage_output
) {
  if (!inherits(model, "bayesma")) {
    cli::cli_abort("{.arg model} must be a {.cls bayesma} object.")
  }
  if (!is.data.frame(data)) {
    cli::cli_abort("{.arg data} must be a data frame.")
  }
  rlang::arg_match(measure, c("OR", "HR", "RR", "IRR", "MD", "SMD"))
  rlang::arg_match(sort_studies_by, c("author", "year", "effect"))
  rlang::arg_match(shrinkage_output, c("density", "pointinterval"))

  # studyvar is required for column renaming
  if (rlang::quo_is_null(rlang::enquo(studyvar))) {
    cli::cli_abort("{.arg studyvar} must be provided (bare column name).")
  }

  studyvar_name <- rlang::as_name(rlang::ensym(studyvar))
  if (!studyvar_name %in% names(data)) {
    cli::cli_abort("Can't find column {.val {studyvar_name}} in {.arg data}.")
  }

  if (!rlang::quo_is_null(rlang::enquo(year))) {
    year_name <- rlang::as_name(rlang::ensym(year))
    if (!year_name %in% names(data)) {
      cli::cli_abort("Can't find column {.val {year_name}} in {.arg data}.")
    }
  }

  # Check that bayesma object has required components
  if (is.null(model$meta)) {
    cli::cli_abort("bayesma object is missing {.code meta} component.")
  }
  if (is.null(model$draws) && is.null(model$fit)) {
    cli::cli_abort("bayesma object must contain either {.code draws} or {.code fit}.")
  }

  # Validate subgroup variable if subgroup = TRUE
  if (isTRUE(subgroup)) {
    if (is.null(model$meta$call_args)) {
      cli::cli_abort(
        c(
          "Subgroup analysis with bayesma requires stored call arguments for refitting.",
          "i" = "Refit your model with the latest version of bayesma which stores call_args."
        )
      )
    }

    if (rlang::quo_is_null(rlang::enquo(subgroup_var))) {
      if (!"Subgroup" %in% names(data)) {
        cli::cli_abort("{.arg subgroup = TRUE} but no {.arg subgroup_var} specified and no {.val Subgroup} column found in {.arg data}.")
      }
    } else {
      subgroup_var_name <- rlang::as_name(rlang::ensym(subgroup_var))
      if (!subgroup_var_name %in% names(data)) {
        cli::cli_abort("Can't find subgroup column {.val {subgroup_var_name}} in {.arg data}.")
      }
    }
  }

  invisible(TRUE)
}


#' Validate Inputs for Sensitivity Plot
#'
#' @description
#' Internal validation function that checks the validity of core inputs to the
#' sensitivity plot function. Accepts bayesma objects.
#'
#' @param model A bayesma object to validate.
#' @param data A data frame to validate.
#' @param priors Priors object(s) — bayesma priors for bayesma.
#' @param measure Character string specifying the effect measure.
#'
#' @return Invisible TRUE if validation passes, otherwise throws an error.
#'
#' @noRd
validate_inputs_sens_plot <- function(
    model,
    data,
    priors,
    measure
) {

  if (!inherits(model, "bayesma")) {
    cli::cli_abort("{.arg model} must be a {.cls bayesma} object.")
  }

  if (!is.data.frame(data)) {
    cli::cli_abort("{.arg data} must be a data frame.")
  }

  rlang::arg_match(measure, c("OR", "HR", "RR", "IRR", "MD", "SMD"))

  validate_priors_bayesma(priors)

  invisible(TRUE)
}

#' Validate Priors Object
#'
#' For bayesma models, priors should be a named list where each element is
#' itself a named list of bayesma_prior objects (e.g., list(mu_prior = normal(0, 10),
#' tau_prior = half_cauchy(0, 0.5))).
#'
#' @noRd
validate_priors_bayesma <- function(priors) {

  if (is.null(priors)) {
    return(invisible(TRUE))
  }

  if (!is.list(priors)) {
    cli::cli_abort(
      "{.arg priors} must be a named list of prior specifications for bayesma models."
    )
  }

  # Each element should be a named list of bayesma_prior objects,
  # or a single bayesma_prior
  purrr::walk(priors, function(p) {
    if (inherits(p, "bayesma_prior")) return(invisible(NULL))
    if (is.list(p)) {
      # Allow metadata fields (name, label, etc.) alongside bayesma_prior objects
      known_metadata <- c("name", "label", "description")
      non_prior <- purrr::keep(p, ~ !inherits(.x, "bayesma_prior") &&
                                 !is.null(.x) &&
                                 !is.character(.x))
      # Also remove known metadata keys
      non_prior <- non_prior[!names(non_prior) %in% known_metadata]
      if (length(non_prior) > 0) {
        cli::cli_abort(
          "Each element in {.arg priors} must be a {.cls bayesma_prior} object or a named list of {.cls bayesma_prior} objects."
        )
      }
    }
  })

  invisible(TRUE)
}

#' Validate sensitivity analysis priors
#'
#' @noRd
validate_sensitivity_priors <- function(priors, call = rlang::caller_env()) {
  if (!is.list(priors) || length(priors) == 0) {
    cli::cli_abort("{.arg priors} must be a non-empty list.", call = call)
  }

  for (nm in names(priors)) {
    p <- priors[[nm]]
    if (!is.list(p)) cli::cli_abort("priors${nm} must be a list.", call = call)

    if (is.null(p$mu_prior)) {
      cli::cli_abort("priors${nm} must contain {.code mu_prior}.", call = call)
    }
    # tau can be optional for CE models, but required if you run RE models
    # (keep it permissive here)
    if (!is.null(p$name) && !is.character(p$name)) {
      cli::cli_abort("priors${nm}$name must be a character string.", call = call)
    }
  }

  priors
}
