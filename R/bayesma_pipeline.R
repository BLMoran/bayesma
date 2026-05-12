# bayesma modular pipeline
#
# Five user-inspectable stages:
#   1. bayesma_spec()       -- validate + extract -> bayesma_spec
#   2. bayesma_stan_code()  -- spec -> named list of Stan blocks + full program
#   3. bayesma_stan_data()  -- spec -> cmdstanr data list
#   4. bayesma_fit()        -- compile + sample
#   5. bayesma_extract()    -- fit + spec -> tidy effect components
#      bayesma_output()     -- assemble final `bayesma` object
#
# `bayesma()` is a thin orchestrator over these stages. Each stage is exported
# so users can pause at any point to inspect, customise, or plug in their own
# Stan program via the `custom_model` argument.


# -----------------------------------------------------------------------------
# Stage 0: Stan block parser (internal)
# -----------------------------------------------------------------------------

#' Split a Stan program into its named blocks
#'
#' Brace-balanced parser that extracts `functions`, `data`, `transformed data`,
#' `parameters`, `transformed parameters`, `model`, and `generated quantities`
#' blocks from a Stan program. Missing blocks return `""`.
#'
#' The returned content excludes the surrounding braces; use
#' `assemble_stan_blocks()` to reconstitute a full program.
#'
#' @param stan_code Character. A Stan program.
#' @return Named list of block bodies (character scalars).
#' @noRd
parse_stan_blocks <- function(stan_code) {
  blocks <- list(
    functions              = "",
    data                   = "",
    transformed_data       = "",
    parameters             = "",
    transformed_parameters = "",
    model                  = "",
    generated_quantities   = ""
  )

  clean <- gsub("//[^\n]*",                      "", stan_code)
  clean <- gsub("/\\*([^*]|\\*(?!/))*\\*/",      "", clean, perl = TRUE)

  headers <- c(
    transformed_parameters = "transformed\\s+parameters",
    transformed_data       = "transformed\\s+data",
    generated_quantities   = "generated\\s+quantities",
    functions              = "functions",
    parameters             = "parameters",
    model                  = "model",
    data                   = "data"
  )

  n   <- nchar(clean)
  pos <- 1L

  while (pos <= n) {
    best <- NULL
    for (key in names(headers)) {
      pattern <- paste0("(^|[\\s;}])", headers[[key]], "\\s*\\{")
      rest <- substr(clean, pos, n)
      m <- regexpr(pattern, rest, perl = TRUE)
      if (m != -1L) {
        if (is.null(best) || m < best$relpos) {
          best <- list(relpos = as.integer(m),
                       length = attr(m, "match.length"),
                       key    = key)
        }
      }
    }
    if (is.null(best)) break

    brace_open <- pos + best$relpos - 1L + best$length - 1L
    depth <- 1L
    i <- brace_open + 1L
    while (i <= n && depth > 0L) {
      ch <- substr(clean, i, i)
      if (ch == "{")      depth <- depth + 1L
      else if (ch == "}") depth <- depth - 1L
      i <- i + 1L
    }
    brace_close <- i - 1L
    body <- substr(clean, brace_open + 1L, brace_close - 1L)
    blocks[[best$key]] <- trimws(body, which = "both")
    pos <- brace_close + 1L
  }

  blocks
}


#' Reassemble a Stan program from a named list of block bodies
#'
#' @param blocks Named list (output of `parse_stan_blocks()` or a user-built
#'   equivalent). Block bodies must not include surrounding braces.
#' @return Character scalar: a complete Stan program.
#' @noRd
assemble_stan_blocks <- function(blocks) {
  order_keys <- c("functions", "data", "transformed_data", "parameters",
                  "transformed_parameters", "model", "generated_quantities")
  display <- c(
    functions              = "functions",
    data                   = "data",
    transformed_data       = "transformed data",
    parameters             = "parameters",
    transformed_parameters = "transformed parameters",
    model                  = "model",
    generated_quantities   = "generated quantities"
  )
  parts <- purrr::map_chr(order_keys, function(k) {
    body <- blocks[[k]] %||% ""
    if (!nzchar(trimws(body))) return("")
    paste0(display[[k]], " {\n", body, "\n}")
  })
  paste(parts[nzchar(parts)], collapse = "\n\n")
}


# -----------------------------------------------------------------------------
# Stage 1: spec -- validate args, extract data, resolve priors
# -----------------------------------------------------------------------------

#' @noRd
bayesma_spec <- function(
    data,
    studyvar,
    event_ctrl       = NULL,
    event_int        = NULL,
    mean_ctrl        = NULL,
    mean_int         = NULL,
    sd_ctrl          = NULL,
    sd_int           = NULL,
    n_ctrl           = NULL,
    n_int            = NULL,
    likelihood       = c("binomial", "gaussian", "poisson"),
    model_type       = c("random_effect", "common_effect", "bias_corrected",
                         "bc_bnp", "selection_copas", "selection_weight",
                         "pet_peese", "mixture_model"),
    stage            = c("one_stage", "two_stage"),
    re_dist          = c("normal", "t", "skew_normal", "mixture"),
    small_sample     = c("none", "t_approx", "hjsk"),
    multi_arm        = NULL,
    rho_prior        = NULL,
    mu_prior         = NULL,
    tau_prior        = NULL,
    gamma_prior      = NULL,
    nu_prior         = NULL,
    alpha_prior      = NULL,
    mixture_priors   = NULL,
    b_prior          = NULL,
    p_bias_prior     = NULL,
    w_bias_prior     = NULL,
    mu_beta_prior    = NULL,
    tau_beta_prior   = NULL,
    bnp_concentration_max = NULL,
    bnp_K_max        = NULL,
    use_known_bias   = FALSE,
    selection_priors = NULL,
    p_cutoffs        = c(0.025, 0.05),
    n_components     = 2L,
    robust           = FALSE,
    robust_prior     = NULL,
    robust_df        = 4,
    robust_weight    = NULL,
    custom_model     = NULL,
    custom_data      = NULL,
    estimand         = NULL,
    cate_covariate   = NULL,
    baseline_risk    = NULL,
    re_min_k         = NULL
) {
  likelihood   <- rlang::arg_match(likelihood)
  model_type   <- rlang::arg_match(model_type)
  stage        <- rlang::arg_match(stage)
  re_dist      <- rlang::arg_match(re_dist)
  small_sample <- rlang::arg_match(small_sample)

  estimand <- resolve_estimand(estimand, likelihood)
  validate_estimand_args(estimand, likelihood, cate_covariate, baseline_risk,
                         data, stage)

  if (!is.null(custom_model) && !is.character(custom_model)) {
    cli::cli_abort(
      "{.arg custom_model} must be a character scalar containing Stan code.",
      call = rlang::caller_env()
    )
  }
  if (!is.null(custom_data) && !is.list(custom_data)) {
    cli::cli_abort(
      "{.arg custom_data} must be a named list.",
      call = rlang::caller_env()
    )
  }

  validate_bayesma_args(model_type, stage, re_dist, use_known_bias, data)

  opts         <- reconcile_bayesma_options(model_type, stage, re_dist, small_sample)
  re_dist      <- opts$re_dist
  small_sample <- opts$small_sample

  has_multi_arm <- !is.null(multi_arm)
  validate_multi_arm_args(multi_arm, stage, re_dist, rho_prior, data)
  if (has_multi_arm && is.null(rho_prior)) rho_prior <- uniform(-1, 1)

  validate_robust_args(robust, robust_prior, robust_weight, model_type)
  if (robust) {
    if (is.null(robust_prior)) {
      robust_prior <- if (stage == "one_stage" && likelihood == "gaussian") {
        normal(0, 100)
      } else {
        normal(0, 10)
      }
    }
    if (is.null(robust_weight)) robust_weight <- beta(2, 2)
  }

  validate_prior_args(
    stage, model_type, re_dist,
    mu_prior, tau_prior, gamma_prior,
    nu_prior, alpha_prior, mixture_priors,
    b_prior, p_bias_prior, w_bias_prior,
    use_known_bias, selection_priors
  )

  priors <- resolve_priors(
    stage, likelihood, model_type, re_dist,
    mu_prior, tau_prior, gamma_prior,
    nu_prior, alpha_prior, mixture_priors,
    b_prior, p_bias_prior, w_bias_prior,
    selection_priors,
    mu_beta_prior  = mu_beta_prior,
    tau_beta_prior = tau_beta_prior
  )
  if (has_multi_arm) priors$rho <- rho_prior

  validate_required_columns(likelihood, event_ctrl, event_int,
                            mean_ctrl, mean_int, sd_ctrl, sd_int)

  extract_col <- function(d, var_name) {
    if (is.null(var_name)) return(NULL)
    val <- d[[var_name]]
    if (is.null(val)) {
      cli::cli_abort("Variable {.val {var_name}} not found in {.arg data}.",
                     call = rlang::caller_env())
    }
    val
  }

  study_vec    <- extract_col(data, studyvar)
  S            <- length(study_vec)
  study_labels <- if (is.factor(study_vec)) levels(study_vec)
                  else as.character(study_vec)

  n_c <- extract_col(data, n_ctrl)
  n_i <- extract_col(data, n_int)

  if (likelihood %in% c("binomial", "poisson")) {
    outcome_ctrl <- extract_col(data, event_ctrl)
    outcome_int  <- extract_col(data, event_int)
    sd_c <- NULL; sd_i <- NULL
  } else {
    outcome_ctrl <- extract_col(data, mean_ctrl)
    outcome_int  <- extract_col(data, mean_int)
    sd_c <- extract_col(data, sd_ctrl)
    sd_i <- extract_col(data, sd_int)
  }

  es <- compute_effect_sizes(outcome_ctrl, outcome_int, n_c, n_i,
                             sd_c, sd_i, S, likelihood)

  robust_config <- if (robust) {
    list(enabled = TRUE, prior = robust_prior,
         df = robust_df, weight = robust_weight)
  } else {
    list(enabled = FALSE)
  }

  multi_arm_config <- if (has_multi_arm) {
    ma_vec    <- extract_col(data, multi_arm)
    ma_ids    <- as.integer(as.factor(ma_vec))
    list(
      enabled        = TRUE,
      ma_study_id    = ma_ids,
      n_ma_studies   = length(unique(ma_ids)),
      arms_per_study = as.integer(table(ma_ids)),
      rho_prior      = rho_prior
    )
  } else {
    list(enabled = FALSE)
  }

  known_bias_vec <- if (model_type == "bias_corrected" && use_known_bias) {
    as.integer(data[["biased"]])
  } else if (model_type == "bias_corrected") {
    rep(0L, S)
  } else {
    NULL
  }

  bc_bnp_config <- if (model_type == "bc_bnp") {
    a0 <- priors$p_bias$alpha
    a1 <- priors$p_bias$beta
    alpha_max_default <- (1 / 5) * ((S - 1) * a0 - a1) / (a0 + a1)
    if (!is.finite(alpha_max_default) || alpha_max_default <= 0.5) {
      alpha_max_default <- 1
    }
    am <- bnp_concentration_max %||% alpha_max_default
    Km <- bnp_K_max %||% (1L + as.integer(ceiling(5 * am)))
    list(alpha_max = am, K_max = Km, a0 = a0, a1 = a1)
  } else {
    NULL
  }

  call_args <- list(
    studyvar         = studyvar,
    event_ctrl       = event_ctrl,      event_int        = event_int,
    mean_ctrl        = mean_ctrl,       mean_int         = mean_int,
    sd_ctrl          = sd_ctrl,         sd_int           = sd_int,
    n_ctrl           = n_ctrl,          n_int            = n_int,
    likelihood       = likelihood,      model_type       = model_type,
    stage            = stage,           re_dist          = re_dist,
    small_sample     = small_sample,    multi_arm        = multi_arm,
    rho_prior        = rho_prior,       mu_prior         = mu_prior,
    tau_prior        = tau_prior,       gamma_prior      = gamma_prior,
    nu_prior         = nu_prior,        alpha_prior      = alpha_prior,
    mixture_priors   = mixture_priors,  b_prior          = b_prior,
    p_bias_prior     = p_bias_prior,    w_bias_prior     = w_bias_prior,
    use_known_bias   = use_known_bias,  selection_priors = selection_priors,
    p_cutoffs        = p_cutoffs,       n_components     = n_components,
    robust           = robust,          robust_prior     = robust_prior,
    robust_df        = robust_df,       robust_weight    = robust_weight,
    custom_model     = custom_model,    custom_data      = custom_data,
    estimand         = estimand,
    cate_covariate   = cate_covariate,  baseline_risk    = baseline_risk,
    re_min_k         = re_min_k
  )

  spec <- list(
    likelihood       = likelihood,
    model_type       = model_type,
    stage            = stage,
    re_dist          = re_dist,
    small_sample     = small_sample,
    study_vec        = study_vec,
    study_labels     = study_labels,
    S                = S,
    outcome_ctrl     = outcome_ctrl,
    outcome_int      = outcome_int,
    n_c              = n_c,
    n_i              = n_i,
    sd_c             = sd_c,
    sd_i             = sd_i,
    es               = es,
    priors           = priors,
    robust_config    = robust_config,
    multi_arm_config = multi_arm_config,
    has_multi_arm    = has_multi_arm,
    p_cutoffs        = p_cutoffs,
    n_components     = as.integer(n_components),
    use_known_bias   = use_known_bias,
    known_bias_vec   = known_bias_vec,
    bc_bnp_config    = bc_bnp_config,
    custom_model     = custom_model,
    custom_data      = custom_data,
    estimand         = estimand,
    cate_covariate   = cate_covariate,
    baseline_risk    = baseline_risk,
    call_args        = call_args
  )
  class(spec) <- c("bayesma_spec", "list")
  spec
}


#' @export
print.bayesma_spec <- function(x, ...) {
  cat("<bayesma_spec>\n",
      "  likelihood   : ", x$likelihood, "\n",
      "  model_type   : ", x$model_type, "\n",
      "  stage        : ", x$stage, "\n",
      "  re_dist      : ", x$re_dist, "\n",
      "  small_sample : ", x$small_sample, "\n",
      "  studies (S)  : ", x$S, "\n",
      "  robust       : ", isTRUE(x$robust_config$enabled), "\n",
      "  multi_arm    : ", isTRUE(x$multi_arm_config$enabled), "\n",
      "  custom_model : ", !is.null(x$custom_model), "\n",
      sep = "")
  invisible(x)
}


# -----------------------------------------------------------------------------
# Internal: prior-only Stan code injection
# -----------------------------------------------------------------------------

#' @noRd
inject_prior_only <- function(stan_code) {
  # 1. Add `int<lower=0, upper=1> prior_only;` to the data block.
  #    Data blocks contain no braces in their body, so [^}]* is safe.
  stan_code <- sub(
    pattern     = "(data \\{[^}]*?)(\n\\})",
    replacement = "\\1\n  int<lower=0, upper=1> prior_only;\n\\2",
    x           = stan_code,
    perl        = TRUE
  )

  # 2. Wrap the likelihood section of the model block in if (!prior_only).
  #    The generators always emit "  // Likelihood\n" as a separator; we
  #    split there, then re-join with the conditional guard.
  marker <- "// Likelihood\n"
  if (!grepl(marker, stan_code, fixed = TRUE)) return(stan_code)

  # Find the marker (allowing for leading whitespace after auto-format)
  split_pos <- regexpr("\\s*// Likelihood\n", stan_code, perl = TRUE)
  if (split_pos == -1L) return(stan_code)

  before_lik <- substr(stan_code, 1L, split_pos - 1L)
  after_lik  <- substr(stan_code, split_pos + attr(split_pos, "match.length"), nchar(stan_code))

  # The model block's closing `}` is the first `\n}` (no leading spaces).
  end_pos  <- regexpr("\n\\}", after_lik)
  if (end_pos == -1L) return(stan_code)

  lik_body   <- substr(after_lik, 1L, end_pos - 1L)
  after_body <- substr(after_lik, end_pos, nchar(after_lik))

  paste0(
    before_lik,
    "  if (!prior_only) {\n  // Likelihood\n",
    lik_body, "\n  }",
    after_body
  )
}


# -----------------------------------------------------------------------------
# Stage 2: bayesma_stan_code -- spec -> named blocks + full program
# -----------------------------------------------------------------------------

#' @noRd
bayesma_stan_code <- function(spec, format = TRUE) {
  if (!inherits(spec, "bayesma_spec")) {
    cli::cli_abort("{.arg spec} must be a {.cls bayesma_spec} object.")
  }

  raw <- if (!is.null(spec$custom_model)) {
    spec$custom_model
  } else {
    generate_full_stan_program(spec)
  }

  full <- if (isTRUE(format)) format_stan_code(raw) else raw
  blocks <- parse_stan_blocks(full)

  out <- c(blocks, list(full = full))
  class(out) <- c("bayesma_stan_code", "list")
  out
}


#' Pretty-print Stan code with stanc --auto-format
#'
#' Thin wrapper around `cmdstanr::CmdStanModel$format()` that takes a Stan
#' program string and returns the canonically formatted version. Returns the
#' input unchanged when cmdstan is unavailable or the program does not parse.
#'
#' @param stan_code Character scalar. A Stan program.
#' @return Character scalar.
#' @noRd
format_stan_code <- function(stan_code) {
  tryCatch({
    tmp <- cmdstanr::write_stan_file(stan_code)
    on.exit(unlink(tmp), add = TRUE)
    mod <- cmdstanr::cmdstan_model(tmp, compile = FALSE)
    out <- utils::capture.output(
      mod$format(overwrite_file = FALSE, canonicalize = "deprecations")
    )
    paste(out, collapse = "\n")
  }, error = function(e) stan_code)
}


#' @noRd
generate_full_stan_program <- function(spec) {
  mt <- spec$model_type

  if (mt %in% c("random_effect", "common_effect")) {
    if (spec$stage == "two_stage") {
      use_t <- spec$small_sample %in% c("t_approx", "hjsk")
      generate_stan_code_two_stage(
        model_type       = mt,
        re_dist          = spec$re_dist,
        use_t_likelihood = use_t,
        priors           = spec$priors,
        robust_config    = spec$robust_config
      )
    } else {
      generate_stan_code_one_stage(
        likelihood       = spec$likelihood,
        model_type       = mt,
        re_dist          = spec$re_dist,
        priors           = spec$priors,
        robust_config    = spec$robust_config,
        multi_arm_config = spec$multi_arm_config
      )
    }
  } else if (mt == "bias_corrected") {
    generate_stan_code_bias_corrected(spec$priors)
  } else if (mt == "bc_bnp") {
    generate_stan_code_bc_bnp(spec$priors, K = spec$bc_bnp_config$K_max)
  } else if (mt == "mixture_model") {
    generate_stan_code_mixture_model(spec$priors, M = spec$n_components)
  } else if (mt == "selection_copas") {
    generate_stan_code_selection_copas(spec$priors, re_dist = spec$re_dist)
  } else if (mt == "selection_weight") {
    generate_stan_code_selection_weight(spec$priors, length(spec$p_cutoffs) + 1L)
  } else if (mt == "pet_peese") {
    generate_stan_code_pet_peese(spec$priors, predictor = "inv_sqrt_n")
  } else {
    cli::cli_abort("Unknown {.arg model_type}: {.val {mt}}")
  }
}


#' @export
print.bayesma_stan_code <- function(x, ...) {
  cat(x$full, "\n", sep = "")
  invisible(x)
}

#' @export
format.bayesma_stan_code <- function(x, ...) x$full


# -----------------------------------------------------------------------------
# Stage 3: bayesma_stan_data -- spec -> cmdstanr data list
# -----------------------------------------------------------------------------

#' @noRd
bayesma_stan_data <- function(spec) {
  if (!inherits(spec, "bayesma_spec")) {
    cli::cli_abort("{.arg spec} must be a {.cls bayesma_spec} object.")
  }

  sd_list <- switch(spec$model_type,
    random_effect    = build_stan_data_standard(spec),
    common_effect    = build_stan_data_standard(spec),
    bias_corrected   = build_stan_data_bias_corrected(spec),
    bc_bnp           = build_stan_data_bc_bnp(spec),
    mixture_model    = build_stan_data_mixture_model(spec),
    selection_copas  = build_stan_data_selection_copas(spec),
    selection_weight = build_stan_data_selection_weight(spec),
    pet_peese        = build_stan_data_pet_peese(spec),
    cli::cli_abort("Unknown {.arg model_type}: {.val {spec$model_type}}")
  )

  if (!is.null(spec$custom_data)) {
    for (nm in names(spec$custom_data)) sd_list[[nm]] <- spec$custom_data[[nm]]
  }

  sd_list
}


#' @noRd
build_stan_data_standard <- function(spec) {
  is_re <- spec$model_type == "random_effect"

  if (spec$stage == "two_stage") {
    sd_list <- list(S = spec$S, y = spec$es$yi, se = spec$es$sei)
    if (spec$small_sample %in% c("t_approx", "hjsk")) {
      n_total <- as.integer(spec$n_c + spec$n_i)
      if (any(n_total <= 2)) {
        cli::cli_abort("Student-t adjustment requires total sample size > 2.")
      }
      sd_list$df <- as.numeric(n_total - 2)
    }
    if (spec$re_dist == "mixture" && is_re) {
      sd_list$K <- spec$n_components
      sd_list$prior_dirichlet_alpha <-
        rep(spec$priors$mixture$w$alpha, spec$n_components)
    }
    return(sd_list)
  }

  arm_data <- tibble::tibble(
    study_id = rep(seq_len(spec$S), times = 2),
    treat    = rep(c(0L, 1L), each = spec$S),
    outcome  = c(spec$outcome_ctrl, spec$outcome_int),
    n        = c(spec$n_c, spec$n_i)
  )
  if (spec$likelihood == "gaussian") {
    arm_data <- dplyr::mutate(arm_data,
      sd = c(spec$sd_c, spec$sd_i), se = .data$sd / sqrt(.data$n)
    )
  }

  sd_list <- list(
    N = nrow(arm_data), S = spec$S,
    treat = arm_data$treat, study = arm_data$study_id
  )
  if (spec$likelihood == "binomial") {
    sd_list$events <- as.integer(arm_data$outcome)
    sd_list$n      <- as.integer(arm_data$n)
  } else if (spec$likelihood == "gaussian") {
    sd_list$y  <- arm_data$outcome
    sd_list$se <- arm_data$se
  } else if (spec$likelihood == "poisson") {
    sd_list$events   <- as.integer(arm_data$outcome)
    sd_list$exposure <- as.numeric(arm_data$n)
  }
  if (spec$re_dist == "mixture" && is_re) {
    sd_list$K <- spec$n_components
    sd_list$prior_dirichlet_alpha <-
      rep(spec$priors$mixture$w$alpha, spec$n_components)
  }
  if (isTRUE(spec$multi_arm_config$enabled)) {
    sd_list$n_ma_studies <- spec$multi_arm_config$n_ma_studies
    sd_list$comp_to_ma   <- spec$multi_arm_config$ma_study_id
  }

  attr(sd_list, "arm_data") <- arm_data
  sd_list
}


#' @noRd
build_stan_data_twostage_effect_size <- function(spec) {
  list(N = spec$S, y = spec$es$yi, se = spec$es$sei)
}


#' @noRd
build_stan_data_selection_copas <- function(spec) {
  list(
    N     = spec$S,
    y     = spec$es$yi,
    se    = spec$es$sei,
    s_max = max(spec$es$sei)
  )
}


#' @noRd
build_stan_data_mixture_model <- function(spec) {
  list(
    N = spec$S,
    M = as.integer(spec$n_components),
    y = spec$es$yi,
    se = spec$es$sei
  )
}


#' @noRd
build_stan_data_bc_bnp <- function(spec) {
  cfg <- spec$bc_bnp_config
  mb  <- spec$priors$mu_beta
  list(
    N         = spec$S,
    y         = spec$es$yi,
    se        = spec$es$sei,
    K         = as.integer(cfg$K_max),
    alpha_max = cfg$alpha_max,
    a0        = cfg$a0,
    a1        = cfg$a1,
    B_lower   = if (mb$family == "uniform") mb$lower else -15,
    B_upper   = if (mb$family == "uniform") mb$upper else  15
  )
}


#' @noRd
build_stan_data_bias_corrected <- function(spec) {
  list(
    N              = spec$S,
    y              = spec$es$yi,
    se_y           = spec$es$sei,
    use_known_bias = as.integer(spec$use_known_bias),
    known_bias     = spec$known_bias_vec
  )
}


#' @noRd
build_stan_data_selection_weight <- function(spec) {
  list(
    N         = spec$S,
    y         = spec$es$yi,
    se        = spec$es$sei,
    K         = length(spec$p_cutoffs) + 1L,
    p_cutoffs = sort(spec$p_cutoffs)
  )
}


#' @noRd
build_stan_data_pet_peese <- function(spec) {
  list(
    N       = spec$S,
    y       = spec$es$yi,
    se      = spec$es$sei,
    n_total = as.array(spec$n_c + spec$n_i)
  )
}


# -----------------------------------------------------------------------------
# Stage 4: bayesma_fit -- compile + sample
# -----------------------------------------------------------------------------

#' @noRd
bayesma_fit <- function(spec,
                        code          = bayesma_stan_code(spec),
                        stan_data     = bayesma_stan_data(spec),
                        chains        = 4,
                        iter_warmup   = 1000,
                        iter_sampling = 1000,
                        adapt_delta   = 0.95,
                        seed          = 1234,
                        ...) {
  if (!inherits(spec, "bayesma_spec")) {
    cli::cli_abort("{.arg spec} must be a {.cls bayesma_spec} object.")
  }

  stan_program <- if (inherits(code, "bayesma_stan_code")) code$full
                  else as.character(code)

  sample_model <- function(prog, data_list) {
    mod <- get_cmdstan_model_cached(prog)
    mod$sample(
      data = data_list, chains = chains,
      iter_warmup = iter_warmup, iter_sampling = iter_sampling,
      adapt_delta = adapt_delta, seed = seed,
      refresh = 0, show_messages = FALSE, show_exceptions = FALSE, ...
    )
  }

  if (!is.null(spec$custom_model) || spec$model_type != "pet_peese") {
    fit <- sample_model(stan_program, stan_data)
    out <- list(fit = fit, stan_code = code, stan_data = stan_data)
    class(out) <- c("bayesma_fit", "list")
    return(out)
  }

  pet_code  <- generate_stan_code_pet_peese(spec$priors, predictor = "inv_sqrt_n")
  fit_pet   <- sample_model(pet_code, stan_data)

  pet_draws <- as.vector(
    posterior::subset_draws(fit_pet$draws("mu"), variable = "mu")
  )
  direction   <- if (mean(spec$es$yi) < 0) "negative" else "positive"
  threshold   <- 0.10
  prob_effect <- if (direction == "negative") mean(pet_draws < 0)
                 else                         mean(pet_draws > 0)
  use_peese   <- prob_effect > (1 - threshold)

  if (use_peese) {
    peese_code <- generate_stan_code_pet_peese(spec$priors, predictor = "inv_n")
    fit_final  <- sample_model(peese_code, stan_data)
    stan_program <- peese_code
    recommended <- "peese"
  } else {
    fit_final   <- fit_pet
    recommended <- "pet"
  }

  final_blocks <- parse_stan_blocks(stan_program)
  final_code   <- c(final_blocks, list(full = stan_program))
  class(final_code) <- c("bayesma_stan_code", "list")

  out <- list(
    fit       = fit_final,
    stan_code = final_code,
    stan_data = stan_data,
    pet_peese = list(
      recommended = recommended,
      prob_effect = prob_effect,
      direction   = direction,
      threshold   = threshold,
      fit_pet     = fit_pet
    )
  )
  class(out) <- c("bayesma_fit", "list")
  out
}


# -----------------------------------------------------------------------------
# Stage 5: bayesma_extract -- fit + spec -> tidy effect components
# -----------------------------------------------------------------------------

#' @noRd
bayesma_extract <- function(fit, spec) {
  if (!inherits(fit, "bayesma_fit")) {
    cli::cli_abort("{.arg fit} must be a {.cls bayesma_fit} object.")
  }
  if (!inherits(spec, "bayesma_spec")) {
    cli::cli_abort("{.arg spec} must be a {.cls bayesma_spec} object.")
  }

  eff <- if (!is.null(spec$custom_model)) {
    extract_effects_custom(fit$fit, spec)
  } else {
    switch(spec$model_type,
      random_effect    = extract_effects_standard(fit$fit, spec),
      common_effect    = extract_effects_standard(fit$fit, spec),
      bias_corrected   = extract_effects_bias_corrected(fit$fit, spec),
      bc_bnp           = extract_effects_bc_bnp(fit$fit, spec),
      mixture_model    = extract_effects_mixture_model(fit$fit, spec),
      selection_copas  = extract_effects_selection_copas(fit$fit, spec),
      selection_weight = extract_effects_selection_weight(fit$fit, spec),
      pet_peese        = extract_effects_pet_peese(fit$fit, spec)
    )
  }
  class(eff) <- c("bayesma_effects", "list")
  eff
}


#' @noRd
extract_effects_standard <- function(fit, spec) {
  is_re <- spec$model_type == "random_effect"
  S     <- spec$S

  key_vars <- "mu"
  if (is_re) {
    key_vars <- switch(spec$re_dist,
      normal      = c(key_vars, "tau"),
      t           = c(key_vars, "tau", "nu"),
      skew_normal = c(key_vars, "tau", "alpha_skew"),
      mixture     = {
        if (spec$stage == "two_stage") {
          c(key_vars,
            paste0("mu_k[", seq_len(spec$n_components), "]"),
            paste0("tau_k[", seq_len(spec$n_components), "]"),
            paste0("w[",    seq_len(spec$n_components), "]"))
        } else {
          c(key_vars, "tau",
            paste0("delta_k[", seq_len(spec$n_components), "]"),
            paste0("tau_k[",   seq_len(spec$n_components), "]"),
            paste0("w[",       seq_len(spec$n_components), "]"))
        }
      }
    )
    if (isTRUE(spec$multi_arm_config$enabled) && spec$re_dist != "mixture") {
      key_vars <- c(key_vars, "rho")
    }
  }
  if (isTRUE(spec$robust_config$enabled)) key_vars <- c(key_vars, "pi_main")

  summary_tbl <- stan_summary(fit, variables = key_vars)

  draw_vars <- key_vars
  if (is_re && spec$re_dist %in% c("normal", "t", "skew_normal")) {
    study_pars <- if (spec$stage == "two_stage") {
      paste0("theta[",   seq_len(S), "]")
    } else {
      paste0("epsilon[", seq_len(S), "]")
    }
    draw_vars <- c(draw_vars, study_pars)
  } else if (is_re && spec$re_dist == "mixture" && spec$stage == "two_stage") {
    draw_vars <- c(draw_vars, paste0("cluster[", seq_len(S), "]"))
  }
  if (isTRUE(spec$robust_config$enabled) && is_re &&
      spec$re_dist %in% c("normal", "t", "skew_normal")) {
    draw_vars <- c(draw_vars, paste0("prob_outlier[", seq_len(S), "]"))
  }
  if (is_re) {
    if (tryCatch({ fit$draws("mu_new"); TRUE }, error = function(e) FALSE)) {
      draw_vars <- c(draw_vars, "mu_new")
    }
  }
  draws <- posterior::as_draws_df(fit$draws(variables = draw_vars))

  pooled_draws <- as.vector(
    posterior::subset_draws(fit$draws("mu"), variable = "mu")
  )
  pooled_row <- tibble::tibble(
    study    = "Pooled",
    estimate = stats::median(pooled_draws),
    lower    = stats::quantile(pooled_draws, 0.025),
    upper    = stats::quantile(pooled_draws, 0.975),
    type     = "pooled"
  )

  if (is_re && spec$re_dist %in% c("normal", "t", "skew_normal")) {
    var_prefix <- if (spec$stage == "two_stage") "theta" else "epsilon"
    study_rows <- purrr::map(seq_len(S), function(i) {
      vn <- paste0(var_prefix, "[", i, "]")
      d  <- as.vector(posterior::subset_draws(fit$draws(vn), variable = vn))
      eff <- if (var_prefix == "theta") d else pooled_draws + d
      tibble::tibble(
        study    = spec$study_labels[i],
        estimate = stats::median(eff),
        lower    = stats::quantile(eff, 0.025),
        upper    = stats::quantile(eff, 0.975),
        type     = "study"
      )
    }) |> purrr::list_rbind()
  } else {
    study_rows <- tibble::tibble(
      study    = spec$study_labels,
      estimate = spec$es$yi,
      lower    = spec$es$yi - 1.96 * spec$es$sei,
      upper    = spec$es$yi + 1.96 * spec$es$sei,
      type     = "study"
    )
  }

  if (isTRUE(spec$robust_config$enabled) && is_re &&
      spec$re_dist %in% c("normal", "t", "skew_normal")) {
    study_rows$prob_outlier <- extract_per_study_medians(fit, "prob_outlier", S)
  }

  effect_label <- switch(spec$likelihood,
    binomial = "log_or", gaussian = "mean_diff", poisson = "log_rr"
  )
  forest_df <- dplyr::bind_rows(study_rows, pooled_row) |>
    dplyr::mutate(study        = forcats::fct_inorder(.data$study),
                  effect_scale = effect_label)

  list(
    summary      = summary_tbl,
    draws        = draws,
    forest_df    = forest_df,
    effect_label = effect_label
  )
}


#' @noRd
extract_effects_bias_corrected <- function(fit, spec) {
  S        <- spec$S
  key_vars <- c("mu", "B", "tau", "inv_var", "p_bias")
  draw_vars <- c(key_vars, "mu_biased",
                 paste0("prob_biased[", seq_len(S), "]"))
  extract_effects_simple(
    fit = fit, spec = spec,
    key_vars = key_vars, draw_vars = draw_vars,
    pooled_label = "Pooled (unbiased)",
    extra_study_cols = list(
      prob_biased = extract_per_study_medians(fit, "prob_biased", S)
    )
  )
}


#' @noRd
extract_effects_mixture_model <- function(fit, spec) {
  S <- spec$S
  M <- as.integer(spec$n_components)
  key_vars  <- c("mu",
                 paste0("mu_k[", seq_len(M), "]"),
                 paste0("tau_k[", seq_len(M), "]"),
                 paste0("w[", seq_len(M), "]"))
  draw_vars <- key_vars
  extract_effects_simple(
    fit = fit, spec = spec,
    key_vars = key_vars, draw_vars = draw_vars,
    pooled_label = "Pooled (mixture-weighted mean)"
  )
}


#' @noRd
extract_effects_bc_bnp <- function(fit, spec) {
  S         <- spec$S
  key_vars  <- c("mu_theta", "tau_theta", "mu_beta", "tau_beta",
                 "pi_bias", "alpha")
  draw_vars <- c(key_vars,
                 paste0("prob_biased[", seq_len(S), "]"))
  extract_effects_simple(
    fit = fit, spec = spec,
    key_vars = key_vars, draw_vars = draw_vars,
    pooled_label = "Pooled (bias-corrected, BNP)",
    extra_study_cols = list(
      prob_biased = extract_per_study_medians(fit, "prob_biased", S)
    )
  )
}


#' @noRd
extract_effects_selection_copas <- function(fit, spec) {
  S         <- spec$S
  key_vars  <- c("mu", "tau", "gamma0", "gamma1", "rho")
  draw_vars <- c(key_vars, "selection_bias",
                 paste0("prob_published[", seq_len(S), "]"))
  extract_effects_simple(
    fit = fit, spec = spec,
    key_vars = key_vars, draw_vars = draw_vars,
    pooled_label = "Pooled (selection-adjusted)",
    extra_study_cols = list(
      prob_published = extract_per_study_medians(fit, "prob_published", S)
    )
  )
}


#' @noRd
extract_effects_selection_weight <- function(fit, spec) {
  K <- length(spec$p_cutoffs) + 1L
  key_vars <- c("mu", "tau", paste0("omega[", seq_len(K), "]"))
  extract_effects_simple(
    fit = fit, spec = spec,
    key_vars = key_vars, draw_vars = key_vars,
    pooled_label = "Pooled (weight-adjusted)",
    extra_meta = list(p_cutoffs = spec$p_cutoffs)
  )
}


#' @noRd
extract_effects_pet_peese <- function(fit, spec) {
  key_vars <- c("mu", "beta_bias")
  extract_effects_simple(
    fit = fit, spec = spec,
    key_vars = key_vars, draw_vars = key_vars,
    pooled_label = "Pooled (adjusted)"
  )
}


#' @noRd
extract_effects_simple <- function(fit, spec, key_vars, draw_vars,
                                   pooled_label,
                                   extra_meta = list(),
                                   extra_study_cols = list()) {
  summary_tbl <- stan_summary(fit, variables = key_vars)
  draws       <- posterior::as_draws_df(fit$draws(variables = draw_vars))
  mu_draws    <- as.vector(
    posterior::subset_draws(fit$draws("mu"), variable = "mu")
  )

  pooled_row <- tibble::tibble(
    study    = pooled_label,
    estimate = stats::median(mu_draws),
    lower    = stats::quantile(mu_draws, 0.025),
    upper    = stats::quantile(mu_draws, 0.975),
    type     = "pooled"
  )

  study_rows <- tibble::tibble(
    study    = spec$study_labels,
    estimate = spec$es$yi,
    lower    = spec$es$yi - 1.96 * spec$es$sei,
    upper    = spec$es$yi + 1.96 * spec$es$sei,
    type     = "study"
  )
  for (nm in names(extra_study_cols)) study_rows[[nm]] <- extra_study_cols[[nm]]

  forest_df <- dplyr::bind_rows(study_rows, pooled_row) |>
    dplyr::mutate(study = forcats::fct_inorder(.data$study))

  list(
    summary    = summary_tbl,
    draws      = draws,
    forest_df  = forest_df,
    extra_meta = extra_meta
  )
}


#' @noRd
extract_effects_custom <- function(fit, spec) {
  draws_raw <- posterior::as_draws_df(fit$draws())
  draw_vars <- setdiff(names(draws_raw), c(".chain", ".iteration", ".draw"))
  summary_tbl <- tryCatch(
    stan_summary(fit),
    error = function(e) NULL
  )

  forest_df <- NULL
  if ("mu" %in% draw_vars) {
    mu_draws <- as.vector(
      posterior::subset_draws(fit$draws("mu"), variable = "mu")
    )
    pooled_row <- tibble::tibble(
      study    = "Pooled",
      estimate = stats::median(mu_draws),
      lower    = stats::quantile(mu_draws, 0.025),
      upper    = stats::quantile(mu_draws, 0.975),
      type     = "pooled"
    )
    study_rows <- tibble::tibble(
      study    = spec$study_labels,
      estimate = spec$es$yi,
      lower    = spec$es$yi - 1.96 * spec$es$sei,
      upper    = spec$es$yi + 1.96 * spec$es$sei,
      type     = "study"
    )
    forest_df <- dplyr::bind_rows(study_rows, pooled_row) |>
      dplyr::mutate(study = forcats::fct_inorder(.data$study))
  }

  list(
    summary      = summary_tbl,
    draws        = draws_raw,
    forest_df    = forest_df,
    effect_label = "custom"
  )
}


# -----------------------------------------------------------------------------
# Stage 6: bayesma_output -- assemble the final bayesma object
# -----------------------------------------------------------------------------

#' @noRd
bayesma_output <- function(spec, fit, effects) {
  if (!inherits(spec,    "bayesma_spec"))    cli::cli_abort("{.arg spec} must be {.cls bayesma_spec}.")
  if (!inherits(fit,     "bayesma_fit"))     cli::cli_abort("{.arg fit} must be {.cls bayesma_fit}.")
  if (!inherits(effects, "bayesma_effects")) cli::cli_abort("{.arg effects} must be {.cls bayesma_effects}.")

  code_full <- if (inherits(fit$stan_code, "bayesma_stan_code")) fit$stan_code$full
               else as.character(fit$stan_code)

  arm_data <- attr(fit$stan_data, "arm_data")
  if (is.null(arm_data)) {
    arm_data <- tibble::tibble(
      study = spec$study_labels,
      yi    = spec$es$yi,
      sei   = spec$es$sei
    )
  }

  meta <- list(
    likelihood   = spec$likelihood,
    model_type   = spec$model_type,
    re_dist      = spec$re_dist,
    stage        = spec$stage,
    small_sample = spec$small_sample,
    study_labels = spec$study_labels,
    priors       = spec$priors,
    effect_label = effects$effect_label %||% "custom",
    robust       = isTRUE(spec$robust_config$enabled),
    multi_arm    = isTRUE(spec$multi_arm_config$enabled),
    es           = spec$es,
    call_args    = spec$call_args
  )
  if (isTRUE(spec$robust_config$enabled))    meta$robust_config    <- spec$robust_config
  if (isTRUE(spec$multi_arm_config$enabled)) meta$multi_arm_config <- spec$multi_arm_config
  if (spec$model_type == "bias_corrected")   meta$use_known_bias   <- spec$use_known_bias
  if (!is.null(spec$outcome_ctrl))           meta$outcome_ctrl     <- spec$outcome_ctrl
  if (!is.null(spec$n_c))                    meta$n_c              <- spec$n_c

  if (!is.null(effects$extra_meta)) {
    for (nm in names(effects$extra_meta)) meta[[nm]] <- effects$extra_meta[[nm]]
  }
  if (!is.null(fit$pet_peese)) {
    meta$recommended <- fit$pet_peese$recommended
    meta$prob_effect <- fit$pet_peese$prob_effect
    meta$direction   <- fit$pet_peese$direction
    meta$threshold   <- fit$pet_peese$threshold
    meta$fit_pet     <- fit$pet_peese$fit_pet
  }
  if (!is.null(spec$custom_model)) meta$custom_model <- TRUE

  out <- list(
    fit       = fit$fit,
    summary   = effects$summary,
    forest_df = effects$forest_df,
    draws     = effects$draws,
    stan_code = code_full,
    stan_data = fit$stan_data,
    arm_data  = arm_data,
    meta      = meta
  )
  class(out) <- "bayesma"
  out
}


# -----------------------------------------------------------------------------
# Stan code generators for bias / selection / PET-PEESE
# -----------------------------------------------------------------------------

#' @noRd
generate_stan_code_selection_copas <- function(priors, re_dist = "t") {
  if (!re_dist %in% c("normal", "t")) {
    cli::cli_abort(
      "{.arg re_dist} {.val {re_dist}} not supported for {.val selection_copas}. \\
       Use {.val normal} or {.val t}."
    )
  }

  p  <- priors
  sp <- p$selection %||% list()
  if (is.null(sp$gamma0)) sp$gamma0 <- uniform(-2, 2)
  if (is.null(sp$rho))    sp$rho    <- uniform(-1, 1)

  mu_tgt   <- emit_prior_target(p$mu, "mu")
  tau_tgt  <- emit_prior_target(p$tau, "tau")
  g0_tgt   <- emit_prior_target(sp$gamma0, "gamma0")
  g1_tgt   <- emit_prior_target(sp$gamma1, "gamma1")
  rho_tgt  <- emit_prior_target(sp$rho, "rho")
  nu_tgt   <- if (re_dist == "t") emit_prior_target(p$nu, "nu") else ""
  tau_bnds <- emit_prior_bounds(p$tau, default_lower = 0)
  rho_bnds <- emit_prior_bounds(sp$rho, default_lower = -1, default_upper = 1)
  g0_bnds  <- emit_prior_bounds(sp$gamma0)

  nu_par <- if (re_dist == "t") "  real<lower=2> nu;\n" else ""
  u_prior <- if (re_dist == "t") {
    "  u ~ student_t(nu, 0, 1);"
  } else {
    "  u ~ std_normal();"
  }
  mu_new_rng <- if (re_dist == "t") {
    "real mu_new = mu + tau * student_t_rng(nu, 0, 1);"
  } else {
    "real mu_new = normal_rng(mu, tau);"
  }

  paste0(
    "data {
  int<lower=1> N;
  vector[N] y;
  vector<lower=0>[N] se;
  real<lower=0> s_max;
}
parameters {
  real mu;
  real", tau_bnds, " tau;
  real", g0_bnds, " gamma0;
  real<lower=0, upper=s_max> gamma1;
  real", rho_bnds, " rho;
  vector[N] u;
", nu_par, "}
transformed parameters {
  vector[N] a;
  vector[N] mu_i;
  for (i in 1:N) {
    a[i]    = gamma0 + gamma1 / se[i];
    mu_i[i] = mu + tau * u[i];
  }
}
model {
  ", mu_tgt, "
  ", tau_tgt, "
  ", g0_tgt, "
  ", g1_tgt, "
  ", rho_tgt, "
  ", nu_tgt, "
", u_prior, "
  for (i in 1:N) {
    real arg = (a[i] + rho * (y[i] - mu_i[i]) / se[i]) / sqrt(1 - square(rho));
    target += normal_lpdf(y[i] | mu_i[i], se[i]);
    target += normal_lcdf(arg | 0, 1);
    target += -normal_lcdf(a[i] | 0, 1);
  }
}
generated quantities {
  real pooled = mu;
  vector[N] prob_published;
  for (i in 1:N)
    prob_published[i] = Phi(a[i]);
  ", mu_new_rng, "
  real selection_bias;
  {
    real avg_se = mean(se);
    real avg_a  = mean(a);
    selection_bias = rho * avg_se * exp(std_normal_lpdf(avg_a)) / Phi(avg_a);
  }
}")
}


#' @noRd
generate_stan_code_mixture_model <- function(priors, M) {
  if (M < 2L) {
    cli::cli_abort(
      "{.arg n_components} must be >= 2 for {.val mixture_model}.",
      call = rlang::caller_env()
    )
  }
  mp <- priors$mixture
  mu_k_tgt  <- emit_prior_target(mp$mu_k, "mu_k[m]")
  tau_k_tgt <- emit_prior_target(mp$tau_k, "tau_k[m]")
  tau_k_bnds <- emit_prior_bounds(mp$tau_k, default_lower = 0)
  w_alpha <- if (!is.null(mp$w) && mp$w$family == "dirichlet" &&
                 length(mp$w$alpha) == M) {
    paste0("[", paste(mp$w$alpha, collapse = ", "), "]'")
  } else if (!is.null(mp$w) && mp$w$family == "dirichlet" &&
             length(mp$w$alpha) == 1) {
    paste0("rep_vector(", mp$w$alpha, ", M)")
  } else {
    "rep_vector(1.0, M)"
  }

  paste0(
    "data {
  int<lower=2> M;
  int<lower=1> N;
  vector[N] y;
  vector<lower=0>[N] se;
}
parameters {
  ordered[M] mu_k;
  vector", tau_k_bnds, "[M] tau_k;
  simplex[M] w;
}
model {
  for (m in 1:M) {
    ", mu_k_tgt, "
    ", tau_k_tgt, "
  }
  w ~ dirichlet(", w_alpha, ");
  for (i in 1:N) {
    vector[M] log_comp;
    for (m in 1:M) {
      log_comp[m] = log(w[m])
                   + normal_lpdf(y[i] | mu_k[m],
                                 sqrt(square(se[i]) + square(tau_k[m])));
    }
    target += log_sum_exp(log_comp);
  }
}
generated quantities {
  real mu     = sum(w .* mu_k);
  real pooled = mu;
  matrix[N, M] prob_cluster;
  for (i in 1:N) {
    vector[M] log_comp;
    for (m in 1:M) {
      log_comp[m] = log(w[m])
                   + normal_lpdf(y[i] | mu_k[m],
                                 sqrt(square(se[i]) + square(tau_k[m])));
    }
    real lse = log_sum_exp(log_comp);
    for (m in 1:M)
      prob_cluster[i, m] = exp(log_comp[m] - lse);
  }
}")
}


#' @noRd
generate_stan_code_bc_bnp <- function(priors, K) {
  p <- priors

  mu_tgt        <- emit_prior_target(p$mu, "mu_theta")
  tau_tgt       <- emit_prior_target(p$tau, "tau_theta")
  tau_beta_tgt  <- emit_prior_target(p$tau_beta, "tau_beta")
  pi_bias_tgt   <- emit_prior_target(p$p_bias, "pi_bias")
  tau_bnds      <- emit_prior_bounds(p$tau, default_lower = 0)
  tau_beta_bnds <- emit_prior_bounds(p$tau_beta, default_lower = 0)

  paste0(
    "data {
  int<lower=1> N;
  vector[N] y;
  vector<lower=0>[N] se;
  int<lower=1> K;
  real<lower=0.5> alpha_max;
  real<lower=0> a0;
  real<lower=0> a1;
  real B_lower;
  real B_upper;
}
parameters {
  real mu_theta;
  real", tau_bnds, " tau_theta;
  real<lower=B_lower, upper=B_upper> mu_beta;
  real", tau_beta_bnds, " tau_beta;
  real<lower=0, upper=1> pi_bias;
  real<lower=0.5, upper=alpha_max> alpha;
  vector<lower=0, upper=1>[K - 1] v;
  ordered[K] beta_star;
  vector[N] theta_z;
}
transformed parameters {
  simplex[K] w;
  vector[N] theta = mu_theta + tau_theta * theta_z;
  {
    real remaining = 1;
    for (k in 1:(K - 1)) {
      w[k] = v[k] * remaining;
      remaining = remaining * (1 - v[k]);
    }
    w[K] = remaining;
  }
}
model {
  ", mu_tgt, "
  ", tau_tgt, "
  ", tau_beta_tgt, "
  ", pi_bias_tgt, "
  for (k in 1:(K - 1))
    target += beta_lpdf(v[k] | 1, alpha);
  beta_star ~ normal(mu_beta, tau_beta);
  theta_z   ~ std_normal();
  for (i in 1:N) {
    vector[K + 1] log_comp;
    log_comp[1] = log1m(pi_bias) + normal_lpdf(y[i] | theta[i], se[i]);
    for (k in 1:K) {
      log_comp[k + 1] = log(pi_bias) + log(w[k])
                       + normal_lpdf(y[i] | theta[i] + beta_star[k], se[i]);
    }
    target += log_sum_exp(log_comp);
  }
}
generated quantities {
  real mu     = mu_theta;
  real pooled = mu_theta;
  vector[N] prob_biased;
  for (i in 1:N) {
    vector[K + 1] log_comp;
    log_comp[1] = log1m(pi_bias) + normal_lpdf(y[i] | theta[i], se[i]);
    for (k in 1:K) {
      log_comp[k + 1] = log(pi_bias) + log(w[k])
                       + normal_lpdf(y[i] | theta[i] + beta_star[k], se[i]);
    }
    real lse = log_sum_exp(log_comp);
    prob_biased[i] = 1 - exp(log_comp[1] - lse);
  }
  real mu_new = normal_rng(mu_theta, tau_theta);
}")
}


#' @noRd
generate_stan_code_selection_weight <- function(priors, K) {
  p <- priors
  mu_tgt   <- emit_prior_target(p$mu, "mu")
  tau_tgt  <- emit_prior_target(p$tau, "tau")
  tau_bnds <- emit_prior_bounds(p$tau, default_lower = 0)

  paste0(
    "data {
  int<lower=1> N;
  vector[N] y;
  vector<lower=0>[N] se;
  int<lower=2> K;
  vector[K-1] p_cutoffs;
}
transformed data {
  vector[K+1] z_bounds;
  z_bounds[1] = 8.0;
  z_bounds[K+1] = -8.0;
  for (k in 1:(K-1))
    z_bounds[k+1] = inv_Phi(1.0 - p_cutoffs[k]);
}
parameters {
  real mu;
  real", tau_bnds, " tau;
  vector<lower=0.01, upper=0.99>[K-1] omega_raw;
}
transformed parameters {
  vector[K] omega;
  vector[N] sigma;
  omega[1] = 1.0;
  for (k in 1:(K-1))
    omega[k+1] = omega_raw[k];
  for (i in 1:N)
    sigma[i] = sqrt(square(tau) + square(se[i]));
}
model {
  ", mu_tgt, "
  ", tau_tgt, "
  for (k in 1:(K-1))
    target += beta_lpdf(omega_raw[k] | 1, 1);
  for (i in 1:N) {
    real z_i = y[i] / se[i];
    real p_val = 1.0 - Phi(z_i);
    real log_w;
    if (p_val < p_cutoffs[1]) {
      log_w = 0.0;
    } else if (K == 2) {
      log_w = log(omega[2]);
    } else {
      log_w = log(omega[K]);
      for (k in 1:(K-2)) {
        if (p_val >= p_cutoffs[k] && p_val < p_cutoffs[k+1])
          log_w = log(omega[k+1]);
      }
    }
    target += normal_lpdf(y[i] | mu, sigma[i]) + log_w;
    {
      real norm_c = 0;
      for (k in 1:K) {
        real prob_k = Phi((z_bounds[k] * se[i] - mu) / sigma[i])
                    - Phi((z_bounds[k+1] * se[i] - mu) / sigma[i]);
        norm_c += omega[k] * fmax(prob_k, 1e-15);
      }
      target += -log(fmax(norm_c, 1e-15));
    }
  }
}
generated quantities {
  real pooled = mu;
  vector[K] weights = omega;
  real mu_new = normal_rng(mu, tau);
}")
}


#' @noRd
generate_stan_code_pet_peese <- function(priors, predictor = c("inv_sqrt_n", "inv_n")) {
  predictor <- match.arg(predictor)
  p <- priors
  mu_tgt <- emit_prior_target(p$mu, "mu")
  sp <- p$selection %||% list()
  if (is.null(sp$beta_bias))
    sp$beta_bias <- if (predictor == "inv_sqrt_n") normal(0, 1) else normal(0, 2)
  beta_tgt <- emit_prior_target(sp$beta_bias, "beta_bias")

  paste0(
    "data {
  int<lower=1> N;
  vector[N] y;
  vector<lower=0>[N] se;
  vector<lower=0>[N] n_total;
}
transformed data {
  vector[N] inv_sqrt_n;
  vector[N] inv_n;
  for (i in 1:N) {
    inv_sqrt_n[i] = 1.0 / sqrt(n_total[i]);
    inv_n[i] = 1.0 / n_total[i];
  }
}
parameters {
  real mu;
  real beta_bias;
}
model {
  ", mu_tgt, "
  ", beta_tgt, "
  for (i in 1:N)
    target += normal_lpdf(y[i] | mu + beta_bias * ", predictor, "[i], se[i]);
}
generated quantities {
  real pooled = mu;
  real bias_slope = beta_bias;
}")
}
