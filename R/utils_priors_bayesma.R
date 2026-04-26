#' Extract Tau Prior from Model
#'
#' @description
#' Extracts the prior specification for the heterogeneity parameter (tau)
#' from a brms or bayesma model.
#'
#' @param model A fitted brmsfit or bayesma object
#'
#' @return Character string of the tau prior, or NA if not found
#'
#' @keywords internal
#' @noRd
extract_tau_from_model <- function(model) {

  if (inherits(model, "brmsfit")) {
    tau_prior <- model$prior |>
      dplyr::filter(class == "sd") |>
      dplyr::pull(prior)

    if (length(tau_prior) == 0) {
      NA_character_
    } else {
      tau_prior[[1]]
    }

  } else if (inherits(model, "bayesma")) {
    # Extract tau prior from stored priors in meta
    tau_prior <- model$meta$priors$tau
    if (is.null(tau_prior)) return(NA_character_)

    # Convert bayesma_prior to a string representation
    format(tau_prior)

  } else {
    NA_character_
  }
}

#' Convert Prior String to Unicode Format
#'
#' @description
#' Converts prior distribution strings to Unicode mathematical notation
#' for pretty printing in tables. Works with both brms prior strings and
#' bayesma prior format strings.
#'
#' @param prior_string Character string representing a prior distribution
#'
#' @return Unicode-formatted prior string
#'
#' @details
#' Supports conversion of:
#' \itemize{
#'   \item normal() / N() to 𝒩()
#'   \item cauchy() / C() to 𝒞()
#'   \item student_t() / t() to 𝓉()
#'   \item exponential() / Exp() to ℰ()
#'   \item HN() to 𝒩⁺()  (half-normal)
#' }
#'
#' @keywords internal
#' @noRd
prior_to_unicode <- function(prior_string) {
  if (is.null(prior_string) || is.na(prior_string)) {
    return(NA_character_)
  }

  # Unicode math-script letters (STIX Two Math)
  script_N <- "\U0001D4A9"  # 𝒩
  script_C <- "\U0001D49E"  # 𝒞
  script_H <- "\U0001D4D7"  # 𝓗
  script_t <- "\U0001D4C9"  # 𝓉
  script_E <- "\U2130"      # ℰ

  prior_string <- trimws(prior_string)

  # ---- bayesma format strings ----
  # HN(mean, sd) -> half-normal
  if (grepl("^HN\\s*\\(", prior_string)) {
    params <- sub("^HN\\s*\\(([^)]+)\\).*$", "\\1", prior_string)
    params <- gsub("\\s+", "", params)
    return(paste0(script_N, "\u207A(", params, ")\u2009"))
  }

  # N(mean, sd) -> normal (bayesma format)
  if (grepl("^N\\s*\\(", prior_string)) {
    params <- sub("^N\\s*\\(([^)]+)\\).*$", "\\1", prior_string)
    params <- gsub("\\s+", "", params)
    return(paste0(script_N, "(", params, ")\u2009"))
  }

  # HC(location, scale) -> half-cauchy (bayesma format)
  if (grepl("^HC\\s*\\(", prior_string)) {
    params <- sub("^HC\\s*\\(([^)]+)\\).*$", "\\1", prior_string)
    params <- gsub("\\s+", "", params)
    return(paste0(script_H,script_C, "(", params, ")\u2009"))
  }

  # t(df, location, scale) -> half-student-t (bayesma format)
  if (grepl("^t\\s*\\(", prior_string)) {
    params <- sub("^t\\s*\\(([^)]+)\\).*$", "\\1", prior_string)
    parts <- strsplit(params, ",")[[1]]
    parts <- trimws(parts)
    if (length(parts) == 3) {
      df  <- gsub("\\s+", "", parts[1])
      loc <- gsub("\\s+", "", parts[2])
      scl <- gsub("\\s+", "", parts[3])
      if (nchar(df) == 1 && grepl("^[0-9]$", df)) {
        df <- c(
          "0" = "\u2080", "1" = "\u2081", "2" = "\u2082", "3" = "\u2083",
          "4" = "\u2084", "5" = "\u2085", "6" = "\u2086", "7" = "\u2087",
          "8" = "\u2088", "9" = "\u2089"
        )[df]
      }
      return(paste0(script_t, df, "(", loc, ", ", scl, ")\u2009"))
    }
  }

  # Exp(rate) -> exponential (bayesma format)
  if (grepl("^Exp\\s*\\(", prior_string)) {
    params <- sub("^Exp\\s*\\(([^)]+)\\).*$", "\\1", prior_string)
    params <- gsub("\\s+", "", params)
    return(paste0(script_E, "(", params, ")\u2009"))
  }

  # ---- brms format strings (original logic) ----
  # Normal
  if (grepl("^normal\\s*\\(", prior_string)) {
    params <- sub("^normal\\s*\\(([^)]+)\\).*$", "\\1", prior_string)
    params <- gsub("\\s+", "", params)
    return(paste0(script_N, "(", params, ")\u2009"))
  }

  # Cauchy (ignore truncation in display)
  if (grepl("^cauchy\\s*\\(", prior_string)) {
    base <- sub("\\s*\\[.*$", "", prior_string)
    params <- sub("^cauchy\\s*\\(([^)]+)\\).*$", "\\1", base)
    params <- gsub("\\s+", "", params)
    return(paste0(script_C, "(", params, ")\u2009"))
  }

  # Student-t
  if (grepl("^student_t\\s*\\(", prior_string)) {
    clean <- sub("\\[.*$", "", prior_string)
    params <- sub("^student_t\\s*\\(([^)]+)\\).*$", "\\1", clean)
    parts <- strsplit(params, ",")[[1]]
    parts <- trimws(parts)

    if (length(parts) != 3) {
      return(prior_string)
    }

    df  <- gsub("\\s+", "", parts[1])
    loc <- gsub("\\s+", "", parts[2])
    scl <- gsub("\\s+", "", parts[3])

    if (nchar(df) == 1 && grepl("^[0-9]$", df)) {
      df <- c(
        "0" = "\u2080", "1" = "\u2081", "2" = "\u2082", "3" = "\u2083",
        "4" = "\u2084", "5" = "\u2085", "6" = "\u2086", "7" = "\u2087",
        "8" = "\u2088", "9" = "\u2089"
      )[df]
    }

    return(paste0(script_t, df, "(", loc, ", ", scl, ")\u2009"))
  }

  # Exponential
  if (grepl("^exponential\\s*\\(", prior_string)) {
    params <- sub("^exponential\\s*\\(([^)]+)\\).*$", "\\1", prior_string)
    params <- gsub("\\s+", "", params)
    return(paste0(script_E, "(", params, ")\u2009"))
  }

  prior_string
}

#' Extract Mu and Tau Priors
#'
#' @description
#' Extracts and formats the prior specifications for both the mean (mu) and
#' heterogeneity (tau) parameters from a list of priors. Supports both
#' brms priors (brmsprior objects) and bayesma priors (named lists of
#' bayesma_prior objects).
#'
#' @param priors Named list of prior specifications
#' @param model A fitted brmsfit or bayesma object (for extracting default tau prior)
#'
#' @return A tibble with columns: prior, prior_label, mu_prior_unicode, tau_prior_unicode
#'
#' @keywords internal
#' @noRd
extract_mu_tau_priors <- function(priors, model) {

  model_tau <- extract_tau_from_model(model)
  is_bayesma <- inherits(model, "bayesma")

  purrr::imap(priors, function(prior_obj, prior_name) {

    if (is_bayesma) {
      # bayesma priors: named list with bayesma_prior objects
      mu_prior_str <- if (!is.null(prior_obj$mu_prior)) {
        format(prior_obj$mu_prior)
      } else {
        NA_character_
      }

      tau_prior_str <- if (!is.null(prior_obj$tau_prior)) {
        format(prior_obj$tau_prior)
      } else {
        NA_character_
      }

    } else {
      # brms priors: brmsprior data frames
      mu_prior <- prior_obj |>
        dplyr::filter(class == "Intercept") |>
        dplyr::pull(prior)

      tau_prior <- prior_obj |>
        dplyr::filter(class == "sd") |>
        dplyr::pull(prior)

      mu_prior_str <- if (length(mu_prior) > 0) mu_prior[[1]] else NA_character_
      tau_prior_str <- if (length(tau_prior) > 0) tau_prior[[1]] else NA_character_
    }

    tibble::tibble(
      prior = prior_name,
      prior_label = dplyr::recode(
        prior_name,
        vague = "Vague",
        weakreg = "Weakly Regularising",
        informative = "Informative",
        .default = stringr::str_to_sentence(prior_name)
      ),
      mu_prior_unicode =
        if (!is.na(mu_prior_str)) prior_to_unicode(mu_prior_str) else NA_character_,
      tau_prior_unicode =
        if (!is.na(tau_prior_str)) {
          prior_to_unicode(tau_prior_str)
        } else if (!is.na(model_tau)) {
          prior_to_unicode(model_tau)
        } else {
          NA_character_
        }
    )
  }) |> purrr::list_rbind()
}

#' Apply Math Font to GT Table
#'
#' @description
#' Applies a mathematical font (like STIX Two Math) to specific columns
#' in a gt table for proper rendering of mathematical symbols.
#'
#' @param gt_tbl A gt table object
#' @param columns Column names to apply the math font to
#' @param math_font Character string. Name of the math font. Default is "STIX Two Math"
#'
#' @return Modified gt table object
#'
#' @keywords internal
#' @noRd
apply_math_font <- function(
    gt_tbl,
    columns,
    math_font = "STIX Two Math"
) {
  gt_tbl |>
    gt::tab_style(
      style = gt::cell_text(font = math_font),
      locations = gt::cells_body(columns = columns))
}

#' Extract Priors from RoBMA Fit
#'
#' @description
#' Extracts and formats prior specifications from a RoBMA_brms fit object
#' for display in tables.
#'
#' @param robma_fit A RoBMA_brms fit object
#'
#' @return List with mu_prior_unicode and tau_prior_unicode
#'
#' @keywords internal
#' @noRd
extract_priors_from_robma_fit <- function(robma_fit) {

  prior_to_label <- function(p) {
    fam <- tolower(p$family %||% "")
    if (fam %in% c("normal", "gaussian")) {
      sprintf("N(%.2f, %.2f)", p$mean %||% 0, p$sd %||% NA_real_)
    } else if (fam %in% c("student_t", "student-t", "t", "student")) {
      sprintf("t(df=%.1f, loc=%.2f, scale=%.2f)", p$df, p$location, p$scale)
    } else {
      fam
    }
  }

  if (!is.null(robma_fit$meta$priors)) {

    eff <- robma_fit$meta$priors$effect %||% list()
    het <- robma_fit$meta$priors$heterogeneity %||% list()

    mu_labels  <- unique(purrr::map_chr(eff, prior_to_label))
    tau_labels <- unique(purrr::map_chr(het, prior_to_label))

    mu_prior_unicode <- paste(mu_labels, collapse = " / ")
    tau_prior_unicode <- paste(tau_labels, collapse = " / ")

    return(list(
      mu_prior_unicode = mu_prior_unicode,
      tau_prior_unicode = tau_prior_unicode
    ))
  }

  rlang::abort(c(
    "Cannot extract priors from this RoBMA object.",
    "i" = "Expected robma_fit$meta$priors$effect and $heterogeneity."
  ))
}
