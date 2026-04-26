#' Create INSPECT-SR Trustworthiness Plot
#'
#' @param data An \code{inspect_sr} object or a data frame with columns
#'   Study, D1, D2, D3, D4, Overall.
#' @param sort_studies_by \code{"study"} (default, alphabetical) or
#'   \code{"overall"} (by severity).
#' @param add_legend Logical. Add a legend panel (default FALSE).
#' @param incl_checks Logical. If \code{TRUE}, expand the table to show every
#'   INSPECT-SR item (1.1, 1.2, …, 4.11) grouped under domain spanners, with
#'   a per-domain "Overall" and a final study-level "Overall". Requires
#'   \code{data} to be an \code{inspect_sr} object. Default \code{FALSE}.
#' @param font Character. Font family (NULL = default).
#' @param title,subtitle Character. Optional.
#' @param title_align \code{"left"} (default), \code{"center"}/\code{"centre"},
#'   or \code{"right"}.
#'
#' @return A patchwork object.
#' @export
inspect_plot <- function(data,
                         sort_studies_by = "study",
                         add_legend = FALSE,
                         incl_checks = FALSE,
                         font = NULL,
                         title = NULL,
                         title_align = "left",
                         subtitle = NULL) {

  if (is.null(data)) cli::cli_abort("{.arg data} must be provided.")

  if (isTRUE(incl_checks) && !inherits(data, "inspect_sr")) {
    cli::cli_abort(
      "{.arg incl_checks = TRUE} requires an {.cls inspect_sr} object (per-item columns are read from {.code attr(data, 'data')})."
    )
  }

  df <- if (inherits(data, "inspect_sr")) as.data.frame(data) else data

  required_cols <- c("Study", "D1", "D2", "D3", "D4", "Overall")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    cli::cli_abort("Missing columns: {.val {missing_cols}}.")
  }

  # If per-item view requested, merge the per-item columns from the input data
  # plus fill in the four automated items from the details attribute.
  if (isTRUE(incl_checks)) {
    input_df <- as.data.frame(attr(data, "data"))
    details <- attr(data, "details")
    det_by_study <- stats::setNames(details, vapply(details, `[[`, character(1), "study"))

    item_cols <- c(paste0("d1_", 1:3), paste0("d2_", 1:5),
                   paste0("d3_", 1:2), paste0("d4_", 1:11))
    input_df$Study <- input_df$study
    item_df <- input_df[, c("Study", item_cols), drop = FALSE]

    # Fill automated items from details$<check>$judgement.
    auto_map <- list(d4_3 = "carlisle", d4_6 = "n_consistency",
                     d4_8 = "grim",      d4_9 = "pvalue")
    for (col in names(auto_map)) {
      check_name <- auto_map[[col]]
      vals <- vapply(item_df$Study, function(s) {
        j <- det_by_study[[s]][[check_name]]$judgement
        if (is.null(j)) NA_character_ else as.character(j)
      }, character(1))
      item_df[[col]] <- vals
    }

    df <- merge(item_df, df, by = "Study", all.x = TRUE)
  }

  # Sort

  if (sort_studies_by == "overall") {
    sev <- c("No concerns" = 1, "Some concerns" = 2, "Serious concerns" = 3)
    df$sort_key <- sev[df$Overall]
    df$sort_key[is.na(df$sort_key)] <- 0
    df <- df[order(df$sort_key, df$Study), ]
    df$sort_key <- NULL
  } else {
    df <- df[order(df$Study), ]
  }

  if (isTRUE(incl_checks)) {
    item_cols <- c(paste0("d1_", 1:3), paste0("d2_", 1:5),
                   paste0("d3_", 1:2), paste0("d4_", 1:11))
    display_cols <- c("Study", item_cols, "D1", "D2", "D3", "D4", "Overall")
    df <- df[, display_cols, drop = FALSE]
    inspect_table <- inspect_sr_item_table_fn(df, font = font)
  } else {
    df <- df[, required_cols, drop = FALSE]
    inspect_table <- inspect_sr_table_fn(df, font = font)
  }

  if (isTRUE(add_legend)) {
    legend <- inspect_sr_legend_fn(font)
    plot_out <- patchwork::wrap_table(inspect_table, space = "fixed") +
      patchwork::wrap_table(legend, space = "fixed")
  } else {
    plot_out <- patchwork::wrap_table(inspect_table, space = "fixed")
  }

  if (!is.null(title) || !is.null(subtitle)) {
    hjust_val <- switch(title_align, "left" = 0, "center" = 0.5,
                        "centre" = 0.5, "right" = 1, 0.5)
    title_theme <- ggplot2::theme()
    if (!is.null(title)) {
      title_theme <- title_theme + ggplot2::theme(
        plot.title = ggplot2::element_text(
          size = 16, face = "bold", hjust = hjust_val,
          margin = ggplot2::margin(b = if (is.null(subtitle)) 10 else 5),
          family = font
        ))
    }
    if (!is.null(subtitle)) {
      title_theme <- title_theme + ggplot2::theme(
        plot.subtitle = ggplot2::element_text(
          size = 14, hjust = hjust_val,
          margin = ggplot2::margin(b = 10),
          color = "gray30", family = font
        ))
    }
    plot_out <- plot_out + patchwork::plot_annotation(
      title = title, subtitle = subtitle, theme = title_theme
    )
  }

  plot_out
}


#' Create INSPECT-SR Summary Bar Plot
#'
#' @param data An \code{inspect_sr} object or data frame.
#' @param font,title As in \code{\link{inspect_plot}}.
#' @return A ggplot2 object.
#' @export
inspect_summary_plot <- function(data, font = NULL, title = NULL) {

  df <- if (inherits(data, "inspect_sr")) as.data.frame(data) else data

  long_df <- data.frame(domain = character(), judgement = character(),
                        stringsAsFactors = FALSE)
  for (col in c("D1", "D2", "D3", "D4", "Overall")) {
    vals <- df[[col]]
    vals[is.na(vals)] <- "Not assessed"
    long_df <- rbind(long_df, data.frame(domain = col, judgement = vals,
                                         stringsAsFactors = FALSE))
  }

  domain_labels <- c(D1 = "D1: Post-publication notices",
                     D2 = "D2: Conduct & governance",
                     D3 = "D3: Text & figures",
                     D4 = "D4: Results", Overall = "Overall")

  long_df$domain <- factor(
    long_df$domain,
    levels = rev(c("D1", "D2", "D3", "D4", "Overall")),
    labels = rev(domain_labels[c("D1", "D2", "D3", "D4", "Overall")])
  )
  long_df$judgement <- factor(
    long_df$judgement,
    levels = c("No concerns", "Some concerns", "Serious concerns", "Not assessed")
  )

  colours <- c("No concerns" = "#77bb41", "Some concerns" = "#f5ec00",
               "Serious concerns" = "#e32400", "Not assessed" = "#cccccc")

  ggplot2::ggplot(long_df, ggplot2::aes(x = domain, fill = judgement)) +
    ggplot2::geom_bar(position = "fill", width = 0.7) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = colours, name = "Judgement", drop = FALSE) +
    ggplot2::scale_y_continuous(labels = scales::percent_format()) +
    ggplot2::labs(x = NULL, y = "Proportion of studies", title = title) +
    ggplot2::theme_minimal(base_family = font %||% "") +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 11, hjust = 1),
      plot.title = ggplot2::element_text(face = "bold", size = 14)
    )
}
