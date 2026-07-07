# Model performance tab visualisations.


#' Baseline ML model performance plot
#'
#' Grouped box plots of a performance metric by species and molecular scale,
#' restricted to baseline (non-stratified) models. Uses the amRml parquet
#' columns species, feature_type (scale), feature_subtype (data type),
#' drug_or_class, drug_label, and the metric columns (mcc/bal_acc/f1).
#'
#' @param data Performance tibble from loadMLResults() / queryData().
#' @param bug Species code(s) to include.
#' @param model_scale Molecular scale(s) (feature_type) to include.
#' @param data_type Data encoding(s) (feature_subtype) to include.
#' @param metrics Name of the metric column to plot on the y-axis.
#' @param amr_drug_class Selected drug class(es), or "all" for no filter.
#' @param amr_drug Selected drug(s); points for these are overlaid.
#' @return A plotly box-plot figure (an empty placeholder when there is no
#'   matching data).
#' @keywords internal
#' @noRd
makeModelPerformancePlot <- function(
  data, bug, model_scale, data_type, metrics,
  amr_drug_class, amr_drug
) {
  if (is.null(data) || !is.data.frame(data) || !nrow(data)) {
    return(plotly::plot_ly() |>
      plotly::layout(title = list(text = "No data available", x = 0)))
  }

  # Filter to baseline models (strat_label is NA = no country/year
  # stratification); drop cross-tested rows, which ship in their own files.
  df <- data |>
    dplyr::filter(normalize_species(.data$species) %in% normalize_species(bug)) |>
    dplyr::filter(.data$feature_type %in% model_scale) |>
    dplyr::filter(.data$feature_subtype %in% data_type) |>
    dplyr::filter(is.na(.data$strat_label) | !nzchar(.data$strat_label))
  if ("cross_test" %in% names(df)) {
    df <- dplyr::filter(df, !.data$cross_test)
  }

  # Filter by drug class or drug if not "all"
  if (!is.null(amr_drug_class) && length(amr_drug_class) > 0 &&
    !identical(amr_drug_class, "all")) {
    df <- df |>
      dplyr::filter(
        (.data$drug_label == "drug_class" & .data$drug_or_class %in% amr_drug_class) |
          (.data$drug_label == "drug" & .data$drug_or_class %in% amr_drug)
      )
  }

  if (!nrow(df)) {
    return(plotly::plot_ly() |>
      plotly::layout(title = list(text = "No data for current selection", x = 0)))
  }

  # Normalize species and set ESKAPE factor order
  eskape_order <- c("Efa", "Sau", "Kpn", "Aba", "Pae", "Esp")
  df <- df |>
    dplyr::mutate(species = normalize_species(.data$species))
  present <- unique(df$species)
  lvls <- intersect(eskape_order, present)
  if (length(lvls) > 0) {
    df <- df |> dplyr::mutate(species = factor(.data$species, levels = lvls))
  }

  # Order feature_type so legend + colors line up with SCALE_COLORS.
  ft_levels <- intersect(names(SCALE_COLORS), unique(df$feature_type))
  df <- df |>
    dplyr::mutate(feature_type = factor(.data$feature_type, levels = ft_levels))

  has_drug <- !is.null(amr_drug) && length(amr_drug) > 0 && nzchar(amr_drug[1])

  # Build the plot one feature_type at a time. Both the main box trace and the
  # overlay-marker trace (a transparent box with boxpoints = "all") share the
  # same offsetgroup so boxmode = "group" puts them in the same dodge slot.
  p <- plotly::plot_ly() |>
    plotly::layout(
      title = list(text = "Performance metrics", x = 0, font = list(size = 14, color = "#333333")),
      boxmode = "group",
      xaxis = list(
        title = list(text = "Species", font = list(size = 14, color = "#333333")),
        tickfont = list(size = 12, color = "#333333"), tickangle = -45
      ),
      yaxis = list(
        title = list(text = metrics, font = list(size = 14, color = "#333333")),
        tickfont = list(size = 12, color = "#333333"),
        # MCC spans [-1, 1] (0 = random); bal_acc / f1 stay in [0, 1].
        range = if (identical(metrics, "mcc")) c(-1.03, 1.03) else c(0, 1.03)
      ),
      legend = list(
        title = list(text = "Scale", font = list(size = 12, color = "#333333")),
        font = list(size = 10)
      )
    )

  for (ft in ft_levels) {
    sub <- df |> dplyr::filter(.data$feature_type == ft)
    if (!nrow(sub)) next
    color <- unname(SCALE_COLORS[ft])
    fill <- paste0(color, "66") # ~40% opacity, matches old alpha = 0.4

    p <- p |>
      plotly::add_trace(
        type = "box",
        x = sub$species,
        y = sub[[metrics]],
        name = ft,
        legendgroup = ft,
        offsetgroup = ft,
        line = list(color = color),
        fillcolor = fill,
        marker = list(color = color),
        boxpoints = FALSE
      )

    if (has_drug) {
      sub_pts <- sub |> dplyr::filter(.data$drug_or_class %in% amr_drug)
      if (nrow(sub_pts) > 0) {
        p <- p |>
          plotly::add_trace(
            type = "box",
            x = sub_pts$species,
            y = sub_pts[[metrics]],
            boxpoints = "all",
            pointpos = 0,
            jitter = 0,
            line = list(width = 0),
            fillcolor = "rgba(0,0,0,0)",
            marker = list(size = 12, color = color),
            legendgroup = ft,
            offsetgroup = ft,
            showlegend = FALSE,
            hoverinfo = "x+y+name",
            name = ft
          )
      }
    }
  }

  p
}


#' Shared scale-label display map and matching colors
#'
#' @return A list with `labels` (feature_type -> display label), `order` (the
#'   feature_type keys in display order), and `colors` (SCALE_COLORS keyed by
#'   display label).
#' @keywords internal
#' @noRd
.scale_label_map <- function() {
  labels <- c(
    domains = "Domain", genes = "Gene",
    proteins = "Protein", struct = "Struct",
    args = "ARG", cogs = "COG"
  )
  list(
    labels = labels,
    order = names(labels),
    colors = setNames(unname(SCALE_COLORS[names(labels)]), labels)
  )
}


#' Prepare baseline MCC data for the Performance overview plots
#'
#' Drops stratified / cross-test rows so only baseline models remain, and adds a
#' `species_display` column (from `species_label` when present).
#'
#' @param data Performance tibble from loadMLResults().
#' @return A filtered tibble with `species_display`, or NULL when empty.
#' @keywords internal
#' @noRd
.prep_mcc_data <- function(data) {
  if (is.null(data) || !is.data.frame(data) || !nrow(data)) {
    return(NULL)
  }
  df <- data |>
    dplyr::filter(is.na(.data$strat_label) | !nzchar(.data$strat_label)) |>
    dplyr::filter(!.data$cross_test, !is.na(.data$mcc))
  if (!nrow(df)) {
    return(NULL)
  }
  if ("species_label" %in% names(df)) {
    df <- df |>
      dplyr::mutate(species_display = gsub("_", " ", .data$species_label))
  } else {
    df <- df |>
      dplyr::mutate(species_display = .data$species)
  }
  df
}


#' Facetted MCC strip plot (Performance overview)
#'
#' MCC distribution per species and molecular scale, highlighting the selected
#' drug or drug class via point alpha + size. Baseline rows only (via
#' .prep_mcc_data()).
#'
#' @param data Performance tibble from loadMLResults().
#' @param selected_drug_class Drug class to highlight, "all"/NULL for none.
#' @param selected_drug Drug to highlight (takes priority over the class).
#' @return A plotly figure (empty placeholder when there is no matching data).
#' @keywords internal
#' @noRd
makeMCCStripPlot <- function(data, selected_drug_class = NULL,
                             selected_drug = NULL) {
  df <- .prep_mcc_data(data)
  if (is.null(df)) {
    return(plotly::plot_ly() |> plotly::layout(title = "No data available"))
  }

  sl <- .scale_label_map()
  scale_labels <- sl$labels
  scale_order <- sl$order
  scale_colors <- sl$colors

  df <- df |>
    dplyr::filter(.data$feature_type %in% scale_order) |>
    dplyr::mutate(
      scale_label = factor(
        scale_labels[.data$feature_type],
        levels = unname(scale_labels)
      )
    )

  if (!nrow(df)) {
    return(plotly::plot_ly() |>
      plotly::layout(title = "No data for selected filters"))
  }

  # Highlight selected drug or drug class. Priority: specific drug > class.
  # When drug_class == "all", suppress both highlights - no specific class
  # is selected so the accompanying drug selection is not meaningful.
  class_active <- !is.null(selected_drug_class) &&
    nzchar(selected_drug_class) &&
    selected_drug_class != "all"
  use_drug <- class_active && !is.null(selected_drug) && nzchar(selected_drug)
  use_class <- class_active && !use_drug

  highlighted <- if (use_drug) {
    df$drug_or_class == selected_drug & df$drug_label == "drug"
  } else if (use_class) {
    df$drug_or_class == selected_drug_class &
      df$drug_label == "drug_class"
  } else {
    rep(FALSE, nrow(df))
  }
  df <- df |> dplyr::mutate(highlighted = highlighted)

  any_highlighted <- use_drug || use_class

  df <- df |>
    dplyr::mutate(
      pt_alpha = dplyr::if_else(.data$highlighted & any_highlighted, 0.9,
        dplyr::if_else(any_highlighted, 0.15, 0.7)
      ),
      pt_size = dplyr::if_else(.data$highlighted & any_highlighted, 3.5,
        dplyr::if_else(any_highlighted, 1.5, 2.5)
      )
    )

  g <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data$species_display,
      y = .data$mcc,
      color = .data$scale_label,
      fill = .data$scale_label,
      alpha = .data$pt_alpha,
      size = .data$pt_size,
      text = paste0(
        "Drug/class: ", .data$drug_or_class,
        "\nMCC: ", round(.data$mcc, 3),
        "\nEncoding: ", .data$feature_subtype
      )
    )
  ) +
    ggplot2::geom_boxplot(
      alpha = 0.3, outlier.shape = NA,
      width = 0.5, linewidth = 0.4
    ) +
    ggplot2::geom_jitter(width = 0.15) +
    ggplot2::geom_hline(
      yintercept = 0, linetype = "dashed",
      color = "gray50", linewidth = 0.4
    ) +
    ggplot2::facet_grid(scale_label ~ ., switch = "y") +
    ggplot2::scale_color_manual(values = scale_colors) +
    ggplot2::scale_fill_manual(values = scale_colors) +
    ggplot2::scale_alpha_identity() +
    ggplot2::scale_size_identity() +
    ggplot2::coord_cartesian(ylim = c(-1.05, 1.05)) +
    ggplot2::labs(x = NULL, y = "MCC") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      legend.position = "none",
      strip.placement = "outside",
      strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1, color = "#333333"),
      panel.spacing = ggplot2::unit(0.3, "lines"),
      axis.text = ggplot2::element_text(size = 11, color = "#333333"),
      axis.text.x = ggplot2::element_text(angle = 40, hjust = 1, face = "italic")
    )

  plotly::ggplotly(g, tooltip = "text") |>
    plotly::layout(
      title = list(
        text = "MCC per species and molecular scale", x = 0,
        font = list(size = 13, color = "#333333", family = "Arial, sans-serif")
      ),
      margin = list(t = 50)
    )
}


#' Three-panel MCC heatmap (Performance overview)
#'
#' Three sections sharing the drug-class y-axis, showing median MCC by species,
#' molecular scale, and data encoding. Highlights the selected drug-class row
#' across all three sections.
#'
#' @param data Performance tibble from loadMLResults().
#' @param selected_drug_class Drug class row to highlight, "all"/NULL for none.
#' @return A 3-panel plotly subplot (empty placeholder when there is no
#'   drug-class data).
#' @keywords internal
#' @noRd
makeMCCHeatmap <- function(data, selected_drug_class = NULL) {
  df <- .prep_mcc_data(data)
  if (is.null(df)) {
    return(plotly::plot_ly() |> plotly::layout(title = "No data available"))
  }
  df <- df |> dplyr::filter(.data$drug_label == "drug_class")
  if (!nrow(df)) {
    return(plotly::plot_ly() |>
      plotly::layout(title = "No drug-class data available"))
  }

  # Drug-class order: most-represented first (bottom of y = first row)
  drug_order <- df |>
    dplyr::count(.data$drug_or_class) |>
    dplyr::arrange(dplyr::desc(.data$n)) |>
    dplyr::pull(.data$drug_or_class)
  drug_order_rev <- rev(drug_order)

  df <- df |>
    dplyr::mutate(
      drug_or_class = factor(.data$drug_or_class, levels = drug_order_rev)
    )

  # Section 1: species x drug_class (grayscale median MCC)
  spp_summ <- df |>
    dplyr::group_by(.data$drug_or_class, .data$species_display) |>
    dplyr::summarise(
      med_mcc = median(.data$mcc, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      drug_or_class = factor(.data$drug_or_class, levels = drug_order_rev)
    )

  g1 <- ggplot2::ggplot(
    spp_summ,
    ggplot2::aes(
      x = .data$species_display,
      y = .data$drug_or_class,
      fill = .data$med_mcc,
      text = paste0(
        "Species: ", .data$species_display,
        "\nDrug class: ", .data$drug_or_class,
        "\nMedian MCC: ", round(.data$med_mcc, 3)
      )
    )
  ) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient(
      low = "#f7fbff", high = "#08306b",
      limits = c(0, 1.0), name = "MCC", na.value = "white"
    ) +
    ggplot2::labs(x = NULL, y = "Drug class", title = "Species") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5, size = 12,
        color = "#333333", family = "sans"
      ),
      axis.text.x = ggplot2::element_text(angle = 40, hjust = 1, color = "#333333", face = "italic"),
      axis.text.y = ggplot2::element_text(color = "#333333"),
      legend.position = "bottom"
    )

  # Section 2: molecular scale x drug_class (alpha-modulated scale color)
  sl <- .scale_label_map()
  scale_labels <- sl$labels
  scale_order <- sl$order
  scale_colors <- sl$colors

  sc_summ <- df |>
    dplyr::filter(.data$feature_type %in% scale_order) |>
    dplyr::group_by(.data$drug_or_class, .data$feature_type) |>
    dplyr::summarise(
      med_mcc = median(.data$mcc, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      scale_label = factor(
        scale_labels[.data$feature_type],
        levels = unname(scale_labels)
      ),
      drug_or_class = factor(.data$drug_or_class, levels = drug_order_rev),
      alpha_val = pmax(0, pmin(1, .data$med_mcc))
    )

  g2 <- ggplot2::ggplot(
    sc_summ,
    ggplot2::aes(
      x = .data$scale_label,
      y = .data$drug_or_class,
      fill = .data$scale_label,
      alpha = .data$alpha_val,
      text = paste0(
        "Scale: ", .data$scale_label,
        "\nDrug class: ", .data$drug_or_class,
        "\nMedian MCC: ", round(.data$med_mcc, 3)
      )
    )
  ) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_manual(values = scale_colors, guide = "none") +
    ggplot2::scale_alpha_continuous(range = c(0.15, 1.0), guide = "none") +
    ggplot2::labs(x = NULL, y = NULL, title = "Molecular scale") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5, size = 12,
        color = "#333333", family = "sans"
      ),
      axis.text.x = ggplot2::element_text(angle = 40, hjust = 1, colour = "#333333"),
      axis.text.y = ggplot2::element_blank()
    )

  # Section 3: data encoding x drug_class
  subtype_colors <- c(Binary = "#6495ED", Counts = "#BA55D3")

  st_summ <- df |>
    dplyr::group_by(.data$drug_or_class, .data$feature_subtype) |>
    dplyr::summarise(
      med_mcc = median(.data$mcc, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      subtype_label = dplyr::case_when(
        .data$feature_subtype == "binary" ~ "Binary",
        .data$feature_subtype == "counts" ~ "Counts",
        TRUE ~ .data$feature_subtype
      ),
      drug_or_class = factor(.data$drug_or_class, levels = drug_order_rev),
      alpha_val = pmax(0, pmin(1, .data$med_mcc))
    )

  g3 <- ggplot2::ggplot(
    st_summ,
    ggplot2::aes(
      x = .data$subtype_label,
      y = .data$drug_or_class,
      fill = .data$subtype_label,
      alpha = .data$alpha_val,
      text = paste0(
        "Encoding: ", .data$subtype_label,
        "\nDrug class: ", .data$drug_or_class,
        "\nMedian MCC: ", round(.data$med_mcc, 3)
      )
    )
  ) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_manual(values = subtype_colors, guide = "none") +
    ggplot2::scale_alpha_continuous(range = c(0.15, 1.0), guide = "none") +
    ggplot2::labs(x = NULL, y = NULL, title = "Data type") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5, size = 12,
        color = "#333333", family = "sans"
      ),
      axis.text.x = ggplot2::element_text(angle = 40, hjust = 1, colour = "#333333"),
      axis.text.y = ggplot2::element_blank()
    )

  # Highlight the selected drug-class row across all three panels.
  highlight_y <- if (
    !is.null(selected_drug_class) &&
      nzchar(selected_drug_class) &&
      selected_drug_class != "all" &&
      selected_drug_class %in% drug_order_rev
  ) {
    which(drug_order_rev == selected_drug_class)
  } else {
    NULL
  }
  if (!is.null(highlight_y)) {
    highlight_rect <- ggplot2::annotate(
      "rect",
      xmin = -Inf, xmax = Inf,
      ymin = highlight_y - 0.5, ymax = highlight_y + 0.5,
      fill = NA, color = "#FFD700", linewidth = 1.5
    )
    g1 <- g1 + highlight_rect
    g2 <- g2 + highlight_rect
    g3 <- g3 + highlight_rect
  }

  p1 <- plotly::ggplotly(g1, tooltip = "text")
  p2 <- plotly::ggplotly(g2, tooltip = "text")
  p3 <- plotly::ggplotly(g3, tooltip = "text")

  plotly::subplot(p1, p2, p3,
    nrows  = 1,
    shareY = TRUE,
    widths = c(0.45, 0.35, 0.20),
    margin = 0.03
  ) |>
    plotly::layout(
      title = list(
        text = "MCC by drug class, species, molecular scale, and data type",
        x = 0,
        font = list(size = 13, color = "#333333", family = "Arial, sans-serif")
      ),
      legend = list(orientation = "h", y = -0.15)
    )
}
