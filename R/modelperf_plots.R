# Model performance tab visualizations.
#
# Each plot has a table function that does all the filtering and summarizing
# (modelPerfTable(), mccStripTable(), mccHeatmapTable()) and a plot function
# that only draws that table. The tables are what another front end would
# reuse.


# ---- Shared helpers ----

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
    args = "ARG", cogs = "COG", DefenseCas = "Defense/Cas"
  )
  list(
    labels = labels,
    order = names(labels),
    colors = setNames(unname(SCALE_COLORS[names(labels)]), labels)
  )
}


#' Molecular-scale choices for a selector
#'
#' @param present feature_type values to offer, or NULL for every known scale.
#' @return A named character vector (names = display labels, values =
#'   feature_type), in display order.
#' @keywords internal
#' @noRd
.scale_choices <- function(present = NULL) {
  labels <- .scale_label_map()$labels
  if (!is.null(present)) labels <- labels[names(labels) %in% present]
  stats::setNames(names(labels), unname(labels))
}


#' Give a one-color heatmap a valid color scale
#'
#' When every tile has the same value (for example a single drug class and
#' species), ggplotly() writes a color scale with one stop, which plotly
#' ignores and draws in gray. Repeats that stop at both ends instead.
#'
#' @param p A plotly figure from ggplotly().
#' @return The built figure, with each heatmap's color scale fixed.
#' @keywords internal
#' @noRd
.fix_single_value_heatmap <- function(p) {
  p <- plotly::plotly_build(p)
  for (i in seq_along(p$x$data)) {
    trace <- p$x$data[[i]]
    if (!identical(trace$type, "heatmap")) next
    colors <- unlist(trace$colorscale)
    colors <- colors[grepl("^#|^rgb", colors)]
    if (length(colors) == 1) {
      p$x$data[[i]]$colorscale <- list(list(0, colors), list(1, colors))
    }
  }
  p
}


#' Lower end of a metric axis
#'
#' @param values Metric values being plotted.
#' @return -0.05, or a little below the lowest value when it is below 0, so the
#'   0 line stays in view without leaving most of the axis empty.
#' @keywords internal
#' @noRd
.axis_floor <- function(values) {
  min(0, values, na.rm = TRUE) - 0.05
}


#' Species display order: ESKAPE species first, then the rest alphabetically
#'
#' @param species Species codes, one per row.
#' @param species_display Display names, one per row.
#' @return A factor of `species_display` with that level order.
#' @keywords internal
#' @noRd
.species_display_factor <- function(species, species_display) {
  eskape_order <- c("Efa", "Sau", "Kpn", "Aba", "Pae", "Esp")
  pairs <- unique(data.frame(species = species, display = species_display))
  rank <- match(pairs$species, eskape_order, nomatch = length(eskape_order) + 1)
  levels <- unique(pairs$display[order(rank, pairs$display)])
  factor(species_display, levels = levels)
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
    dplyr::filter(!.data$cross_test %in% TRUE, !is.na(.data$mcc))
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
  df |>
    dplyr::mutate(species_display = .species_display_factor(
      normalize_species(.data$species), .data$species_display
    ))
}


# ---- Tables ----

#' Rows behind the model performance box plot
#'
#' Baseline models (not stratified or cross-tested) for the selected species,
#' feature types and encoding. With a drug class selected, keeps that class's
#' models and the selected drugs' models.
#'
#' @param data Performance tibble from loadMLResults() / queryData().
#' @param bug Species code(s) to include.
#' @param model_scale Molecular scale(s) (feature_type) to include.
#' @param data_type Data encoding(s) (feature_subtype) to include.
#' @param metrics Name of the metric column to report as `value`.
#' @param amr_drug_class Selected drug class(es), or "all" for no filter.
#' @param amr_drug Selected drug(s), marked in `highlighted`.
#' @return A tibble with one row per model: `species`, `species_display`
#'   (factor, ESKAPE species first), `feature_type` (factor, in SCALE_COLORS
#'   order), `feature_subtype`, `drug_label`, `drug_or_class`, `value` and
#'   `highlighted`. Zero rows when nothing matches.
#' @keywords internal
#' @noRd
modelPerfTable <- function(data, bug, model_scale, data_type, metrics,
                           amr_drug_class = "all", amr_drug = NULL) {
  empty <- tibble::tibble(
    species = character(0), species_display = factor(),
    feature_type = factor(), feature_subtype = character(0),
    drug_label = character(0), drug_or_class = character(0),
    value = numeric(0), highlighted = logical(0)
  )
  if (is.null(data) || !is.data.frame(data) || !nrow(data) ||
    !metrics %in% names(data)) {
    return(empty)
  }

  # Baseline models only (strat_label NA, not cross-tested), without NA
  # metrics so plotly's box trace doesn't warn about dropped observations.
  df <- data |>
    dplyr::filter(normalize_species(.data$species) %in% normalize_species(bug)) |>
    dplyr::filter(.data$feature_type %in% model_scale) |>
    dplyr::filter(.data$feature_subtype %in% data_type) |>
    dplyr::filter(is.na(.data$strat_label) | !nzchar(.data$strat_label)) |>
    dplyr::filter(!is.na(.data[[metrics]]))
  if ("cross_test" %in% names(df)) {
    df <- dplyr::filter(df, !.data$cross_test %in% TRUE)
  }

  if (!is.null(amr_drug_class) && length(amr_drug_class) > 0 &&
    !identical(amr_drug_class, "all")) {
    df <- df |>
      dplyr::filter(
        (.data$drug_label == "drug_class" & .data$drug_or_class %in% amr_drug_class) |
          (.data$drug_label == "drug" & .data$drug_or_class %in% amr_drug)
      )
  }
  if (!nrow(df)) {
    return(empty)
  }

  # Prefer the full species name (species_label) over the 3-letter code.
  species <- normalize_species(df$species)
  display <- if ("species_label" %in% names(df)) {
    gsub("_", " ", df$species_label)
  } else {
    species
  }
  ft_levels <- intersect(names(SCALE_COLORS), unique(df$feature_type))

  tibble::tibble(
    species = species,
    species_display = .species_display_factor(species, display),
    feature_type = factor(df$feature_type, levels = ft_levels),
    feature_subtype = df$feature_subtype,
    drug_label = df$drug_label,
    drug_or_class = df$drug_or_class,
    value = df[[metrics]],
    highlighted = df$drug_or_class %in% amr_drug
  )
}


#' Rows behind the MCC strip plot
#'
#' Baseline models with a known scale, marking the selected drug or, when no
#' drug is selected, the selected drug class. Nothing is marked when the class
#' is "all".
#'
#' @param data Performance tibble from loadMLResults().
#' @param selected_drug_class Drug class to mark, "all"/NULL for none.
#' @param selected_drug Drug to mark (takes priority over the class).
#' @return A tibble with `species_display` (factor), `scale_label` (factor),
#'   `drug_label`, `drug_or_class`, `feature_subtype`, `mcc` and
#'   `highlighted`, or NULL when there is no baseline data.
#' @keywords internal
#' @noRd
mccStripTable <- function(data, selected_drug_class = NULL,
                          selected_drug = NULL) {
  df <- .prep_mcc_data(data)
  if (is.null(df)) {
    return(NULL)
  }
  sl <- .scale_label_map()
  df <- df |>
    dplyr::filter(.data$feature_type %in% sl$order) |>
    dplyr::mutate(
      scale_label = factor(
        sl$labels[.data$feature_type],
        levels = intersect(unname(sl$labels), sl$labels[.data$feature_type])
      )
    )

  # Priority: specific drug > class. With class "all" the drug selection is
  # not meaningful, so nothing is marked.
  class_active <- !is.null(selected_drug_class) &&
    nzchar(selected_drug_class) &&
    selected_drug_class != "all"
  use_drug <- class_active && !is.null(selected_drug) && nzchar(selected_drug)
  use_class <- class_active && !use_drug
  highlighted <- if (use_drug) {
    df$drug_or_class == selected_drug & df$drug_label == "drug"
  } else if (use_class) {
    df$drug_or_class == selected_drug_class & df$drug_label == "drug_class"
  } else {
    rep(FALSE, nrow(df))
  }

  tibble::tibble(
    species_display = df$species_display,
    scale_label = df$scale_label,
    drug_label = df$drug_label,
    drug_or_class = df$drug_or_class,
    feature_subtype = df$feature_subtype,
    mcc = df$mcc,
    highlighted = highlighted
  )
}


#' Median MCC per drug class, by species, scale and data type
#'
#' The three heatmap panels as one long table. Drug-class models only.
#'
#' @param data Performance tibble from loadMLResults().
#' @return A tibble with `panel` ("Species", "Molecular scale" or "Data
#'   type"), `column` (the species, scale or data type), `drug_or_class`
#'   (factor, most-modeled class last), `median_mcc` and `n_models`, or NULL
#'   when there is no drug-class data.
#' @keywords internal
#' @noRd
mccHeatmapTable <- function(data) {
  df <- .prep_mcc_data(data)
  if (is.null(df)) {
    return(NULL)
  }
  df <- df |> dplyr::filter(.data$drug_label == "drug_class")
  if (!nrow(df)) {
    return(NULL)
  }

  # Most-modeled drug class first, shown at the top of the y-axis.
  drug_order <- df |>
    dplyr::count(.data$drug_or_class) |>
    dplyr::arrange(dplyr::desc(.data$n)) |>
    dplyr::pull(.data$drug_or_class)

  sl <- .scale_label_map()
  subtype_labels <- c(binary = "Binary", counts = "Counts")
  by_column <- function(df, panel, column) {
    df |>
      dplyr::mutate(column = column) |>
      dplyr::group_by(.data$drug_or_class, .data$column) |>
      dplyr::summarise(
        median_mcc = stats::median(.data$mcc, na.rm = TRUE),
        n_models = dplyr::n(),
        .groups = "drop"
      ) |>
      dplyr::mutate(panel = panel)
  }
  scales <- df |> dplyr::filter(.data$feature_type %in% sl$order)
  subtypes <- dplyr::coalesce(
    unname(subtype_labels[df$feature_subtype]), df$feature_subtype
  )

  dplyr::bind_rows(
    by_column(df, "Species", as.character(df$species_display)),
    by_column(scales, "Molecular scale", unname(sl$labels[scales$feature_type])),
    by_column(df, "Data type", subtypes)
  ) |>
    dplyr::mutate(
      drug_or_class = factor(.data$drug_or_class, levels = rev(drug_order))
    ) |>
    dplyr::select("panel", "column", "drug_or_class", "median_mcc", "n_models")
}


# ---- Plots ----

#' Baseline ML model performance plot
#'
#' Grouped box plots of a performance metric by species and molecular scale,
#' drawn from modelPerfTable(). Points mark the selected drugs.
#'
#' @inheritParams modelPerfTable
#' @return A plotly box-plot figure (an empty placeholder when there is no
#'   matching data).
#' @keywords internal
#' @noRd
makeModelPerformancePlot <- function(
  data, bug, model_scale, data_type, metrics,
  amr_drug_class, amr_drug
) {
  if (is.null(data) || !is.data.frame(data) || !nrow(data)) {
    return(plotly_placeholder("No data available"))
  }
  df <- modelPerfTable(
    data, bug, model_scale, data_type, metrics, amr_drug_class, amr_drug
  )
  if (!nrow(df)) {
    return(plotly_placeholder("No data for current selection"))
  }

  sl <- .scale_label_map()
  metric_titles <- c(mcc = "MCC", bal_acc = "Balanced accuracy", f1 = "F1")

  # Build the plot one feature_type at a time. Both the main box trace and the
  # overlay-marker trace (a transparent box with boxpoints = "all") share the
  # same offsetgroup so boxmode = "group" puts them in the same dodge slot.
  p <- plotly::plot_ly() |>
    plotly::layout(
      title = list(text = "Performance metrics", x = 0, font = list(size = 14)),
      boxmode = "group",
      xaxis = list(
        title = list(text = "Species", font = list(size = 10)),
        tickfont = list(size = 10), tickangle = -45
      ),
      yaxis = list(
        title = list(
          text = unname(metric_titles[metrics]) %||% metrics,
          font = list(size = 10)
        ),
        tickfont = list(size = 10),
        # From 0 (chance for MCC) to 1, extended below 0 only when a model
        # scores lower, so the models aren't squeezed into the top of [-1, 1].
        range = c(.axis_floor(df$value), 1.03)
      ),
      legend = list(
        title = list(text = "Scale", font = list(size = 12)),
        font = list(size = 10)
      )
    )

  for (ft in levels(df$feature_type)) {
    sub <- df |> dplyr::filter(.data$feature_type == ft)
    if (!nrow(sub)) next
    color <- unname(SCALE_COLORS[ft])
    fill <- paste0(color, "66") # ~40% opacity, matches old alpha = 0.4
    label <- unname(sl$labels[ft]) %||% ft

    p <- p |>
      plotly::add_trace(
        type = "box",
        x = sub$species_display,
        y = sub$value,
        name = label,
        legendgroup = ft,
        offsetgroup = ft,
        line = list(color = color),
        fillcolor = fill,
        marker = list(color = color),
        boxpoints = FALSE
      )

    sub_pts <- sub |> dplyr::filter(.data$highlighted)
    if (nrow(sub_pts) > 0) {
      p <- p |>
        plotly::add_trace(
          type = "box",
          x = sub_pts$species_display,
          y = sub_pts$value,
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
          name = label
        )
    }
  }

  p
}


#' Facetted MCC strip plot (Performance overview)
#'
#' MCC distribution per species and molecular scale, drawn from
#' mccStripTable(); the marked drug or drug class is drawn larger and the rest
#' faded.
#'
#' @inheritParams mccStripTable
#' @return A plotly figure (empty placeholder when there is no matching data).
#' @keywords internal
#' @noRd
makeMCCStripPlot <- function(data, selected_drug_class = NULL,
                             selected_drug = NULL) {
  df <- mccStripTable(data, selected_drug_class, selected_drug)
  if (is.null(df)) {
    return(plotly_placeholder("No data available"))
  }
  if (!nrow(df)) {
    return(plotly_placeholder("No data for selected filters"))
  }

  any_highlighted <- any(df$highlighted)
  df <- df |>
    dplyr::mutate(
      pt_alpha = dplyr::if_else(.data$highlighted, 0.9,
        dplyr::if_else(any_highlighted, 0.15, 0.7)
      ),
      pt_size = dplyr::if_else(.data$highlighted, 3.5,
        dplyr::if_else(any_highlighted, 1.5, 2.5)
      )
    )

  # Keep alpha/size on geom_jitter only. In the top-level aes() they'd also
  # apply to geom_boxplot's line width, triggering ggplot2 3.4's size-for-lines
  # deprecation.
  g <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data$species_display,
      y = .data$mcc,
      color = .data$scale_label,
      fill = .data$scale_label,
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
    ggplot2::geom_jitter(
      ggplot2::aes(alpha = .data$pt_alpha, size = .data$pt_size),
      width = 0.15
    ) +
    ggplot2::geom_hline(
      yintercept = 0, linetype = "dashed",
      color = "gray50", linewidth = 0.4
    ) +
    ggplot2::facet_grid(scale_label ~ ., switch = "y") +
    ggplot2::scale_color_manual(values = .scale_label_map()$colors) +
    ggplot2::scale_fill_manual(values = .scale_label_map()$colors) +
    ggplot2::scale_alpha_identity() +
    ggplot2::scale_size_identity() +
    ggplot2::coord_cartesian(ylim = c(.axis_floor(df$mcc), 1.05)) +
    ggplot2::labs(x = NULL, y = "MCC") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position   = "none",
      strip.placement   = "outside",
      strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1),
      panel.spacing     = ggplot2::unit(0.3, "lines")
    )

  plotly::ggplotly(g, tooltip = "text") |>
    plotly::layout(
      title = list(
        text = "MCC per species and molecular scale", x = 0,
        font = list(size = 13, color = "#333333", family = "Arial, sans-serif")
      ),
      xaxis = list(
        tickfont = list(size = 10), tickangle = -45
      ),
      margin = list(t = 50)
    )
}


#' Three-panel MCC heatmap (Performance overview)
#'
#' Three sections sharing the drug-class y-axis, showing median MCC by species,
#' molecular scale, and data encoding, drawn from mccHeatmapTable(). Highlights
#' the selected drug-class row across all three sections.
#'
#' @param data Performance tibble from loadMLResults().
#' @param selected_drug_class Drug class row to highlight, "all"/NULL for none.
#' @return A 3-panel plotly subplot (empty placeholder when there is no
#'   drug-class data).
#' @keywords internal
#' @noRd
makeMCCHeatmap <- function(data, selected_drug_class = NULL) {
  if (is.null(.prep_mcc_data(data))) {
    return(plotly_placeholder("No data available"))
  }
  tbl <- mccHeatmapTable(data)
  if (is.null(tbl)) {
    return(plotly_placeholder("No drug-class data available"))
  }
  drug_levels <- levels(tbl$drug_or_class)
  panel <- function(name) {
    tbl |>
      dplyr::filter(.data$panel == name) |>
      dplyr::mutate(
        alpha_val = pmax(0, pmin(1, .data$median_mcc)),
        text = paste0(
          name, ": ", .data$column,
          "\nDrug class: ", .data$drug_or_class,
          "\nMedian MCC: ", round(.data$median_mcc, 3),
          "\nModels: ", .data$n_models
        )
      )
  }
  panel_theme <- function(show_y) {
    ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          hjust = 0.5, size = 12,
          color = "#333333", family = "sans"
        ),
        axis.text.x = ggplot2::element_text(angle = 40, hjust = 1),
        axis.text.y = if (show_y) NULL else ggplot2::element_blank(),
        legend.position = "bottom"
      )
  }

  # Section 1: species x drug_class (grayscale median MCC)
  g1 <- ggplot2::ggplot(
    panel("Species"),
    ggplot2::aes(
      x = .data$column, y = .data$drug_or_class,
      fill = .data$median_mcc, text = .data$text
    )
  ) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient(
      low = "#f7fbff", high = "#08306b",
      limits = c(0, 1.0), name = "MCC", na.value = "white"
    ) +
    ggplot2::labs(x = NULL, y = "Drug class", title = "Species") +
    panel_theme(show_y = TRUE)

  # Section 2: molecular scale x drug_class (alpha-modulated scale color)
  sl <- .scale_label_map()
  scales <- panel("Molecular scale") |>
    dplyr::mutate(column = factor(
      .data$column,
      levels = intersect(unname(sl$labels), .data$column)
    ))
  g2 <- ggplot2::ggplot(
    scales,
    ggplot2::aes(
      x = .data$column, y = .data$drug_or_class,
      fill = .data$column, alpha = .data$alpha_val, text = .data$text
    )
  ) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_manual(values = sl$colors, guide = "none") +
    ggplot2::scale_alpha_continuous(range = c(0.15, 1.0), guide = "none") +
    ggplot2::labs(x = NULL, y = NULL, title = "Molecular scale") +
    panel_theme(show_y = FALSE)

  # Section 3: data encoding x drug_class
  subtype_colors <- c(Binary = "#6495ED", Counts = "#BA55D3")
  g3 <- ggplot2::ggplot(
    panel("Data type"),
    ggplot2::aes(
      x = .data$column, y = .data$drug_or_class,
      fill = .data$column, alpha = .data$alpha_val, text = .data$text
    )
  ) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_manual(values = subtype_colors, guide = "none") +
    ggplot2::scale_alpha_continuous(range = c(0.15, 1.0), guide = "none") +
    ggplot2::labs(x = NULL, y = NULL, title = "Data type") +
    panel_theme(show_y = FALSE)

  # Highlight the selected drug-class row across all three panels.
  highlight_y <- if (
    !is.null(selected_drug_class) &&
      nzchar(selected_drug_class) &&
      selected_drug_class != "all" &&
      selected_drug_class %in% drug_levels
  ) {
    which(drug_levels == selected_drug_class)
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

  p1 <- .fix_single_value_heatmap(plotly::ggplotly(g1, tooltip = "text"))
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
