# Cross-model holdout visualisations.


# makeCrossModelFeatureImportancePlot: heatmap of top features for holdout models.
# top_data: pre-loaded top-features tibble (country or year stratified rows).
# amRml column mapping:
#   drug_or_class -> drug/class abbreviation
#   strat_label   -> "country" or "year"
#   strat_value   -> trained-on country/year
makeCrossModelFeatureImportancePlot <- function(
  top_data, bug, drug, cross_model, top_n_features,
  annotated_dir = NULL
) {
  if (is.null(top_data) || !is.data.frame(top_data) || !nrow(top_data)) {
    return(NULL)
  }

  top_n_features <- ifelse(
    identical(top_n_features, "all"), "all",
    as.numeric(top_n_features)
  )

  strat <- if (cross_model == "country") "country" else "year"
  strat_col <- "strat_value" # trained-on value in new schema

  features_df <- top_data |>
    dplyr::filter(.data$species %in% bug) |>
    dplyr::filter(.data$drug_or_class %in% drug) |>
    dplyr::filter(.data$strat_label == strat) |>
    dplyr::filter(!.data$cross_test)

  if (!nrow(features_df)) {
    return(NULL)
  }

  if (!identical(top_n_features, "all")) {
    features_df <- features_df |>
      dplyr::group_by(!!rlang::sym(strat_col)) |>
      dplyr::slice_max(order_by = .data$Importance, n = top_n_features) |>
      dplyr::ungroup()
  }

  vi_wider <- features_df |>
    dplyr::select(.data$Variable, .data$Importance, !!rlang::sym(strat_col)) |>
    tidyr::pivot_wider(
      names_from = strat_col,
      values_from = "Importance",
      values_fn = mean
    )

  if (!nrow(vi_wider)) {
    return(NULL)
  }

  vi_mat <- vi_wider |>
    tibble::column_to_rownames("Variable") |>
    as.matrix()

  # column-wise min-max normalisation
  vi_mat <- apply(vi_mat, 2, function(x) {
    rng <- range(x, na.rm = TRUE)
    if (diff(rng) == 0) {
      return(x)
    }
    (x - rng[1]) / diff(rng)
  })

  max_val <- max(vi_mat, na.rm = TRUE)
  min_val <- min(vi_mat, na.rm = TRUE)

  plotly::plot_ly(
    x = colnames(vi_mat),
    y = rownames(vi_mat),
    z = vi_mat,
    type = "heatmap",
    colorscale = list(c(0, "#c8e8e8"), c(1, "#4e9a9a")),
    colorbar = list(title = "Importance"),
    hovertemplate = paste0(
      "<b>Feature:</b> %{y}<br>",
      "<b>Group:</b> %{x}<br>",
      "<b>Importance:</b> %{z:.3f}<extra></extra>"
    )
  ) |>
    plotly::layout(
      xaxis = list(title = "", tickangle = -45, side = "top"),
      yaxis = list(title = "", autorange = "reversed"),
      margin = list(l = 200, b = 20, t = 80)
    )
}


# makeCrossModelRidgePlot: balanced accuracy by drug class for holdout models,
# coloured by Same (self-eval) vs Different (cross-eval).
# cross_model: "country" or "time"
makeCrossModelRidgePlot <- function(perf_data, bug, cross_model) {
  if (is.null(perf_data) || !is.data.frame(perf_data) || !nrow(perf_data)) {
    return(plotly::plot_ly() |>
      plotly::layout(title = list(text = "No data available", x = 0)))
  }

  strat <- if (cross_model == "country") "country" else "year"

  df <- perf_data |>
    dplyr::filter(
      normalize_species(.data$species) %in% normalize_species(bug)
    ) |>
    dplyr::filter(.data$strat_label == strat) |>
    dplyr::filter(.data$drug_label == "drug_class")

  if (!nrow(df)) {
    return(plotly::plot_ly() |>
      plotly::layout(title = list(text = "No data for selection", x = 0)))
  }

  # Label: Same = trained & tested on same stratum, Different = cross-tested
  df <- df |>
    dplyr::mutate(
      test_type = dplyr::if_else(
        .data$cross_test, "Different", "Same"
      )
    )

  strat_label <- if (cross_model == "country") "Countries" else "Time periods"

  colors <- c(Same = "#d4872a", Different = "#5b8db8")

  p <- plotly::plot_ly()

  for (tt in c("Same", "Different")) {
    sub <- df[df$test_type == tt, ]
    if (!nrow(sub)) next
    col <- colors[[tt]]
    p <- p |>
      plotly::add_trace(
        type = "box",
        x = sub$bal_acc,
        y = sub$drug_or_class,
        name = tt,
        orientation = "h",
        boxpoints = "all",
        jitter = 0.4,
        pointpos = 0,
        marker = list(
          color = col, opacity = 0.7, size = 6
        ),
        line = list(color = col),
        fillcolor = paste0(col, "44"),
        hovertemplate = paste0(
          "<b>Drug class:</b> %{y}<br>",
          "<b>Bal. Accuracy:</b> %{x:.3f}<br>",
          "<b>Test:</b> ", tt, "<extra></extra>"
        )
      )
  }

  p |> plotly::layout(
    title = list(
      text = paste("Balanced accuracy by drug class -", strat_label),
      x = 0, font = list(size = 13)
    ),
    xaxis = list(title = "Balanced accuracy", range = c(0, 1.05)),
    yaxis = list(title = ""),
    legend = list(title = list(text = "Test")),
    boxmode = "group",
    margin = list(l = 100, t = 50)
  )
}


# makeCrossModelPerformancePlot: heatmap of balanced accuracy for holdout models.
# perf_data: pre-loaded performance tibble (country or year stratified rows).
# amRml column mapping:
#   drug_or_class   -> drug/class abbreviation
#   strat_label     -> "country" or "year"
#   strat_value     -> trained-on country/year
#   strat_value_test -> tested-on country/year (NA for self-evaluation)
makeCrossModelPerformancePlot <- function(perf_data, bug, drug, cross_model) {
  if (is.null(perf_data) || !is.data.frame(perf_data) || !nrow(perf_data)) {
    return(NULL)
  }

  strat <- if (cross_model == "country") "country" else "year"

  df <- perf_data |>
    dplyr::filter(normalize_species(.data$species) %in% normalize_species(bug)) |>
    dplyr::filter(.data$drug_or_class %in% drug) |>
    dplyr::filter(.data$strat_label == strat)

  if (!nrow(df)) {
    return(NULL)
  }

  # For self-evaluation rows (strat_value_test is NA), set tested = trained
  df <- df |>
    dplyr::mutate(
      strat_value_test = dplyr::if_else(
        is.na(.data$strat_value_test),
        .data$strat_value,
        .data$strat_value_test
      )
    )

  # Aggregate bal_acc (mean across scales/encodings for same train/test pair)
  models_performance <- df |>
    dplyr::group_by(.data$strat_value, .data$strat_value_test) |>
    dplyr::summarise(bal_acc = mean(.data$bal_acc, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(
      names_from = "strat_value",
      values_from = "bal_acc"
    ) |>
    tibble::column_to_rownames("strat_value_test") |>
    as.matrix()

  if (!length(models_performance)) {
    return(NULL)
  }

  min_val <- min(models_performance, na.rm = TRUE)
  max_val <- max(models_performance, na.rm = TRUE)

  plotly::plot_ly(
    x = colnames(models_performance),
    y = rownames(models_performance),
    z = models_performance,
    type = "heatmap",
    colorscale = list(
      c(0, "#f2f0f7"), c(0.25, "#cbc9e2"),
      c(0.5, "#9e9ac8"), c(0.75, "#756bb1"),
      c(1, "#54278f")
    ),
    zmin = 0.5, zmax = 1.0,
    colorbar = list(title = "Balanced\naccuracy"),
    hovertemplate = paste0(
      "<b>Train data:</b> %{x}<br>",
      "<b>Test data:</b> %{y}<br>",
      "<b>Bal. Accuracy:</b> %{z:.3f}<extra></extra>"
    )
  ) |>
    plotly::layout(
      title = list(
        text = if (cross_model == "country") {
          "Cross-country performance"
        } else {
          "Cross-time performance"
        },
        x = 0, font = list(size = 13)
      ),
      xaxis = list(title = "Train data", side = "bottom"),
      yaxis = list(title = "Test data", autorange = "reversed"),
      margin = list(l = 100, b = 60, t = 50)
    )
}


# .buildCrossDrugMatrix: shared data prep for the cross-drug heatmaps.
# Builds an ordered trained-on (rows) x tested-on (cols) matrix of median nMCC,
# restricted to individual drugs (drug_label == "drug"); drug classes are
# excluded. The diagonal (train == test) is filled from the within-drug
# baseline models. Adapted from amRml::plotCrossDrug() to amRviz's parquet
# schema (ref_drug/test_drug/nmcc in place of drug_or_class/tested_on/mcc).
#
# perf_data: pre-loaded performance tibble (baseline + cross rows) as returned
#   by loadMLResults(). Cross-drug rows have cross_test == TRUE with ref_drug /
#   test_drug populated; baseline rows have cross_test == FALSE, strat_label NA.
# meta: optional tibble mapping drug abbreviation -> class abbreviation, with
#   columns drug_abbr and class_abbr (e.g. from a *_metadata.parquet); used to
#   group rows/columns by class.
# bug: optional species filter (matched via normalize_species()).
#
# Returns list(mat, class_of) with the ordered matrix and a label -> class
# lookup, or NULL when there are no cross-drug rows to plot.
.buildCrossDrugMatrix <- function(perf_data, meta = NULL, bug = NULL) {
  if (is.null(perf_data) || !is.data.frame(perf_data) || !nrow(perf_data)) {
    return(NULL)
  }
  required_cols <- c(
    "cross_test", "ref_drug", "test_drug", "nmcc", "drug_or_class",
    "drug_label", "strat_label"
  )
  if (!all(required_cols %in% names(perf_data))) {
    return(NULL)
  }

  # Restrict to the selected species. Cross-drug rows are stored with
  # species == "cross", so a plain species filter would drop them; resolve the
  # chosen species code to its species_label (the directory name shared by the
  # baseline and cross rows) and filter on that when available.
  df <- perf_data
  if (!is.null(bug)) {
    if ("species_label" %in% names(df)) {
      labs <- df |>
        dplyr::filter(normalize_species(.data$species) %in% normalize_species(bug)) |>
        dplyr::pull(.data$species_label) |>
        unique()
      if (length(labs)) {
        df <- df |> dplyr::filter(.data$species_label %in% labs)
      }
    } else {
      df <- df |>
        dplyr::filter(normalize_species(.data$species) %in% normalize_species(bug))
    }
  }
  if (!nrow(df)) {
    return(NULL)
  }

  # Cross-drug cells: trained on ref_drug, tested on test_drug. Restricted to
  # individual drugs (drug_label == "drug"); drug classes are excluded.
  cross <- df |>
    dplyr::filter(
      .data$cross_test %in% TRUE,
      .data$drug_label == "drug",
      !is.na(.data$ref_drug), !is.na(.data$test_drug)
    ) |>
    dplyr::group_by(trained_on = .data$ref_drug, tested_on = .data$test_drug) |>
    dplyr::summarise(
      median_nmcc = stats::median(.data$nmcc, na.rm = TRUE),
      .groups = "drop"
    )
  if (!nrow(cross)) {
    return(NULL)
  }

  # Self-evaluation diagonal from the within-drug baseline models.
  drugs <- base::union(cross$trained_on, cross$tested_on)
  diag_df <- df |>
    dplyr::filter(
      .data$cross_test %in% FALSE, is.na(.data$strat_label),
      .data$drug_label == "drug",
      .data$drug_or_class %in% drugs
    ) |>
    dplyr::group_by(.data$drug_or_class) |>
    dplyr::summarise(
      median_nmcc = stats::median(.data$nmcc, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::transmute(
      trained_on = .data$drug_or_class,
      tested_on = .data$drug_or_class, .data$median_nmcc
    )

  # Collapse any duplicate train/test pairs so the matrix has one value per cell.
  heatmap_df <- dplyr::bind_rows(cross, diag_df) |>
    dplyr::group_by(.data$trained_on, .data$tested_on) |>
    dplyr::summarise(
      median_nmcc = stats::median(.data$median_nmcc, na.rm = TRUE),
      .groups = "drop"
    )

  mat <- heatmap_df |>
    tidyr::pivot_wider(names_from = "tested_on", values_from = "median_nmcc") |>
    tibble::column_to_rownames("trained_on") |>
    as.matrix()
  if (!length(mat)) {
    return(NULL)
  }

  # Map each drug label to its class for ordering + annotation.
  class_map <- NULL
  if (!is.null(meta) && all(c("drug_abbr", "class_abbr") %in% names(meta))) {
    class_map <- meta |> dplyr::distinct(.data$drug_abbr, .data$class_abbr)
  }
  labels <- base::union(rownames(mat), colnames(mat))
  class_of <- stats::setNames(rep(NA_character_, length(labels)), labels)
  if (!is.null(class_map)) {
    hit <- match(labels, class_map$drug_abbr)
    class_of[!is.na(hit)] <- class_map$class_abbr[hit[!is.na(hit)]]
  }

  # Order rows/columns by class, then alphabetically within class.
  row_order <- rownames(mat)[order(class_of[rownames(mat)], rownames(mat))]
  col_order <- colnames(mat)[order(class_of[colnames(mat)], colnames(mat))]
  mat <- mat[row_order, col_order, drop = FALSE]

  list(mat = mat, class_of = class_of)
}

# .crossDrugColorScale: shared diverging colour scaling for both cross-drug
# heatmaps. Centres an RdBu ramp (low nMCC red -> high blue) on 0.5 (random
# performance) with a symmetric radius derived from the finite matrix values,
# so the static export and the interactive panel are calibrated identically.
# n = number of colour steps. Returns list(colors, lo, hi), or NULL when the
# matrix has no finite values.
.crossDrugColorScale <- function(mat, n = 100) {
  finite_vals <- mat[is.finite(mat)]
  if (!length(finite_vals)) {
    return(NULL)
  }
  rad <- max(abs(range(finite_vals) - 0.5))
  if (!is.finite(rad) || rad == 0) {
    rad <- 0.5
  }
  colors <- grDevices::colorRampPalette(
    RColorBrewer::brewer.pal(11, "RdBu")
  )(n)
  list(colors = colors, lo = 0.5 - rad, hi = 0.5 + rad)
}

# makeCrossDrugHeatmap: static ComplexHeatmap of cross-drug generalization, used
# for the publication-quality PDF export. This is the package's ComplexHeatmap
# (Bioconductor) visualization. Arguments match .buildCrossDrugMatrix(). Returns
# a ComplexHeatmap::Heatmap object, or NULL when there are no cross-drug rows.
makeCrossDrugHeatmap <- function(perf_data, meta = NULL, bug = NULL) {
  prep <- .buildCrossDrugMatrix(perf_data, meta = meta, bug = bug)
  if (is.null(prep)) {
    return(NULL)
  }
  mat <- prep$mat
  class_of <- prep$class_of
  row_order <- rownames(mat)
  col_order <- colnames(mat)

  # Diverging colour scale centred on 0.5 (random performance for nMCC).
  scale <- .crossDrugColorScale(mat)
  if (is.null(scale)) {
    return(NULL)
  }
  col_fun <- circlize::colorRamp2(
    seq(scale$lo, scale$hi, length.out = length(scale$colors)),
    scale$colors
  )

  # Optional class annotations on rows and columns, sharing one colour key.
  left_anno <- NULL
  top_anno <- NULL
  classes <- sort(unique(stats::na.omit(class_of)))
  if (length(classes)) {
    class_colors <- stats::setNames(
      grDevices::hcl.colors(length(classes), palette = "Dark 3"),
      classes
    )
    left_anno <- ComplexHeatmap::rowAnnotation(
      class = class_of[row_order],
      col = list(class = class_colors),
      show_annotation_name = FALSE, show_legend = FALSE, na_col = "grey80"
    )
    top_anno <- ComplexHeatmap::HeatmapAnnotation(
      class = class_of[col_order],
      col = list(class = class_colors),
      show_annotation_name = FALSE, show_legend = TRUE, na_col = "grey80"
    )
  }

  ComplexHeatmap::Heatmap(
    mat,
    name = "median\nnMCC",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_order = row_order,
    column_order = col_order,
    left_annotation = left_anno,
    top_annotation = top_anno,
    show_row_names = TRUE,
    show_column_names = TRUE,
    column_title = "Tested on",
    row_title = "Trained on",
    column_title_side = "bottom",
    row_title_side = "right",
    row_names_gp = grid::gpar(fontsize = 12),
    column_names_gp = grid::gpar(fontsize = 12),
    column_names_rot = 0,
    rect_gp = grid::gpar(col = NA),
    na_col = "grey95",
    show_heatmap_legend = TRUE
  )
}

# makeCrossDrugHeatmapPlotly: interactive plotly rendering of the same cross-drug
# matrix for the dashboard panel (the static ComplexHeatmap above is reserved for
# the export). Rows/columns stay grouped by class; hover shows the train/test
# pair and median nMCC. Arguments match .buildCrossDrugMatrix(). Returns a plotly
# object, or NULL when there are no cross-drug rows.
makeCrossDrugHeatmapPlotly <- function(perf_data, meta = NULL, bug = NULL) {
  prep <- .buildCrossDrugMatrix(perf_data, meta = meta, bug = bug)
  if (is.null(prep)) {
    return(NULL)
  }
  mat <- prep$mat

  # Diverging RdBu colourscale (low nMCC red -> high blue), matching the export.
  scale <- .crossDrugColorScale(mat)
  if (is.null(scale)) {
    return(NULL)
  }
  stops <- seq(0, 1, length.out = length(scale$colors))
  colorscale <- Map(function(s, col) list(s, col), stops, scale$colors)

  plotly::plot_ly(
    x = colnames(mat),
    y = rownames(mat),
    z = mat,
    type = "heatmap",
    colorscale = colorscale,
    zmin = scale$lo, zmax = scale$hi,
    colorbar = list(title = "median\nnMCC"),
    hovertemplate = paste0(
      "<b>Trained on:</b> %{y}<br>",
      "<b>Tested on:</b> %{x}<br>",
      "<b>Median nMCC:</b> %{z:.3f}<extra></extra>"
    )
  ) |>
    plotly::layout(
      title = list(
        text = "Cross-drug performance", x = 0, font = list(size = 13)
      ),
      xaxis = list(title = "Tested on", side = "bottom"),
      yaxis = list(title = "Trained on", autorange = "reversed"),
      margin = list(l = 80, b = 60, t = 50)
    )
}
