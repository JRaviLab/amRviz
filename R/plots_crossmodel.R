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

