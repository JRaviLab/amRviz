# Cross-model holdout visualizations.
#
# Models trained on one country or time period and tested on the same one
# ("Same") or on another ("Different"). Each plot draws one table, built by
# holdoutAccuracyTable(), holdoutGridTable() or holdoutFeatureTable().


# Metric the holdout plots show, and its axis title. Balanced accuracy is 0.5
# for a model no better than chance.
.HOLDOUT_METRIC <- "bal_acc"
.HOLDOUT_METRIC_TITLE <- c(bal_acc = "Balanced accuracy", mcc = "MCC")


#' Strat label for the holdout grouping
#' @param cross_model "country" or "time".
#' @keywords internal
#' @noRd
.holdout_strat <- function(cross_model) {
  if (identical(cross_model, "country")) "country" else "year"
}


# ---- Tables ----

#' Holdout model scores per drug class, by whether the test group matched
#'
#' @param perf_data Performance tibble from loadMLResults().
#' @param bug Species code(s) to include.
#' @param cross_model "country" or "time".
#' @param metric Performance column to report as `value`.
#' @return A tibble with one row per model: `drug_class`, `test` ("Same" or
#'   "Different") and `value`. Drug-class models only; models without the
#'   metric are left out.
#' @keywords internal
#' @noRd
holdoutAccuracyTable <- function(perf_data, bug, cross_model,
                                 metric = .HOLDOUT_METRIC) {
  empty <- tibble::tibble(
    drug_class = character(0), test = character(0), value = numeric(0)
  )
  if (is.null(perf_data) || !is.data.frame(perf_data) || !nrow(perf_data) ||
    !metric %in% names(perf_data)) {
    return(empty)
  }
  df <- perf_data |>
    dplyr::filter(
      normalize_species(.data$species) %in% normalize_species(bug),
      .data$strat_label %in% .holdout_strat(cross_model),
      .data$drug_label == "drug_class",
      !is.na(.data[[metric]])
    )
  if (!nrow(df)) {
    return(empty)
  }
  tibble::tibble(
    drug_class = df$drug_or_class,
    test = ifelse(df$cross_test %in% TRUE, "Different", "Same"),
    value = df[[metric]]
  )
}


#' Mean score for each pair of training and test groups
#'
#' Models tested on the group they were trained on count as tested on that
#' group. Scores are averaged over feature types, encodings and seeds.
#'
#' @inheritParams holdoutAccuracyTable
#' @param drug Drug or drug class to include.
#' @return A tibble with `train`, `test`, `value` (the mean) and `n_models`.
#' @keywords internal
#' @noRd
holdoutGridTable <- function(perf_data, bug, drug, cross_model,
                             metric = .HOLDOUT_METRIC) {
  empty <- tibble::tibble(
    train = character(0), test = character(0),
    value = numeric(0), n_models = integer(0)
  )
  if (is.null(perf_data) || !is.data.frame(perf_data) || !nrow(perf_data) ||
    !metric %in% names(perf_data)) {
    return(empty)
  }
  # `.env$drug`: perf_data may have a column named `drug`, which would
  # otherwise shadow the argument inside the data mask.
  df <- perf_data |>
    dplyr::filter(
      normalize_species(.data$species) %in% normalize_species(bug),
      .data$drug_or_class %in% .env$drug,
      .data$strat_label %in% .holdout_strat(cross_model),
      !is.na(.data[[metric]])
    )
  if (!nrow(df)) {
    return(empty)
  }
  # coalesce(as.character(...)) so an all-NA logical strat_value_test column
  # (no cross-tested rows) does not raise a type error.
  tibble::tibble(
    train = as.character(df$strat_value),
    test = dplyr::coalesce(as.character(df$strat_value_test), as.character(df$strat_value)),
    value = df[[metric]]
  ) |>
    dplyr::group_by(.data$train, .data$test) |>
    dplyr::summarise(
      n_models = dplyr::n(),
      value = mean(.data$value),
      .groups = "drop"
    ) |>
    dplyr::select("train", "test", "value", "n_models")
}


#' Most important features of each group's holdout models
#'
#' Stratified (not cross-tested) models trained on each group, with each
#' feature's importance averaged over seeds, then the `top_n` features per
#' group, scaled 0 to 1 within the group.
#'
#' @param top_data Top-features tibble from loadTopFeat().
#' @param bug Species code(s) to include.
#' @param drug Drug or drug class to include.
#' @param cross_model "country" or "time".
#' @param top_n Number of features per group, or "all".
#' @return A tibble with `feature`, `feature_type`, `group` and `importance`.
#' @keywords internal
#' @noRd
holdoutFeatureTable <- function(top_data, bug, drug, cross_model, top_n = 10) {
  empty <- tibble::tibble(
    feature = character(0), feature_type = character(0),
    group = character(0), importance = numeric(0)
  )
  if (is.null(top_data) || !is.data.frame(top_data) || !nrow(top_data)) {
    return(empty)
  }
  df <- top_data |>
    dplyr::filter(
      normalize_species(.data$species) %in% normalize_species(bug),
      .data$drug_or_class %in% .env$drug,
      .data$strat_label %in% .holdout_strat(cross_model),
      !.data$cross_test %in% TRUE
    )
  if (!nrow(df)) {
    return(empty)
  }
  tbl <- df |>
    dplyr::group_by(
      group = as.character(.data$strat_value),
      feature = .data$Variable, .data$feature_type
    ) |>
    dplyr::summarise(importance = mean(.data$Importance), .groups = "drop")
  if (!identical(top_n, "all")) {
    tbl <- tbl |>
      dplyr::group_by(.data$group) |>
      dplyr::slice_max(.data$importance, n = as.numeric(top_n), with_ties = FALSE) |>
      dplyr::ungroup()
  }
  tbl |>
    dplyr::group_by(.data$group) |>
    dplyr::mutate(importance = {
      rng <- range(.data$importance)
      if (diff(rng) == 0) rep(1, dplyr::n()) else (.data$importance - rng[1]) / diff(rng)
    }) |>
    dplyr::ungroup() |>
    dplyr::select("feature", "feature_type", "group", "importance")
}


# ---- Plots ----

#' Holdout-model feature-importance heatmap
#'
#' Draws holdoutFeatureTable(): the top features of each group's holdout
#' models, scaled within each group. Gray cells are features not among that
#' group's top features.
#'
#' @param top_data Top-features tibble from loadTopFeat().
#' @param bug Species code(s) to include.
#' @param drug Drug/class identifier(s) to include.
#' @param cross_model "country" or "time".
#' @param top_n_features Number of features per stratum, or "all".
#' @param annotated_dir Currently unused; kept for call-site consistency.
#' @return A plotly heatmap, or NULL when there is nothing to plot.
#' @keywords internal
#' @noRd
makeCrossModelFeatureImportancePlot <- function(
  top_data, bug, drug, cross_model, top_n_features,
  annotated_dir = NULL
) {
  tbl <- holdoutFeatureTable(top_data, bug, drug, cross_model, top_n_features)
  if (!nrow(tbl)) {
    return(NULL)
  }
  # Feature rows in order of first appearance by group, so each group's top
  # features sit together.
  tbl <- tbl[order(tbl$group, -tbl$importance), ]
  row_id <- paste0(tbl$feature, " (", tbl$feature_type, ")")
  rows <- unique(row_id)
  groups <- sort(unique(tbl$group))
  mat <- matrix(NA_real_, length(rows), length(groups), dimnames = list(rows, groups))
  mat[cbind(match(row_id, rows), match(tbl$group, groups))] <- tbl$importance

  plotly::plot_ly(
    x = groups,
    y = rows,
    z = mat,
    type = "heatmap",
    colorscale = list(c(0, "#c8e8e8"), c(1, "#4e9a9a")),
    colorbar = list(title = "Importance\n(scaled)"),
    hovertemplate = paste0(
      "<b>Feature:</b> %{y}<br>",
      "<b>Trained on:</b> %{x}<br>",
      "<b>Importance:</b> %{z:.3f}<extra></extra>"
    )
  ) |>
    plotly::layout(
      xaxis = list(title = "", tickangle = -45, side = "top"),
      yaxis = list(title = "", autorange = "reversed"),
      plot_bgcolor = "#eeeeee",
      margin = list(l = 260, b = 20, t = 80)
    )
}


#' Holdout model scores by drug class
#'
#' Draws holdoutAccuracyTable(): one box per drug class and test type, with
#' every model as a point. The dashed line marks chance.
#'
#' @param perf_data Performance tibble from loadMLResults().
#' @param bug Species code(s) to include.
#' @param cross_model "country" or "time".
#' @return A horizontal plotly box plot (empty placeholder when there is no
#'   matching data).
#' @keywords internal
#' @noRd
makeCrossModelRidgePlot <- function(perf_data, bug, cross_model) {
  if (is.null(perf_data) || !is.data.frame(perf_data) || !nrow(perf_data)) {
    return(plotly_placeholder("No data available"))
  }
  df <- holdoutAccuracyTable(perf_data, bug, cross_model)
  if (!nrow(df)) {
    return(plotly_placeholder("No data for selection"))
  }

  strat_label <- if (cross_model == "country") "Countries" else "Time periods"
  metric_title <- .HOLDOUT_METRIC_TITLE[[.HOLDOUT_METRIC]]
  chance <- if (identical(.HOLDOUT_METRIC, "mcc")) 0 else 0.5
  colors <- c(Same = "#d4872a", Different = "#5b8db8")

  p <- plotly::plot_ly()
  for (tt in c("Same", "Different")) {
    sub <- df[df$test == tt, ]
    if (!nrow(sub)) next
    col <- colors[[tt]]
    p <- p |>
      plotly::add_trace(
        type = "box",
        x = sub$value,
        y = sub$drug_class,
        name = tt,
        orientation = "h",
        boxpoints = "all",
        jitter = 0.4,
        pointpos = 0,
        marker = list(color = col, opacity = 0.7, size = 6),
        line = list(color = col),
        fillcolor = paste0(col, "44"),
        hovertemplate = paste0(
          "<b>Drug class:</b> %{y}<br>",
          "<b>", metric_title, ":</b> %{x:.3f}<br>",
          "<b>Test:</b> ", tt, "<extra></extra>"
        )
      )
  }

  p |> plotly::layout(
    title = list(
      text = paste(metric_title, "by drug class -", strat_label),
      x = 0, font = list(size = 13)
    ),
    xaxis = list(
      title = metric_title,
      range = c(min(chance, df$value) - 0.05, 1.03)
    ),
    yaxis = list(title = ""),
    legend = list(title = list(text = "Test group")),
    shapes = list(list(
      type = "line", x0 = chance, x1 = chance, y0 = 0, y1 = 1, yref = "paper",
      line = list(color = "gray", dash = "dash", width = 1)
    )),
    boxmode = "group",
    margin = list(l = 100, t = 50)
  )
}


#' Cross-stratum score heatmap
#'
#' Draws holdoutGridTable(): the mean score of models trained on each group
#' (columns) and tested on each group (rows), printed in each cell. Gray
#' cells are pairs with no model.
#'
#' @param perf_data Performance tibble from loadMLResults().
#' @param bug Species code(s) to include.
#' @param drug Drug/class identifier(s) to include.
#' @param cross_model "country" or "time".
#' @return A plotly heatmap, or an empty placeholder plot (with an explanatory
#'   title) when there is nothing to compare.
#' @keywords internal
#' @noRd
makeCrossModelPerformancePlot <- function(perf_data, bug, drug, cross_model) {
  if (is.null(perf_data) || !is.data.frame(perf_data) || !nrow(perf_data)) {
    return(plotly_placeholder("No data available"))
  }
  tbl <- holdoutGridTable(perf_data, bug, drug, cross_model)
  if (!nrow(tbl)) {
    return(plotly_placeholder("No data for selection"))
  }

  # A cross-stratum grid needs at least two trained-on strata to compare.
  n_strata <- dplyr::n_distinct(tbl$train)
  if (n_strata < 2) {
    unit <- if (cross_model == "country") "countries" else "time periods"
    return(plotly_placeholder(paste0(
      "Not enough data to compare: needs models trained on 2+ ", unit,
      " (this selection has ", n_strata, ")."
    )))
  }

  train <- sort(unique(tbl$train))
  test <- sort(unique(tbl$test))
  idx <- cbind(match(tbl$test, test), match(tbl$train, train))
  value <- matrix(NA_real_, length(test), length(train), dimnames = list(test, train))
  models <- value
  value[idx] <- tbl$value
  models[idx] <- tbl$n_models
  # Cell labels, blank where there is no model (plotly would print 0.00).
  label <- ifelse(is.na(value), "", sprintf("%.2f", value))
  metric_title <- .HOLDOUT_METRIC_TITLE[[.HOLDOUT_METRIC]]

  plotly::plot_ly(
    x = train,
    y = test,
    z = value,
    customdata = models,
    type = "heatmap",
    colorscale = list(
      c(0, "#f2f0f7"), c(0.25, "#cbc9e2"),
      c(0.5, "#9e9ac8"), c(0.75, "#756bb1"),
      c(1, "#54278f")
    ),
    zmin = if (identical(.HOLDOUT_METRIC, "mcc")) 0 else 0.5,
    zmax = 1.0,
    text = label,
    texttemplate = "%{text}",
    colorbar = list(title = gsub(" ", "\n", metric_title)),
    hovertemplate = paste0(
      "<b>Train data:</b> %{x}<br>",
      "<b>Test data:</b> %{y}<br>",
      "<b>", metric_title, ":</b> %{z:.3f}<br>",
      "<b>Models:</b> %{customdata}<extra></extra>"
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
      plot_bgcolor = "#dddddd",
      margin = list(l = 100, b = 60, t = 50)
    )
}
