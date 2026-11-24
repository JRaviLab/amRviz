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


# ---- Names ----

#' Full names for the drug, drug class and country codes
#'
#' Read from the ORB metadata each row points to (`metadata_path`), which
#' pairs each code with its name: `drug_abbr` with `genome_drug.antibiotic`,
#' `class_abbr` with `drug_class`, and `country_abbr` with
#' `genome.isolation_country`.
#'
#' @param data Performance or top-features tibble with a `metadata_path`
#'   column.
#' @return A list of two named character vectors, code -> name: `drug` (drugs
#'   and drug classes) and `country`. Empty when there is no metadata.
#' @keywords internal
#' @noRd
holdoutNames <- function(data) {
  out <- list(drug = character(0), country = character(0))
  if (!is.data.frame(data) || !"metadata_path" %in% names(data)) {
    return(out)
  }
  paths <- unique(stats::na.omit(as.character(data$metadata_path)))
  paths <- paths[file.exists(paths)]
  if (!length(paths)) {
    return(out)
  }
  cols <- c(
    "genome_drug.antibiotic", "drug_abbr", "drug_class", "class_abbr",
    "genome.isolation_country", "country_abbr"
  )
  meta <- dplyr::bind_rows(lapply(paths, function(p) {
    tryCatch(
      arrow::read_parquet(p, col_select = dplyr::any_of(cols)),
      error = function(e) tibble::tibble()
    )
  }))
  pairs <- function(name, code) {
    if (!all(c(name, code) %in% names(meta))) {
      return(character(0))
    }
    d <- meta[!is.na(meta[[name]]) & !is.na(meta[[code]]), c(name, code)]
    d <- d[!duplicated(d[[code]]), ]
    stats::setNames(gsub("_", " ", d[[name]]), d[[code]])
  }
  list(
    drug = c(
      pairs("genome_drug.antibiotic", "drug_abbr"),
      pairs("drug_class", "class_abbr")
    ),
    country = pairs("genome.isolation_country", "country_abbr")
  )
}


#' Label codes with their full names, as "Australia (AUS)"
#'
#' @param code Codes to label.
#' @param lookup Named character vector, code -> name, from holdoutNames().
#' @return Labels; a code with no known name is returned as is.
#' @keywords internal
#' @noRd
.holdout_label <- function(code, lookup) {
  code <- as.character(code)
  name <- unname(lookup[code])
  ifelse(is.na(name), code, paste0(name, " (", code, ")"))
}


#' Name lookup for a holdout grouping: countries have names, years don't
#' @keywords internal
#' @noRd
.holdout_group_names <- function(names, cross_model) {
  if (identical(cross_model, "country")) names$country else character(0)
}


#' Drug menu for the holdouts tab, grouped into drug classes and drugs
#'
#' @param perf_data Performance tibble from loadMLResults().
#' @param bug Species code to include.
#' @param names Name lookups from holdoutNames().
#' @return A list for `selectInput(choices = )`: "Drug classes" and "Drugs"
#'   groups of codes, each labeled like "aminoglycosides (AMG)". Empty groups
#'   are left out; `character(0)` when there are no holdout models.
#' @keywords internal
#' @noRd
holdoutDrugMenu <- function(perf_data, bug, names = holdoutNames(perf_data)) {
  codes <- getHoldoutsDrugChoices(perf_data, bug)
  if (!length(codes)) {
    return(character(0))
  }
  classes <- perf_data$drug_or_class[perf_data$drug_label %in% "drug_class"]
  is_class <- codes %in% classes
  # as.list() so a group with one entry still shows as a group.
  group <- function(x) as.list(stats::setNames(x, .holdout_label(x, names$drug)))
  menu <- list(
    "Drug classes" = group(codes[is_class]),
    "Drugs" = group(codes[!is_class])
  )
  menu[lengths(menu) > 0]
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
#' Stratified (not cross-tested) models trained on each group. Each feature
#' type (scale) is its own model, so features are ranked within each group and
#' feature type: the `top_n` features of each, each divided by the most
#' important one of that group and type (so that one is 1). Dividing, rather
#' than stretching to 0-1, keeps the least important feature shown above 0:
#' it is still a top feature.
#'
#' Importance is averaged over seeds with its sign, so a feature whose seeds
#' disagree on direction counts for less. amRml's `Sign` "NEG" points to
#' Resistant and "POS" to Susceptible.
#'
#' @param top_data Top-features tibble from loadTopFeat().
#' @param bug Species code(s) to include.
#' @param drug Drug or drug class to include.
#' @param cross_model "country" or "time".
#' @param top_n Number of features per group and feature type, or "all".
#' @return A tibble with `feature`, `feature_type`, `group`, `importance` and
#'   `direction` ("Resistant", "Susceptible", or NA when there is no `Sign`).
#' @keywords internal
#' @noRd
holdoutFeatureTable <- function(top_data, bug, drug, cross_model, top_n = 5) {
  empty <- tibble::tibble(
    feature = character(0), feature_type = character(0),
    group = character(0), importance = numeric(0), direction = character(0)
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
  # +1 toward Resistant, -1 toward Susceptible, NA without a sign.
  sign <- if ("Sign" %in% names(df)) {
    dplyr::case_when(df$Sign == "NEG" ~ 1, df$Sign == "POS" ~ -1)
  } else {
    rep(NA_real_, nrow(df))
  }
  df$signed <- df$Importance * dplyr::coalesce(sign, 1)
  df$has_sign <- !is.na(sign)

  tbl <- df |>
    dplyr::group_by(
      group = as.character(.data$strat_value),
      feature = .data$Variable, .data$feature_type
    ) |>
    dplyr::summarise(
      signed = mean(.data$signed),
      has_sign = all(.data$has_sign),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      importance = abs(.data$signed),
      direction = dplyr::case_when(
        !.data$has_sign ~ NA_character_,
        .data$signed > 0 ~ "Resistant",
        .data$signed < 0 ~ "Susceptible"
      )
    )
  if (!identical(top_n, "all")) {
    tbl <- tbl |>
      dplyr::group_by(.data$group, .data$feature_type) |>
      dplyr::slice_max(.data$importance, n = as.numeric(top_n), with_ties = FALSE) |>
      dplyr::ungroup()
  }
  tbl |>
    dplyr::group_by(.data$group, .data$feature_type) |>
    dplyr::mutate(importance = {
      top <- max(.data$importance)
      if (top > 0) .data$importance / top else rep(1, dplyr::n())
    }) |>
    dplyr::ungroup() |>
    dplyr::select("feature", "feature_type", "group", "importance", "direction")
}


#' Holdout top-feature rows for the table beside the heatmap
#'
#' The rows behind holdoutFeatureTable(), one per model and feature, without
#' the columns that the filter makes constant or that only identify the run
#' (`task_id`, `shuffled`, `cross_test`, `strat_value_test`, `species_label`,
#' `metadata_path`). `Sign` is shown as `direction` ("Resistant" or
#' "Susceptible").
#'
#' @inheritParams holdoutFeatureTable
#' @return A tibble, empty when nothing matches.
#' @keywords internal
#' @noRd
holdoutFeatureRows <- function(top_data, bug, drug, cross_model) {
  if (is.null(top_data) || !is.data.frame(top_data) || !nrow(top_data)) {
    return(tibble::tibble())
  }
  df <- top_data |>
    dplyr::filter(
      normalize_species(.data$species) %in% normalize_species(bug),
      .data$drug_or_class %in% .env$drug,
      .data$strat_label %in% .holdout_strat(cross_model),
      !.data$cross_test %in% TRUE
    )
  if ("Sign" %in% names(df)) {
    df <- df |>
      dplyr::mutate(Sign = dplyr::case_when(
        .data$Sign == "NEG" ~ "Resistant",
        .data$Sign == "POS" ~ "Susceptible"
      )) |>
      dplyr::rename(direction = "Sign")
  }
  df |>
    dplyr::select(-dplyr::any_of(c(
      "task_id", "shuffled", "cross_test", "strat_value_test",
      "species_label", "metadata_path"
    )))
}


#' Why a drug's holdout top features are empty
#'
#' @param drug Drug or drug class shown.
#' @param cross_model "country" or "time".
#' @return A one-line message for the empty plot.
#' @keywords internal
#' @noRd
holdoutFeatureEmptyMessage <- function(drug, cross_model) {
  group <- if (identical(cross_model, "country")) "country" else "time period"
  paste0(
    "No ", drug, " models were trained and tested within a single ", group,
    ", so there are no top features to show"
  )
}


# ---- Plots ----

#' Holdout-model feature-importance heatmap
#'
#' Draws holdoutFeatureTable(): the top features of each group's holdout
#' models for each feature type, relative to the most important feature of
#' that group and type. Color shows direction (amber toward resistant, blue
#' toward susceptible) and depth shows importance. Gray cells are features not
#' among that group's top features.
#'
#' @param top_data Top-features tibble from loadTopFeat().
#' @param bug Species code(s) to include.
#' @param drug Drug/class identifier(s) to include.
#' @param cross_model "country" or "time".
#' @param top_n_features Number of features per stratum and feature type, or
#'   "all".
#' @param annotated_dir Currently unused; kept for call-site consistency.
#' @param names Name lookups from holdoutNames(), for country names in the
#'   hover text.
#' @return A plotly heatmap, or NULL when there is nothing to plot.
#' @keywords internal
#' @noRd
makeCrossModelFeatureImportancePlot <- function(
  top_data, bug, drug, cross_model, top_n_features,
  annotated_dir = NULL, names = holdoutNames(top_data)
) {
  tbl <- holdoutFeatureTable(top_data, bug, drug, cross_model, top_n_features)
  if (!nrow(tbl)) {
    return(NULL)
  }
  # Rows grouped by feature type, then by group, most important first, so each
  # group's top features of each type sit together.
  tbl <- tbl[order(tbl$feature_type, tbl$group, -tbl$importance), ]
  row_id <- paste0(tbl$feature, " (", tbl$feature_type, ")")
  rows <- unique(row_id)
  groups <- sort(unique(tbl$group))
  idx <- cbind(match(row_id, rows), match(tbl$group, groups))

  signed <- any(!is.na(tbl$direction))
  z <- matrix(NA_real_, length(rows), length(groups), dimnames = list(rows, groups))
  z[idx] <- if (signed) {
    ifelse(tbl$direction %in% "Susceptible", -tbl$importance, tbl$importance)
  } else {
    tbl$importance
  }

  group_label <- .holdout_label(tbl$group, .holdout_group_names(names, cross_model))
  hover <- matrix(NA_character_, length(rows), length(groups))
  hover[idx] <- paste0(
    "<b>Feature:</b> ", row_id, "<br>",
    "<b>Trained on:</b> ", group_label, "<br>",
    "<b>Importance:</b> ", sprintf("%.2f", tbl$importance),
    " of the top ", tbl$feature_type, " feature",
    ifelse(is.na(tbl$direction), "",
      paste0("<br><b>Points to:</b> ", tolower(tbl$direction))
    )
  )

  # Fixed color range, so color shows importance relative to the top feature
  # rather than stretching the lowest one shown to the palest color.
  color <- if (signed) {
    list(
      colorscale = list(
        c(0, PHENOTYPE_COLORS[["Susceptible"]]), c(0.5, "#f7f7f7"),
        c(1, PHENOTYPE_COLORS[["Resistant"]])
      ),
      zmin = -1, zmax = 1,
      colorbar = list(
        title = "Importance<br>(relative to top<br>feature of its type)",
        tickvals = c(-1, -0.5, 0, 0.5, 1),
        ticktext = c("1 susceptible", "0.5", "0", "0.5", "1 resistant")
      )
    )
  } else {
    list(
      colorscale = list(c(0, "#c8e8e8"), c(1, "#4e9a9a")),
      zmin = 0, zmax = 1,
      colorbar = list(title = "Importance<br>(relative to top<br>feature of its type)")
    )
  }

  plotly::plot_ly(
    x = groups,
    y = rows,
    z = z,
    hovertext = hover,
    type = "heatmap",
    colorscale = color$colorscale,
    zmin = color$zmin, zmax = color$zmax,
    colorbar = color$colorbar,
    hoverongaps = FALSE,
    hoverinfo = "text"
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
#' @param names Name lookups from holdoutNames(), for drug class names in the
#'   hover text.
#' @return A horizontal plotly box plot (empty placeholder when there is no
#'   matching data).
#' @keywords internal
#' @noRd
makeCrossModelRidgePlot <- function(perf_data, bug, cross_model,
                                    names = holdoutNames(perf_data)) {
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
        text = .holdout_label(sub$drug_class, names$drug),
        name = tt,
        orientation = "h",
        boxpoints = "all",
        jitter = 0.4,
        pointpos = 0,
        marker = list(color = col, opacity = 0.7, size = 6),
        line = list(color = col),
        fillcolor = paste0(col, "44"),
        hovertemplate = paste0(
          "<b>Drug class:</b> %{text}<br>",
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
#' @param names Name lookups from holdoutNames(), for country names in the
#'   hover text.
#' @return A plotly heatmap, or an empty placeholder plot (with an explanatory
#'   title) when there is nothing to compare.
#' @keywords internal
#' @noRd
makeCrossModelPerformancePlot <- function(perf_data, bug, drug, cross_model,
                                          names = holdoutNames(perf_data)) {
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
  value[idx] <- tbl$value
  # Cell labels, blank where there is no model (plotly would print 0.00).
  label <- ifelse(is.na(value), "", sprintf("%.2f", value))
  metric_title <- .HOLDOUT_METRIC_TITLE[[.HOLDOUT_METRIC]]

  lookup <- .holdout_group_names(names, cross_model)
  hover <- matrix(NA_character_, length(test), length(train))
  hover[idx] <- paste0(
    "<b>Train data:</b> ", .holdout_label(tbl$train, lookup), "<br>",
    "<b>Test data:</b> ", .holdout_label(tbl$test, lookup), "<br>",
    "<b>", metric_title, ":</b> ", sprintf("%.3f", tbl$value), "<br>",
    "<b>Models:</b> ", tbl$n_models
  )

  plotly::plot_ly(
    x = train,
    y = test,
    z = value,
    hovertext = hover,
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
    hoverongaps = FALSE,
    hoverinfo = "text"
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
