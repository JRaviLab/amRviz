# Model Performance tab: cascading selectors, boxplots, and nMCC overview heatmaps.

serverModelPerf <- function(input, output, session, core, results_root) {
  # Refresh the drug-class dropdown to classes present for the selected bug.
  shiny::observeEvent(input$bug_ml_perf_id, {
    data <- core$queryData() |>
      dplyr::filter(
        normalize_species(.data$species) %in%
          normalize_species(input$bug_ml_perf_id)
      ) |>
      dplyr::filter(is.na(.data$strat_label) | !nzchar(.data$strat_label))
    drug_class_vec <- data |>
      dplyr::filter(.data$drug_label == "drug_class") |>
      dplyr::pull(.data$drug_or_class) |>
      unique() |>
      sort()

    shiny::updateSelectInput(
      session,
      inputId  = "drug_class_ml_perf_id",
      choices  = c("all", drug_class_vec),
      selected = "all"
    )
  })

  # Refresh the drug dropdown; restrict via metadata class -> drug mapping unless "all".
  shiny::observeEvent(input$drug_class_ml_perf_id, {
    shiny::req(input$drug_class_ml_perf_id)
    base_data <- core$queryData() |>
      dplyr::filter(
        normalize_species(.data$species) %in%
          normalize_species(input$bug_ml_perf_id)
      ) |>
      dplyr::filter(is.na(.data$strat_label) | !nzchar(.data$strat_label))
    drug_vec <- base_data |>
      dplyr::filter(.data$drug_label == "drug") |>
      dplyr::pull(.data$drug_or_class) |>
      unique() |>
      sort()

    if (input$drug_class_ml_perf_id != "all") {
      sp_codes <- normalize_species(input$bug_ml_perf_id)
      sp_dirs <- core$queryData() |>
        dplyr::filter(normalize_species(.data$species) %in% sp_codes) |>
        dplyr::pull(.data$species_label) |>
        unique()
      meta <- dplyr::bind_rows(lapply(sp_dirs, function(sp_dir) {
        fp <- get_metadata_path(sp_dir, results_root)
        if (!is.null(fp)) .read_parquet_safe(fp, verbose = FALSE) else tibble::tibble()
      }))
      if (nrow(meta) && all(c("class_abbr", "drug_abbr") %in% names(meta))) {
        drugs_in_class <- meta |>
          dplyr::filter(.data$class_abbr %in% input$drug_class_ml_perf_id) |>
          dplyr::pull(.data$drug_abbr) |>
          unique()
        drug_vec <- intersect(drug_vec, drugs_in_class)
      }
      shiny::updateSelectInput(
        session,
        inputId  = "drug_ml_perf_id",
        choices  = drug_vec,
        selected = if (length(drug_vec)) drug_vec[1] else NULL
      )
    } else {
      shiny::updateSelectInput(
        session,
        inputId = "drug_ml_perf_id",
        choices = drug_vec,
        selected = if ("GEN" %in% drug_vec) {
          "GEN"
        } else if (length(drug_vec)) {
          drug_vec[1]
        } else {
          NULL
        }
      )
    }
  })

  output$model_perfomance_plot <- plotly::renderPlotly({
    makeModelPerformancePlot(
      core$queryData(),
      input$bug_ml_perf_id,
      input$model_scale,
      input$data_type,
      input$model_metrics,
      input$drug_class_ml_perf_id,
      input$drug_ml_perf_id
    )
  })

  output$mcc_strip_plot <- plotly::renderPlotly({
    makeMCCStripPlot(
      core$queryData(),
      selected_drug_class = input$drug_class_ml_perf_id,
      selected_drug       = input$drug_ml_perf_id
    )
  })

  output$mcc_heatmap <- plotly::renderPlotly({
    makeMCCHeatmap(
      core$queryData(),
      selected_drug_class = input$drug_class_ml_perf_id
    )
  })
}
