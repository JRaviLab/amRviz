# Model holdouts tab: country and time stratified comparisons.

serverCrossModel <- function(input, output, session, core) {
  # Refresh drug options to strata actually present for this bug's holdout models.
  shiny::observeEvent(
    input$bug_cross_model_comparison_id,
    {
      drugs_vec <- getHoldoutsDrugChoices(
        perf_data = core$queryData(),
        bug = input$bug_cross_model_comparison_id
      )
      shiny::updateSelectInput(
        session,
        inputId = "drug_cross_model_comparison_id",
        choices = drugs_vec,
        selected = if (length(drugs_vec)) drugs_vec[1] else NULL
      )
    },
    ignoreInit = FALSE
  )

  output$cross_model_ridge_country <- plotly::renderPlotly({
    shiny::req(input$bug_cross_model_comparison_id)
    makeCrossModelRidgePlot(
      core$queryData(),
      input$bug_cross_model_comparison_id,
      "country"
    )
  })

  output$cross_model_ridge_time <- plotly::renderPlotly({
    shiny::req(input$bug_cross_model_comparison_id)
    makeCrossModelRidgePlot(
      core$queryData(),
      input$bug_cross_model_comparison_id,
      "time"
    )
  })

  output$cross_model_perf_country <- plotly::renderPlotly({
    shiny::req(
      input$bug_cross_model_comparison_id,
      input$drug_cross_model_comparison_id
    )
    makeCrossModelPerformancePlot(
      core$queryData(),
      input$bug_cross_model_comparison_id,
      input$drug_cross_model_comparison_id,
      "country"
    )
  })

  output$cross_model_perf_time <- plotly::renderPlotly({
    shiny::req(
      input$bug_cross_model_comparison_id,
      input$drug_cross_model_comparison_id
    )
    makeCrossModelPerformancePlot(
      core$queryData(),
      input$bug_cross_model_comparison_id,
      input$drug_cross_model_comparison_id,
      "time"
    )
  })

  output$cross_model_feature_importance_plot <- plotly::renderPlotly({
    shiny::req(
      input$bug_cross_model_comparison_id,
      input$drug_cross_model_comparison_id,
      input$cross_model_comparison
    )
    makeCrossModelFeatureImportancePlot(
      core$topFeatures(),
      input$bug_cross_model_comparison_id,
      input$drug_cross_model_comparison_id,
      input$cross_model_comparison,
      input$cross_model_top_n_features
    )
  })

  output$cross_model_feature_importance_table <- DT::renderDataTable({
    strat <- if (isTRUE(input$cross_model_comparison == "country")) {
      "country"
    } else {
      "year"
    }
    tf <- core$topFeatures() |>
      dplyr::filter(
        normalize_species(.data$species) %in%
          normalize_species(input$bug_cross_model_comparison_id)
      ) |>
      dplyr::filter(.data$drug_or_class %in% input$drug_cross_model_comparison_id) |>
      dplyr::filter(.data$strat_label == strat) |>
      dplyr::filter(!isTRUE(.data$cross_test))
    if (!nrow(tf)) {
      return(NULL)
    }
    makeFeatureImportTable(tf)
  })
}
