# Model holdouts tab: country and time stratified comparisons.

#' Model Holdouts tab server
#'
#' @param id Module id; must match the id passed to `crossModelComparisonUI()`.
#' @param core Shared data from setupServerCore().
#' @return NULL, invisibly.
#' @keywords internal
#' @noRd
crossModelComparisonServer <- function(id, core) {
  shiny::moduleServer(id, function(input, output, session) {
    # Full names for drug, class and country codes, read once per data load.
    names <- shiny::reactive(holdoutNames(core$queryData()))

    # Species choices come from the loaded ML results.
    shiny::observe({
      choices <- core$available_species()
      sel <- if (length(choices)) choices[[1]] else NULL
      shiny::updateSelectizeInput(session, "bug_cross_model_comparison_id",
        choices = choices, selected = sel
      )
    })

    # Refresh drug options to strata actually present for this bug's holdout models.
    shiny::observeEvent(
      input$bug_cross_model_comparison_id,
      {
        menu <- holdoutDrugMenu(
          perf_data = core$queryData(),
          bug = input$bug_cross_model_comparison_id,
          names = names()
        )
        shiny::updateSelectInput(
          session,
          inputId = "drug_cross_model_comparison_id",
          choices = menu,
          selected = if (length(menu)) unlist(menu)[[1]] else NULL
        )
      },
      ignoreInit = FALSE
    )

    output$cross_model_ridge_country <- plotly::renderPlotly({
      shiny::req(input$bug_cross_model_comparison_id)
      makeCrossModelRidgePlot(
        core$queryData(),
        input$bug_cross_model_comparison_id,
        "country",
        names = names()
      )
    })

    output$cross_model_ridge_time <- plotly::renderPlotly({
      shiny::req(input$bug_cross_model_comparison_id)
      makeCrossModelRidgePlot(
        core$queryData(),
        input$bug_cross_model_comparison_id,
        "time",
        names = names()
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
        "country",
        names = names()
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
        "time",
        names = names()
      )
    })

    output$cross_model_feature_importance_plot <- plotly::renderPlotly({
      shiny::req(
        input$bug_cross_model_comparison_id,
        input$drug_cross_model_comparison_id,
        input$cross_model_comparison
      )
      p <- makeCrossModelFeatureImportancePlot(
        core$topFeatures(),
        input$bug_cross_model_comparison_id,
        input$drug_cross_model_comparison_id,
        input$cross_model_comparison,
        input$cross_model_top_n_features,
        names = names()
      )
      # The drug menu also lists drugs with only cross-tested models, which the
      # performance grid shows but which have no top features; say so rather
      # than drawing an empty plot.
      p %||% plotly_placeholder(holdoutFeatureEmptyMessage(
        .holdout_label(input$drug_cross_model_comparison_id, names()$drug),
        input$cross_model_comparison
      ))
    })

    output$cross_model_feature_importance_table <- DT::renderDataTable({
      rows <- holdoutFeatureRows(
        core$topFeatures(),
        input$bug_cross_model_comparison_id,
        input$drug_cross_model_comparison_id,
        input$cross_model_comparison
      )
      if (!nrow(rows)) {
        return(NULL)
      }
      makeFeatureImportTable(rows)
    })

    invisible(NULL)
  })
}
