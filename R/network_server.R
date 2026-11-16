# Drug-feature network tab.

#' Drug-Feature Network tab server
#'
#' @param id Module id; must match the id passed to `networkUI()`.
#' @param core Shared data from setupServerCore().
#' @param results_root User results root, or NULL for the packaged demo data.
#'   Unused: the network reads each model's own ORB.
#' @return NULL, invisibly.
#' @keywords internal
#' @noRd
networkServer <- function(id, core, results_root = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    # Species choices come from the loaded ML results.
    shiny::observe({
      choices <- core$available_species()
      sel <- if (length(choices)) choices[[1]] else NULL
      shiny::updateSelectizeInput(session, "network_bug_id",
        choices = choices, selected = sel
      )
    })

    # The dyad rank score only means something while dyads are shown.
    shiny::observe({
      shinyjs::toggleState(
        "network_resc_dyad",
        condition = isTRUE(input$network_include_dyads)
      )
    })

    # Rescored features are drawn when asked for and the run has them.
    use_rescoring <- shiny::reactive({
      identical(input$network_rescored, "yes") && !is.null(core$rescoring())
    })

    output$network_rescoring_note <- shiny::renderUI({
      if (identical(input$network_rescored, "yes") && is.null(core$rescoring())) {
        shiny::tags$span(
          "The selected runs have no rescoring.",
          style = "font-size: 12px; color: #a33b52;"
        )
      }
    })

    output$drug_feature_network <- networkD3::renderForceNetwork({
      shiny::req(input$network_bug_id, input$network_top_n)
      if (use_rescoring()) {
        shiny::req(
          input$network_resc_mcc, input$network_resc_delta,
          input$network_resc_quantile, input$network_resc_core,
          input$network_resc_cumulative, input$network_resc_dyad
        )
        makeDrugFeatureNetwork(
          NULL,
          input$network_bug_id,
          top_n = input$network_top_n,
          include_dyads = isTRUE(input$network_include_dyads),
          rescoring = core$rescoring(),
          thresholds = list(
            mcc = input$network_resc_mcc,
            delta_mcc = input$network_resc_delta,
            quantile = input$network_resc_quantile,
            core = input$network_resc_core,
            cumulative = input$network_resc_cumulative,
            both_subtypes = isTRUE(input$network_resc_both),
            dyad_rank_score = input$network_resc_dyad
          )
        )
      } else {
        makeDrugFeatureNetwork(
          core$topFeatures(),
          input$network_bug_id,
          top_n = input$network_top_n,
          include_dyads = isTRUE(input$network_include_dyads)
        )
      }
    })

    invisible(NULL)
  })
}
