# Drug-feature network tab.

#' Drug-Feature Network tab server
#'
#' @param id Module id; must match the id passed to `networkUI()`.
#' @param core Shared data from setupServerCore().
#' @param results_root User results root, or NULL for the packaged demo data.
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

    output$drug_feature_network <- networkD3::renderForceNetwork({
      shiny::req(input$network_bug_id, input$network_top_n)
      makeDrugFeatureNetwork(
        core$topFeatures(),
        input$network_bug_id,
        top_n = input$network_top_n,
        include_clusters = isTRUE(input$network_include_clusters),
        # include_cogs = isTRUE(input$network_include_cogs),
        results_root = results_root
      )
    })

    invisible(NULL)
  })
}
