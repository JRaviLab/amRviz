# Drug-feature network tab.

serverNetwork <- function(input, output, session, core, results_root) {
  output$drug_feature_network <- networkD3::renderForceNetwork({
    shiny::req(input$network_bug_id, input$network_top_n)
    makeDrugFeatureNetwork(
      core$topFeatures(),
      input$network_bug_id,
      top_n = input$network_top_n,
      include_clusters = isTRUE(input$network_include_clusters),
      include_cogs = isTRUE(input$network_include_cogs),
      results_root = results_root
    )
  })
}
