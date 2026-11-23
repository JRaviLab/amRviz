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

    # Rescored features are drawn when asked for and the run has them.
    use_rescoring <- shiny::reactive({
      if (!identical(input$network_rescored, "yes")) return(FALSE)
      resc <- core$rescoring()
      !is.null(resc) && any(normalize_species(resc$features[["species"]]) %in%
                             normalize_species(input$network_bug_id))
    })
    dyad_view <- shiny::reactive({
      use_rescoring() && identical(input$network_view, "dyads")
    })

    # Keep disabled controls' selections, enable controls that still apply
    shiny::observe({
      rescored <- use_rescoring()
      features_visible <- !dyad_view() || isTRUE(input$network_dyad_features)
      shinyjs::toggleState("network_include_dyads", condition = !dyad_view())
      shinyjs::toggleClass("network_include_dyads_hint", "text-muted", condition = dyad_view())
      shinyjs::toggleState("network_feature_types", condition = features_visible)
      shinyjs::toggleClass("network_feature_types_panel", "text-muted", condition = !features_visible)
      for (id in c("network_view", "network_resc_quantile", "network_resc_both")) {
        shinyjs::toggleState(id, condition = rescored)
      }
      shinyjs::toggleState("network_dyad_features", condition = dyad_view())
      shinyjs::toggleState("network_resc_dyad", condition = rescored &&
                            (dyad_view() || isTRUE(input$network_include_dyads)))
      has_models <- rescored && any(normalize_species(core$rescoring()$models[["species"]]) %in%
                                     normalize_species(input$network_bug_id))
      shinyjs::toggleState("network_resc_delta", condition = has_models)
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
      if (!dyad_view()) shiny::validate(shiny::need(length(input$network_feature_types) > 0,
                                                   "Please select at least one feature type."))
      if (use_rescoring()) {
        shiny::req(
          input$network_resc_delta, input$network_resc_quantile,
          input$network_resc_dyad
        )
        makeDrugFeatureNetwork(
          NULL,
          input$network_bug_id,
          top_n = input$network_top_n,
          include_dyads = isTRUE(input$network_include_dyads),
          rescoring = core$rescoring(),
          view = if (identical(input$network_view, "dyads")) "dyad" else "feature",
          show_features = isTRUE(input$network_dyad_features),
          feature_types = input$network_feature_types %||% character(),
          weighted_edges = isTRUE(input$network_weighted_edges),
          # One slider per category, and the subtype check; the rest keep everything.
          thresholds = utils::modifyList(rescoringAll(), list(
            delta_mcc = input$network_resc_delta,
            quantile = input$network_resc_quantile,
            both_subtypes = isTRUE(input$network_resc_both),
            dyad_rank_score = input$network_resc_dyad
          ))
        )
      } else {
        makeDrugFeatureNetwork(
          core$topFeatures(),
          input$network_bug_id,
          top_n = input$network_top_n,
          include_dyads = isTRUE(input$network_include_dyads),
          feature_types = input$network_feature_types %||% character(),
          weighted_edges = isTRUE(input$network_weighted_edges)
        )
      }
    })

    invisible(NULL)
  })
}
