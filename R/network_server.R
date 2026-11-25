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
    available_feature_types <- shiny::reactive({
      rescored <- use_rescoring()
      rows <- if (rescored) core$rescoring()$features else core$topFeatures()
      if (is.null(rows) || !nrow(rows)) return(character())
      rows <- rows[normalize_species(rows$species) %in% normalize_species(input$network_bug_id), ]
      if (!rescored) {
        rows <- rows[is.na(rows$strat_label) | !nzchar(rows$strat_label), ]
        rows <- rows[!rows$cross_test %in% TRUE, ]
      }
      unique(rows$feature_type)
    })
    # A numeric field can be empty while someone replaces its value.
    # Keep the last valid limit until a positive whole number is entered.
    top_dyads <- shiny::reactiveVal(10)
    shiny::observeEvent(input$network_top_dyads, {
      value <- input$network_top_dyads
      if (is.numeric(value) && length(value) == 1L && is.finite(value) &&
          value >= 1 && value == floor(value)) top_dyads(value)
    }, ignoreNULL = FALSE)

    shiny::observeEvent(input$network_reset, {
      defaults <- rescoringAll()
      top_dyads(10)
      shiny::updateSliderInput(session, "network_top_n", value = 2)
      shiny::updateNumericInput(session, "network_top_dyads", value = 10)
      shiny::updateCheckboxGroupInput(session, "network_feature_types", selected = names(SCALE_COLORS))
      shiny::updateCheckboxInput(session, "network_include_dyads", value = FALSE)
      shiny::updateCheckboxInput(session, "network_weighted_edges", value = TRUE)
      shiny::updateCheckboxInput(session, "network_dyad_features", value = FALSE)
      shiny::updateRadioButtons(session, "network_rescored", selected = "no")
      shiny::updateRadioButtons(session, "network_view", selected = "dyads")
      shiny::updateSliderInput(session, "network_resc_delta", value = defaults$delta_mcc)
      shiny::updateSliderInput(session, "network_resc_quantile", value = defaults$quantile)
      shiny::updateSliderInput(session, "network_resc_dyad", value = defaults$dyad_rank_score)
      shiny::updateCheckboxInput(session, "network_resc_both", value = defaults$both_subtypes)
    })

    network_thresholds <- shiny::reactive(utils::modifyList(rescoringAll(), list(
      delta_mcc = input$network_resc_delta %||% 0,
      quantile = input$network_resc_quantile %||% 0,
      both_subtypes = isTRUE(input$network_resc_both)
    )))
    dyad_range <- shiny::reactive({
      if (!use_rescoring() || (!dyad_view() && !isTRUE(input$network_include_dyads))) return(NULL)
      # Ignore this slider when finding its range. Include other filters and 
      # top-dyad limit so the range describes the graph that people  can actually see
      net <- networkTables(NULL, input$network_bug_id,
        top_n = input$network_top_n %||% 2, include_dyads = TRUE,
        rescoring = core$rescoring(), thresholds = network_thresholds(),
        view = if (dyad_view()) "dyad" else "feature",
        feature_types = input$network_feature_types %||% character(),
        top_dyads = top_dyads())
      attr(net, "dyad_rank_range")
    })
    dyad_cutoff <- shiny::reactive({
      scores <- dyad_range()
      if (length(scores) != 2L || any(!is.finite(scores))) return(0)
      scores[1] + (input$network_resc_dyad %||% 0) * diff(scores)
    })
    output$network_dyad_cutoff <- shiny::renderText({
      if (is.null(dyad_range())) return("No matching dyads.")
      sprintf("Rank score \u2265 %.3f%%", 100 * dyad_cutoff())
    })

    # Keep disabled controls' selections, enable controls that still apply
    shiny::observe({
      shinyjs::toggleClass("network_workspace", "network-controls-left",
        condition = identical(input$network_control_layout, "left"))
      rescored <- use_rescoring()
      features_visible <- !dyad_view() || isTRUE(input$network_dyad_features)
      shinyjs::toggleState("network_include_dyads", condition = !dyad_view())
      shinyjs::toggleClass("network_include_dyads_hint", "text-muted", condition = dyad_view())
      # Apply both restrictions together so changing views cannot re-enable
      # types without model-derived features. Preserve their selections.
      shinyjs::runjs(sprintf("{
        var group = document.getElementById(%s);
        if (group) group.querySelectorAll('input[type=checkbox]').forEach(function(input) {
          var supported = %s.indexOf(input.value) !== -1;
          input.disabled = !%s || !supported;
          var label = input.closest('label');
          label.classList.toggle('text-muted', input.disabled);
          label.title = !supported ? 'No model support for this feature type in the selected runs/species.' :
            (input.disabled ? 'Enable Show connected features to select feature types.' : '');
        });
      }", jsonlite::toJSON(session$ns("network_feature_types"), auto_unbox = TRUE),
        jsonlite::toJSON(as.character(available_feature_types()), auto_unbox = FALSE),
        jsonlite::toJSON(features_visible, auto_unbox = TRUE)))
      shinyjs::toggleClass("network_feature_types_panel", "text-muted", condition = !features_visible)
      for (id in c("network_view", "network_resc_quantile", "network_resc_both")) {
        shinyjs::toggleState(id, condition = rescored)
      }
      shinyjs::toggleState("network_dyad_features", condition = dyad_view())
      shinyjs::toggleState("network_top_dyads", condition = dyad_view())
      scores <- dyad_range()
      shinyjs::toggleState("network_resc_dyad", condition = rescored &&
        length(scores) == 2L && all(is.finite(scores)) && diff(scores) > 0)
      has_models <- rescored && any(normalize_species(core$rescoring()$models[["species"]]) %in%
                                     normalize_species(input$network_bug_id))
      shinyjs::toggleState("network_resc_delta", condition = has_models)
    })

    output$network_rescoring_note <- shiny::renderUI({
      if (identical(input$network_rescored, "yes") && !use_rescoring()) {
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
      net <- if (use_rescoring()) {
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
          top_dyads = top_dyads(),
          # One slider per category, and the subtype check; the rest keep everything.
          thresholds = utils::modifyList(network_thresholds(), list(
            dyad_rank_score = dyad_cutoff()
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
      shiny::validate(shiny::need(!is.null(net), if (dyad_view()) {
        "No dyads match the current filters."
      } else {
        "No features match the selected types and filters."
      }))
      htmlwidgets::onRender(net, .NETWORK_LAYOUT_JS, data = list(
        reset = input$network_reset %||% 0,
        ready = identical(input$network_rescored, "no") &&
          isTRUE(input$network_top_n == 2) &&
          setequal(input$network_feature_types, names(SCALE_COLORS)) &&
          !isTRUE(input$network_include_dyads) && isTRUE(input$network_weighted_edges)
      ))
    })

    invisible(NULL)
  })
}
