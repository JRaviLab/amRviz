# Bug/Drug feature comparison tab.

serverFeatureImportance <- function(input, output, session, core,
                                    results_root, amrdata_root) {
  output$drug_toggle_ui <- shiny::renderUI({
    choices <- switch(input$drug_toggle_top,
      "bug"   = c("Drug Class" = "class", "Drug" = "drug"),
      "class" = c("Bug" = "bug", "Drug" = "drug"),
      "drug"  = c("Bug" = "bug", "Drug Class" = "class")
    )
    shiny::radioButtons(
      inputId  = "drug_toggle",
      label    = shiny::tags$label("Select by", style = "font-size: 15px;"),
      choices  = choices,
      selected = names(choices)[1],
      inline   = TRUE
    )
  })

  # Refresh the across-bug selector when the user toggles drug vs drug-class.
  shiny::observeEvent(input$across_bug_id, {
    choices <- core$available_species()
    shiny::updateSelectInput(
      session,
      inputId  = "bug_search_amr_across_bug",
      choices  = choices,
      selected = choices
    )
  })

  # Recompute drug/class options for the across-bug view when any filter changes.
  shiny::observeEvent(
    c(
      input$bug_search_amr_across_bug, input$across_bug_id,
      input$bug_drug_comp_model_scale, input$feature_data_type
    ),
    {
      shiny::req(input$bug_search_amr_across_bug)
      message("Across bug feature comparison: ")

      base_tf <- core$topFeatures() |>
        dplyr::filter(normalize_species(.data$species) %in% core$bug_norm_input()) |>
        dplyr::filter(.data$feature_type %in% input$bug_drug_comp_model_scale) |>
        dplyr::filter(.data$feature_subtype %in% input$feature_data_type) |>
        dplyr::filter(is.na(.data$strat_label) | !nzchar(.data$strat_label))

      if (identical(input$across_bug_id, "drug")) {
        drugs_vec <- base_tf |>
          dplyr::filter(.data$drug_label == "drug") |>
          dplyr::pull(.data$drug_or_class) |>
          unique() |>
          sort()
        sel <- intersect("GEN", drugs_vec)
        if (!length(sel)) sel <- utils::head(drugs_vec, 1)
        shiny::updateSelectInput(
          session, "amr_drug_across_bug",
          choices = drugs_vec, selected = sel
        )
      } else {
        drugs_class_vec <- base_tf |>
          dplyr::filter(.data$drug_label == "drug_class") |>
          dplyr::pull(.data$drug_or_class) |>
          unique() |>
          sort()
        sel <- utils::head(drugs_class_vec, 1)
        shiny::updateSelectInput(
          session, "amr_drug_class_across_bug",
          choices = drugs_class_vec, selected = sel
        )
      }
    }
  )

  # Same cascade for the across-drug view; preserves the previous selection where possible.
  shiny::observeEvent(
    c(
      input$bug_search_amr_across_drug, input$across_drug_id,
      input$bug_drug_comp_model_scale, input$feature_data_type
    ),
    {
      shiny::req(input$bug_search_amr_across_drug)
      message("Across drug feature comparison: ")

      bn <- normalize_species(input$bug_search_amr_across_drug)

      base_tf <- core$topFeatures() |>
        dplyr::filter(normalize_species(.data$species) %in% bn) |>
        dplyr::filter(.data$feature_type %in% input$bug_drug_comp_model_scale) |>
        dplyr::filter(.data$feature_subtype %in% input$feature_data_type) |>
        dplyr::filter(is.na(.data$strat_label) | !nzchar(.data$strat_label))

      if (identical(input$across_drug_id, "drug")) {
        drugs_vec <- base_tf |>
          dplyr::filter(.data$drug_label == "drug") |>
          dplyr::pull(.data$drug_or_class) |>
          unique() |>
          sort()
        prev <- shiny::isolate(input$amr_drug_across_drug)
        sel <- prev[prev %in% drugs_vec]
        if (!length(sel)) {
          pref <- c("OXA", "PEN", "MET")
          sel <- intersect(pref, drugs_vec)
          if (!length(sel)) sel <- utils::head(drugs_vec, min(3, length(drugs_vec)))
        }
        shiny::updateSelectInput(
          session, "amr_drug_across_drug",
          choices = drugs_vec, selected = sel
        )
      } else {
        drugs_class_vec <- base_tf |>
          dplyr::filter(.data$drug_label == "drug_class") |>
          dplyr::pull(.data$drug_or_class) |>
          unique() |>
          sort()
        prev <- shiny::isolate(input$amr_drug_class_across_drug)
        sel <- prev[prev %in% drugs_class_vec]
        if (!length(sel)) {
          pref <- c("CEP", "LIN", "MAC", "PEN")
          sel <- intersect(pref, drugs_class_vec)
          if (!length(sel)) sel <- utils::head(drugs_class_vec, min(4, length(drugs_class_vec)))
        }
        shiny::updateSelectInput(
          session, "amr_drug_class_across_drug",
          choices = drugs_class_vec, selected = sel
        )
      }
    }
  )

  output$across_bug_feature_importance_plot <- plotly::renderPlotly({
    if (is.null(input$across_bug_id)) {
      return(NULL)
    }
    amr_drug <- if (input$across_bug_id == "drug") {
      input$amr_drug_across_bug
    } else {
      input$amr_drug_class_across_bug
    }
    makeFeatureImportancePlot(
      core$topFeatures(),
      input$bug_search_amr_across_bug,
      amr_drug,
      input$bug_drug_comp_model_scale,
      input$feature_data_type,
      input$top_n_features,
      input$feature_importance_tabset,
      amrdata_root = amrdata_root,
      results_root = results_root
    )
  })

  output$across_drug_feature_importance_plot <- plotly::renderPlotly({
    if (is.null(input$across_drug_id)) {
      return(NULL)
    }
    amr_drug <- if (input$across_drug_id == "drug") {
      input$amr_drug_across_drug
    } else {
      input$amr_drug_class_across_drug
    }
    makeFeatureImportancePlot(
      core$topFeatures(),
      input$bug_search_amr_across_drug,
      amr_drug,
      input$bug_drug_comp_model_scale,
      input$feature_data_type,
      input$top_n_features,
      input$feature_importance_tabset,
      amrdata_root = amrdata_root,
      results_root = results_root
    )
  })

  # Annotation-joined top features feeding the table, cluster barplot, and ego network.
  enriched_across_bug <- shiny::reactive({
    amr_drug <- if (!is.null(input$across_bug_id) &&
      input$across_bug_id == "drug_class") {
      input$amr_drug_class_across_bug
    } else {
      input$amr_drug_across_bug
    }
    bug <- input$bug_search_amr_across_bug
    tf <- core$filtered_top_features() |>
      dplyr::filter(normalize_species(.data$species) %in% normalize_species(bug)) |>
      dplyr::filter(.data$drug_or_class %in% amr_drug)
    if (!nrow(tf)) {
      return(NULL)
    }
    dplyr::bind_rows(lapply(unique(tf$species), function(sp) {
      enrich_with_annotations(
        tf[tf$species == sp, ],
        species_code = sp,
        results_root = results_root
      )
    }))
  })

  enriched_across_drug <- shiny::reactive({
    amr_drug <- if (!is.null(input$across_drug_id) &&
      input$across_drug_id == "drug_class") {
      input$amr_drug_class_across_drug
    } else {
      input$amr_drug_across_drug
    }
    bug <- input$bug_search_amr_across_drug
    tf <- core$filtered_top_features() |>
      dplyr::filter(normalize_species(.data$species) %in% normalize_species(bug)) |>
      dplyr::filter(.data$drug_or_class %in% amr_drug)
    if (!nrow(tf)) {
      return(NULL)
    }
    enrich_with_annotations(
      tf,
      species_code = bug, results_root = results_root
    )
  })

  output$across_bug_feature_importance_table <- DT::renderDataTable({
    tf <- enriched_across_bug()
    if (is.null(tf) || !nrow(tf)) {
      return(NULL)
    }
    makeFeatureImportTable(tf)
  })

  output$across_drug_feature_importance_table <- DT::renderDataTable({
    tf <- enriched_across_drug()
    if (is.null(tf) || !nrow(tf)) {
      return(NULL)
    }
    makeFeatureImportTable(tf)
  })

  output$across_bug_cog_barplot <- plotly::renderPlotly({
    makeClusterBarChart(enriched_across_bug())
  })
  output$across_drug_cog_barplot <- plotly::renderPlotly({
    makeClusterBarChart(enriched_across_drug())
  })

  output$across_bug_ego_network <- networkD3::renderForceNetwork({
    tf <- enriched_across_bug()
    sel <- input$across_bug_feature_importance_table_rows_selected
    if (is.null(tf) || !nrow(tf) || is.null(sel) || !length(sel)) {
      return(NULL)
    }
    makeFeatureEgoNetwork(tf, tf$Variable[sel])
  })
  output$across_drug_ego_network <- networkD3::renderForceNetwork({
    tf <- enriched_across_drug()
    sel <- input$across_drug_feature_importance_table_rows_selected
    if (is.null(tf) || !nrow(tf) || is.null(sel) || !length(sel)) {
      return(NULL)
    }
    makeFeatureEgoNetwork(tf, tf$Variable[sel])
  })
}
