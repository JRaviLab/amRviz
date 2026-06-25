# Metadata tab: plots, sankey, and per-bug summary stats.

# Read and row-bind metadata parquets for one or more bug ids.
.read_metadata_for_bug <- function(bug_ids, results_root) {
  purrr::map_dfr(bug_ids, function(x) {
    fp <- get_metadata_path(x, results_root)
    if (!is.null(fp) && file.exists(fp)) {
      arrow::read_parquet(fp) |> dplyr::mutate(species = x)
    } else {
      tibble::tibble()
    }
  })
}

# Coarsen the lab method into Laboratory vs Computational evidence.
.add_evidence_column <- function(df) {
  df |>
    dplyr::mutate(
      genome_drug.evidence = dplyr::case_when(
        genome_drug.laboratory_typing_method %in%
          c("Disk diffusion", "MIC", "Broth dilution", "Agar dilution") ~
          "Laboratory Method",
        genome_drug.laboratory_typing_method ==
          "Computational Prediction" ~ "Computational Method",
        TRUE ~ genome_drug.laboratory_typing_method
      )
    )
}


serverMetadata <- function(input, output, session, results_root) {
  output$geo_isolate_plot_ui <- shiny::renderUI({
    shiny::fluidRow(shiny::div(shinydashboard::box(
      title = "", width = 12,
      plotly::plotlyOutput("geo_isolate_plot")
    )))
  })

  output$r_s_across_time_ui <- shiny::renderUI({
    shiny::fluidRow(shinydashboard::box(
      title = "", width = 12,
      plotly::plotlyOutput("r_s_across_time_plot")
    ))
  })

  output$host_isolate_plot_ui <- shiny::renderUI({
    shiny::fluidRow(shinydashboard::box(
      title = "", width = 12,
      plotly::plotlyOutput("host_isolate_plot")
    ))
  })

  output$isolation_source_ui <- shiny::renderUI({
    shiny::fluidRow(shinydashboard::box(
      title = "", width = 12,
      plotly::plotlyOutput("isolation_source_plot")
    ))
  })

  output$resistance_vs_susceptible_ui <- shiny::renderUI({
    shiny::fluidRow(shinydashboard::box(
      title = "", width = 12,
      plotly::plotlyOutput("resistance_vs_susceptible_plot")
    ))
  })

  output$isolation_source_header <- shiny::renderUI({
    title <- if (!is.null(input$isolation_source_tabset) &&
      input$isolation_source_tabset == "Hosts") {
      "Distribution of genomes across hosts"
    } else {
      "Distribution of genomes across isolation sources"
    }
    shiny::div(
      class = "plot-header",
      style = paste0(
        "text-align: center; font-family: 'Arial', sans-serif;",
        " font-size: 14px; margin-top: 8px; margin-bottom: 8px;"
      ),
      title
    )
  })

  # Single-bug metadata used by the sankey and its drug-class selector.
  metadata_for_bug <- shiny::reactive({
    shiny::req(input$bug_metadata_id)
    fp <- get_metadata_path(input$bug_metadata_id, results_root)
    if (is.null(fp) || !file.exists(fp)) {
      return(NULL)
    }
    arrow::read_parquet(fp) |>
      dplyr::mutate(species = input$bug_metadata_id)
  })

  # Default sankey selection: the 3 most-represented drug classes.
  shiny::observe({
    meta <- metadata_for_bug()
    if (is.null(meta) || !nrow(meta) || !"drug_class" %in% names(meta)) {
      return()
    }
    classes <- sort(unique(meta$drug_class[!is.na(meta$drug_class)]))
    top3 <- meta |>
      dplyr::filter(!is.na(.data$drug_class)) |>
      dplyr::count(.data$drug_class, name = "n") |>
      dplyr::arrange(dplyr::desc(.data$n)) |>
      dplyr::slice_head(n = 3) |>
      dplyr::pull(.data$drug_class)
    shiny::updateSelectInput(
      session, "metadata_sankey_classes",
      choices = classes, selected = top3
    )
  })

  output$metadata_sankey <- networkD3::renderSankeyNetwork({
    meta <- metadata_for_bug()
    if (is.null(meta) || !nrow(meta)) {
      return(NULL)
    }
    makeMetadataSankey(meta, drug_classes = input$metadata_sankey_classes)
  })

  shiny::observeEvent(input$bug_metadata_id, {
    output$quick_metadata_stats <- shiny::renderUI({
      makeQuickStats(.read_metadata_for_bug(input$bug_metadata_id, results_root))
    })
  })

  output$resistance_vs_susceptible_plot <- plotly::renderPlotly({
    shiny::req(input$bug_metadata_id)
    makeDatAvailabilityPlot(
      .read_metadata_for_bug(input$bug_metadata_id, results_root)
    )
  })

  output$geo_isolate_plot <- plotly::renderPlotly({
    data <- .read_metadata_for_bug(input$bug_metadata_id, results_root) |>
      dplyr::filter(genome.isolation_country != "") |>
      .add_evidence_column() |>
      dplyr::group_by(genome.isolation_country, genome_drug.antibiotic) |>
      dplyr::summarize(count = dplyr::n(), .groups = "drop") |>
      dplyr::group_by(genome.isolation_country) |>
      dplyr::summarise(count = sum(count), .groups = "drop") |>
      dplyr::collect()
    makeGeoChloroPlot(data)
  })

  output$r_s_across_time_plot <- plotly::renderPlotly({
    data <- .read_metadata_for_bug(input$bug_metadata_id, results_root) |>
      .add_evidence_column() |>
      dplyr::filter(!is.na(genome.collection_year)) |>
      dplyr::group_by(
        genome_drug.antibiotic,
        genome_drug.resistant_phenotype,
        genome.isolation_country,
        genome.collection_year
      ) |>
      dplyr::summarize(n = dplyr::n(), .groups = "drop") |>
      dplyr::collect()
    makeTimeSeriesAMRPlot(data, input$amr_drug_search)
  })

  output$host_isolate_plot <- plotly::renderPlotly({
    data <- .read_metadata_for_bug(input$bug_metadata_id, results_root) |>
      dplyr::filter(genome.host_common_name != "") |>
      .add_evidence_column() |>
      dplyr::collect()
    makeHostIsolatePlot(data)
  })

  output$isolation_source_plot <- plotly::renderPlotly({
    data <- .read_metadata_for_bug(input$bug_metadata_id, results_root) |>
      dplyr::filter(genome.host_common_name != "") |>
      .add_evidence_column() |>
      dplyr::collect()
    makeIsolationSourcesPlot(data)
  })
}
