# Metadata tab server: plots, sankey, and per-bug summary stats.

#' Metadata tab server
#'
#' Fills the species selector from the metadata parquets under `results_root`
#' (so the tab works without ML results), reads the selected species' metadata
#' once, and renders every output on the tab from that one read.
#'
#' @param id Module id; must match the id passed to `metadataUI()`.
#' @param results_root User results root, or NULL for the packaged demo data.
#' @return NULL, invisibly.
#' @keywords internal
#' @noRd
metadataServer <- function(id, results_root = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    choices <- listMetadataSpecies(results_root)
    shiny::updateSelectizeInput(session, "bug_metadata_id",
      choices = choices,
      selected = if (length(choices)) choices[[1]] else NULL
    )

    # One read per species change, shared by every output below.
    meta <- shiny::reactive({
      shiny::req(input$bug_metadata_id)
      readSpeciesMetadata(input$bug_metadata_id, results_root)
    })

    output$geo_isolate_plot_ui <- shiny::renderUI({
      shiny::fluidRow(shiny::div(shinydashboard::box(
        title = "", width = 12,
        plotly::plotlyOutput(ns("geo_isolate_plot"))
      )))
    })

    output$r_s_across_time_ui <- shiny::renderUI({
      shiny::fluidRow(shinydashboard::box(
        title = "", width = 12,
        plotly::plotlyOutput(ns("r_s_across_time_plot"))
      ))
    })

    output$host_isolate_plot_ui <- shiny::renderUI({
      shiny::fluidRow(shinydashboard::box(
        title = "", width = 12,
        plotly::plotlyOutput(ns("host_isolate_plot"))
      ))
    })

    output$isolation_source_ui <- shiny::renderUI({
      shiny::fluidRow(shinydashboard::box(
        title = "", width = 12,
        plotly::plotlyOutput(ns("isolation_source_plot"))
      ))
    })

    output$resistance_vs_susceptible_ui <- shiny::renderUI({
      shiny::fluidRow(shinydashboard::box(
        title = "", width = 12,
        plotly::plotlyOutput(ns("resistance_vs_susceptible_plot"))
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

    # Default sankey selection: the 3 most-represented drug classes.
    shiny::observe({
      m <- meta()
      if (!nrow(m) || !"drug_class" %in% names(m)) {
        return()
      }
      shiny::updateSelectInput(
        session, "metadata_sankey_classes",
        choices = sort(unique(m$drug_class[!is.na(m$drug_class)])),
        selected = topDrugClasses(m)
      )
    })

    output$metadata_sankey <- networkD3::renderSankeyNetwork({
      m <- meta()
      if (!nrow(m)) {
        return(NULL)
      }
      makeMetadataSankey(m, drug_classes = input$metadata_sankey_classes)
    })

    output$quick_metadata_stats <- shiny::renderUI({
      makeQuickStats(meta())
    })

    output$resistance_vs_susceptible_plot <- plotly::renderPlotly({
      makeDataAvailabilityPlot(meta())
    })

    output$geo_isolate_plot <- plotly::renderPlotly({
      makeGeoChloroPlot(metadataByCountry(meta()))
    })

    output$r_s_across_time_plot <- plotly::renderPlotly({
      makeTimeSeriesAMRPlot(metadataByYear(meta()), "all")
    })

    output$host_isolate_plot <- plotly::renderPlotly({
      makeHostIsolatePlot(metadataWithHost(meta()))
    })

    output$isolation_source_plot <- plotly::renderPlotly({
      makeIsolationSourcesPlot(metadataWithHost(meta()))
    })

    invisible(NULL)
  })
}
