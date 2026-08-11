#' Metadata UI Tab
#' @return A tabPanel for metadata visualization
#' @keywords internal
metadataUI <- function() {
  shiny::tabPanel(
    title = "Metadata",
    icon = shiny::icon("globe"),
    shiny::fluidPage(
      # Main Content Section
      shiny::column(
        width = 12,
        amr_select(
          "bug_metadata_id",
          "Bug",
          choices = character(0),
          multiple = FALSE,
          selected = NULL
        )
      ),
      shiny::fluidRow(
        # Quick Stats Section
        shiny::column(
          width = 12,
          shiny::div(
            class = "quick-stats-container",
            style = "padding: 10px; margin-bottom: 20px; border: 1px solid lightgray; border-radius: 5px;",
            shiny::uiOutput("quick_metadata_stats")
          )
        )
      ),
      shiny::fluidRow(
        shiny::column(
          width = 6,
          shiny::div(
            class = "plot-container",
            shiny::div(
              class = "plot-header",
              style = "text-align: center; font-family: 'Arial', sans-serif; font-size: 14px; margin-top: 8px; margin-bottom: 8px;",
              "Distribution of AMR phenotypes"
            ),
            styledBox("resistance_vs_susceptible_ui")
          )
        ),
        shiny::column(
          width = 6,
          shiny::div(
            class = "plot-container",
            shiny::div(
              class = "plot-header",
              style = "text-align: center; font-family: 'Arial', sans-serif; font-size: 14px; margin-top: 8px; margin-bottom: 8px;",
              "Global distribution of resistant phenotypes"
            ),
            styledBox("geo_isolate_plot_ui")
          )
        )
      ),
      shiny::fluidRow(
        shiny::column(
          width = 6,
          shiny::div(
            class = "plot-container",
            shiny::div(
              class = "plot-header",
              style = "text-align: center; font-family: 'Arial', sans-serif; font-size: 14px; margin-top: 8px; margin-bottom: 8px;",
              "Distribution of AMR phenotypes by year"
            ),
            styledBox("r_s_across_time_ui")
          )
        ),
        shiny::column(
          width = 6,
          shiny::div(
            class = "plot-container",
            shiny::uiOutput("isolation_source_header"),
            tabBox(
              id = "isolation_source_tabset",
              width = 12,
              shiny::tabPanel(
                "Isolation sources",
                styledBox("isolation_source_ui")
              ),
              shiny::tabPanel(
                "Hosts",
                styledBox("host_isolate_plot_ui")
              )
            )
          )
        )
      ),
      shiny::fluidRow(
        shiny::column(
          width = 12,
          shiny::div(
            class = "plot-container",
            shiny::div(
              class = "plot-header",
              style = paste0(
                "text-align: center; font-family: 'Arial', sans-serif;",
                " font-size: 14px; margin-top: 10px; margin-bottom: 8px;"
              ),
              paste(
                "Resistance flow:",
                "phenotype -> drug class -> drug ->",
                "country -> host -> isolation source"
              )
            ),
            shiny::div(
              style = "padding: 0 10px;",
              shiny::selectInput(
                inputId = "metadata_sankey_classes",
                label = shiny::tags$label(
                  "Drug classes (top 3 by default)",
                  style = "font-size: 13px;"
                ),
                choices = NULL,
                selected = NULL,
                multiple = TRUE,
                selectize = TRUE,
                width = "100%"
              )
            ),
            shinycssloaders::withSpinner(
              networkD3::sankeyNetworkOutput(
                "metadata_sankey",
                height = "550px"
              )
            )
          )
        )
      )
    )
  )
}
