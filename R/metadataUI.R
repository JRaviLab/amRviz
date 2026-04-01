#' Metadata UI Tab
#' @return A tabPanel for metadata visualization
#' @keywords internal
metadataUI <- function() {
  tabPanel(
    title = "Metadata",
    icon = icon("globe"),
    fluidPage(
      # Main Content Section
      column(
        width = 12,
        amr_select(
          "bug_metadata_id",
          "Bug",
          choices = character(0),
          multiple = FALSE,
          selected = NULL
        )
      ),
      fluidRow(
        # Quick Stats Section
        column(
          width = 12,
          div(
            class = "quick-stats-container",
            style = "padding: 10px; margin-bottom: 20px; border: 1px solid lightgray; border-radius: 5px;",
            uiOutput("quick_metadata_stats")
          )
        )
      ),
      fluidRow(
        column(
          width = 6,
          div(
            class = "plot-container",
            div(
              class = "plot-header",
              style = "text-align: center; font-family: 'Arial', sans-serif; font-size: 14px; margin-top: 4px; margin-bottom: 4px;",
              "Distribution of AMR phenotypes"
            ),
            plotly::plotlyOutput("resistance_vs_susceptible_plot", height = "400px")
          )
        ),
        column(
          width = 6,
          div(
            class = "plot-container",
            div(
              class = "plot-header",
              style = "text-align: center; font-family: 'Arial', sans-serif; font-size: 14px; margin-top: 4px; margin-bottom: 4px;",
              "Global distribution of resistant phenotypes"
            ),
            plotly::plotlyOutput("geo_isolate_plot", height = "300px")
          )
        )
      ),
      fluidRow(
        column(
          width = 6,
          div(
            class = "plot-container",
            div(
              class = "plot-header",
              style = paste0(
                "text-align: center; font-family: 'Arial', sans-serif;",
                " font-size: 14px; margin-top: 10px; margin-bottom: 20px;"
              ),
              "Distribution of AMR phenotypes by year"
            ),
            plotly::plotlyOutput("r_s_across_time_plot", height = "300px")
          )
        ),
        column(
          width = 6,
          div(
            class = "plot-container",
            uiOutput("isolation_source_header"),
            tabBox(
              id = "isolation_source_tabset",
              width = 12,
              tabPanel(
                "Isolation sources",
                plotly::plotlyOutput("isolation_source_plot", height = "400px")
              ),
              tabPanel(
                "Hosts",
                plotly::plotlyOutput("host_isolate_plot", height = "400px")
              )
            )
          )
        )
      )
    )
  )
}
