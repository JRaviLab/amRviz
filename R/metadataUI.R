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
              style = "text-align: center; font-family: 'Arial', sans-serif; font-size: 14px; margin-top: 30px; margin-bottom: 30px;",
              "Distribution of AMR phenotypes"
            ),
            styledBox("resistance_vs_susceptible_ui")
          )
        ),
        column(
          width = 6,
          div(
            class = "plot-container",
            div(
              class = "plot-header",
              style = "text-align: center; font-family: 'Arial', sans-serif; font-size: 14px; margin-top: 30px; margin-bottom: 30px;",
              "Global distribution of resistant phenotypes"
            ),
            styledBox("geo_isolate_plot_ui")
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
              style = "text-align: center; font-family: 'Arial', sans-serif; font-size: 14px; margin-top: 10px; margin-bottom: 20px;",
              "Distribution of AMR phenotypes by year"
            ),
            styledBox("r_s_across_time_ui")
          )
        ),
        column(
          width = 6,
          div(
            class = "plot-container",
            div(
              class = "plot-header",
              style = "text-align: center; font-family: 'Arial', sans-serif; font-size: 14px; margin-top: 10px; margin-bottom: 30px;",
              "Distribution of genomes across isolation sources"
            ),
            tabBox(
              width = 12,
              tabPanel(
                "Isolation sources",
                styledBox("isolation_source_ui")
              ),
              tabPanel(
                "Hosts",
                div(
                  class = "plot-container",
                  div(
                    class = "plot-header",
                    style = "text-align: center; font-family: 'Arial', sans-serif; font-size: 14px; margin-top: 10px; margin-bottom: 0px;",
                    "Distribution of genomes across hosts"
                  ),
                  styledBox("host_isolate_plot_ui")
                )
              )
            )
          )
        )
      )
    )
  )
}
