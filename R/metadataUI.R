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
              style = "text-align: center; font-family: 'Arial', sans-serif; font-size: 14px; margin-top: 8px; margin-bottom: 8px;",
              "Distribution of AMR phenotypes"
            ),
            .exportBtn("resistance_vs_susceptible_plot"),
            styledBox("resistance_vs_susceptible_ui")
          )
        ),
        column(
          width = 6,
          div(
            class = "plot-container",
            div(
              class = "plot-header",
              style = "text-align: center; font-family: 'Arial', sans-serif; font-size: 14px; margin-top: 8px; margin-bottom: 8px;",
              "Global distribution of resistant phenotypes"
            ),
            .exportBtn("geo_isolate_plot"),
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
              style = "text-align: center; font-family: 'Arial', sans-serif; font-size: 14px; margin-top: 8px; margin-bottom: 8px;",
              "Distribution of AMR phenotypes by year"
            ),
            .exportBtn("r_s_across_time_plot"),
            styledBox("r_s_across_time_ui")
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
                .exportBtn("isolation_source_plot"),
                styledBox("isolation_source_ui")
              ),
              tabPanel(
                "Hosts",
                .exportBtn("host_isolate_plot"),
                styledBox("host_isolate_plot_ui")
              )
            )
          )
        )
      ),
      fluidRow(
        column(
          width = 12,
          div(
            class = "plot-container",
            div(
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
            div(
              style = "padding: 0 10px;",
              selectInput(
                inputId = "metadata_sankey_classes",
                label = tags$label(
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
            .exportBtn("metadata_sankey"),
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
