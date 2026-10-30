#' Query Data UI Tab
#' @param id Module id; must match the id passed to `queryDataServer()`.
#' @param source_ui UI shown above the tables, e.g. the results-folder selector
#'   from `resultsSourceUI()`, or NULL.
#' @return A tabPanel for querying and downloading data
#' @keywords internal
queryDataUI <- function(id, source_ui = NULL) {
  ns <- shiny::NS(id)
  shiny::tabPanel(
    title = "Query Data",
    icon = shiny::icon("table", lib = "font-awesome"),
    value = "query_datatable",
    shiny::fluidPage(
      shiny::h2("Query Data"),
      shiny::p("The data table provides the ability to download the raw data for model performance metrics and top features. The table can be customized by adding/removing specific columns."),
      source_ui,
      shiny::tabsetPanel(
        id = ns("query_tabset"),

        # Performance Metrics Tab
        shiny::tabPanel(
          title = "Performance Metrics",
          value = "performance_metrics_tab",
          shiny::fluidRow(
            shiny::column(
              width = 4,
              shiny::selectizeInput(
                inputId = ns("query_data_columns"),
                label = "Add/Remove Column(s)",
                choices = NULL, # Dynamically updated in the server
                selected = NULL, # Default: None selected
                multiple = TRUE,
                options = list(
                  placeholder = "Select columns...",
                  maxOptions = 40,
                  plugins = list("remove_button"), # Allow remove buttons
                  dropdownParent = "body"
                )
              )
            ),
            shiny::column(
              width = 4,
              offset = 4,
              align = "right",
              # Download button for the customized table
              shiny::downloadButton(
                outputId = ns("query_data_download"),
                label = "Download Selected Data (CSV)",
                class = "btn btn-success"
              )
            )
          ),
          shiny::fluidRow(
            shiny::column(
              shinycssloaders::withSpinner(DT::dataTableOutput(outputId = ns("queryDataTable"))),
              width = 12,
              style = "overflow-x: auto; border: 1px solid lightgray;" # Add styling
            )
          )
        ),

        # Top Features Tab
        shiny::tabPanel(
          title = "Top Features",
          value = "top_features_tab",
          shiny::fluidRow(
            shiny::column(
              width = 4,
              shiny::selectizeInput(
                inputId = ns("top_features_columns"),
                label = "Add/Remove Column(s)",
                choices = NULL, # Dynamically updated in the server
                selected = NULL, # Default: None selected
                multiple = TRUE,
                options = list(
                  placeholder = "Select columns...",
                  maxOptions = 40,
                  plugins = list("remove_button"), # Allow remove buttons
                  dropdownParent = "body"
                )
              )
            ),
            shiny::column(
              width = 4,
              offset = 4,
              align = "right",
              # Download button for Top Features table
              shiny::downloadButton(
                outputId = ns("top_features_download"),
                label = "Download Selected Data (CSV)",
                class = "btn btn-success"
              )
            )
          ),
          shiny::fluidRow(
            shiny::column(
              shinycssloaders::withSpinner(DT::dataTableOutput(outputId = ns("topFeaturesTable"))),
              width = 12,
              style = "overflow-x: auto; border: 1px solid lightgray;" # Add styling
            )
          )
        )
      )
    )
  )
}
