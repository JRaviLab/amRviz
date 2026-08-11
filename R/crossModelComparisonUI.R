#' Cross Model Comparison UI Tab
#' @return A tabPanel for cross model comparison visualization
#' @keywords internal
crossModelComparisonUI <- function() {
  shiny::tabPanel(
    title = "Model holdouts",
    icon = shiny::icon("clock"),
    shiny::fluidRow(
      shiny::h3("Model Holdouts Comparison",
        style = "margin-top: 15px; margin-bottom: 15px; font-weight: bold;"
      ),
      style = "padding: 10px;",
      shiny::tagList(
        shiny::column(
          width = 4,
          style = "height: 80px; display: flex; align-items: center;",
          shiny::selectInput(
            inputId = "bug_cross_model_comparison_id",
            label = shiny::tags$label("Bug", style = "font-size: 15px;"),
            choices = character(0),
            multiple = FALSE,
            selectize = TRUE,
            width = "100%"
          )
        ),
        shiny::column(
          width = 4,
          style = "height: 80px; display: flex; align-items: center;",
          shiny::selectInput(
            inputId = "drug_cross_model_comparison_id",
            label = shiny::tags$label("Drug/Drug class", style = "font-size: 15px;"),
            choices = NULL,
            multiple = FALSE,
            selectize = TRUE,
            width = "100%"
          )
        ),
        shiny::column(
          width = 4,
          style = "height: 80px; display: flex; align-items: center;",
          shiny::radioButtons(
            inputId = "cross_model_comparison",
            label = shiny::tags$label("Cross-train models across",
              style = "font-size: 15px;"
            ),
            choices = c(
              "Countries" = "country",
              "Time (5 yr intervals)" = "time"
            ),
            selected = "country",
            inline = TRUE
          )
        )
      )
    ),
    shiny::column(
      width = 12, offset = 0,
      shiny::mainPanel(
        width = 12,
        style = "padding: 0;",
        shiny::tabsetPanel(
          id = "cross_model_comparison_tabset",
          type = "tabs",
          # Tab 1: Accuracy distributions
          shiny::tabPanel(
            "Accuracy distributions",
            shiny::fluidRow(
              style = "padding: 10px;",
              shiny::column(
                width = 6,
                plotly::plotlyOutput("cross_model_ridge_country",
                  height = "400px"
                ),
                style = "padding: 5px; border: 1px solid lightgray;"
              ),
              shiny::column(
                width = 6,
                plotly::plotlyOutput("cross_model_ridge_time",
                  height = "400px"
                ),
                style = "padding: 5px; border: 1px solid lightgray;"
              )
            )
          ),
          # Tab 2: Model performance heatmaps
          shiny::tabPanel(
            "Model performance",
            shiny::fluidRow(
              style = "padding: 10px;",
              shiny::column(
                width = 6,
                plotly::plotlyOutput("cross_model_perf_country",
                  height = "400px"
                ),
                style = "padding: 5px; border: 1px solid lightgray;"
              ),
              shiny::column(
                width = 6,
                plotly::plotlyOutput("cross_model_perf_time",
                  height = "400px"
                ),
                style = "padding: 5px; border: 1px solid lightgray;"
              )
            )
          ),
          # Tab 3: Top features
          shiny::tabPanel(
            "Top features",
            shiny::tagList(
              shiny::fluidRow(
                shiny::column(
                  width = 4,
                  shiny::div(
                    style = "padding: 10px;",
                    shiny::sliderInput(
                      inputId = "cross_model_top_n_features",
                      label = "Top features",
                      min = 0, max = 100, value = 10
                    )
                  )
                )
              ),
              shiny::fluidRow(
                shiny::column(
                  width = 6,
                  plotly::plotlyOutput(
                    "cross_model_feature_importance_plot",
                    height = "100%"
                  ),
                  style = paste0(
                    "padding: 0; height: 400px;",
                    " border: 1px solid lightgray;"
                  )
                ),
                shiny::column(
                  width = 6,
                  DT::dataTableOutput(
                    "cross_model_feature_importance_table"
                  ),
                  style = paste0(
                    "padding: 0; min-height: 400px;",
                    " margin-bottom: 20px;",
                    " border: 1px solid lightgray;"
                  )
                )
              )
            )
          )
        )
      )
    )
  )
}
