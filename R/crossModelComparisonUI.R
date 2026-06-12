#' Cross Model Comparison UI Tab
#' @return A tabPanel for cross model comparison visualization
#' @keywords internal
crossModelComparisonUI <- function() {
  tabPanel(
    title = "Model holdouts",
    icon = icon("clock"),
    fluidRow(
      h3("Model Holdouts Comparison",
        style = "margin-top: 15px; margin-bottom: 15px; font-weight: bold;"
      ),
      style = "padding: 10px;",
      tagList(
        column(
          width = 4,
          style = "height: 80px; display: flex; align-items: center;",
          selectInput(
            inputId = "bug_cross_model_comparison_id",
            label = tags$label("Bug", style = "font-size: 15px;"),
            choices = character(0),
            multiple = FALSE,
            selectize = TRUE,
            width = "100%"
          )
        ),
        column(
          width = 4,
          style = "height: 80px; display: flex; align-items: center;",
          selectInput(
            inputId = "drug_cross_model_comparison_id",
            label = tags$label("Drug/Drug class", style = "font-size: 15px;"),
            choices = NULL,
            multiple = FALSE,
            selectize = TRUE,
            width = "100%"
          )
        ),
        column(
          width = 4,
          style = "height: 80px; display: flex; align-items: center;",
          radioButtons(
            inputId = "cross_model_comparison",
            label = tags$label("Cross-train models across",
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
    column(
      width = 12, offset = 0,
      mainPanel(
        width = 12,
        style = "padding: 0;",
        tabsetPanel(
          id = "cross_model_comparison_tabset",
          type = "tabs",
          # Tab 1: Accuracy distributions
          tabPanel(
            "Accuracy distributions",
            fluidRow(
              style = "padding: 10px;",
              column(
                width = 6,
                .exportBtn("cross_model_ridge_country"),
                plotly::plotlyOutput("cross_model_ridge_country",
                  height = "370px"
                ),
                style = "padding: 5px; border: 1px solid lightgray;"
              ),
              column(
                width = 6,
                .exportBtn("cross_model_ridge_time"),
                plotly::plotlyOutput("cross_model_ridge_time",
                  height = "370px"
                ),
                style = "padding: 5px; border: 1px solid lightgray;"
              )
            )
          ),
          # Tab 2: Model performance heatmaps
          tabPanel(
            "Model performance",
            fluidRow(
              style = "padding: 10px;",
              column(
                width = 6,
                div(
                  style = "text-align: right; padding-bottom: 4px;",
                  downloadButton(
                    "cross_model_perf_country_download", "Export (PDF)",
                    class = "btn btn-export"
                  )
                ),
                plotly::plotlyOutput("cross_model_perf_country",
                  height = "360px"
                ),
                style = "padding: 5px; border: 1px solid lightgray;"
              ),
              column(
                width = 6,
                div(
                  style = "text-align: right; padding-bottom: 4px;",
                  downloadButton(
                    "cross_model_perf_time_download", "Export (PDF)",
                    class = "btn btn-export"
                  )
                ),
                plotly::plotlyOutput("cross_model_perf_time",
                  height = "360px"
                ),
                style = "padding: 5px; border: 1px solid lightgray;"
              )
            )
          ),
          # Tab 3: Top features
          tabPanel(
            "Top features",
            tagList(
              fluidRow(
                column(
                  width = 4,
                  div(
                    style = "padding: 10px;",
                    sliderInput(
                      inputId = "cross_model_top_n_features",
                      label = "Top features",
                      min = 0, max = 100, value = 10
                    )
                  )
                )
              ),
              fluidRow(
                column(
                  width = 6,
                  div(
                    style = "text-align: right; padding: 4px;",
                    downloadButton(
                      "cross_model_fi_download", "Export (PDF)",
                      class = "btn btn-export"
                    )
                  ),
                  plotly::plotlyOutput(
                    "cross_model_feature_importance_plot",
                    height = "92%"
                  ),
                  style = paste0(
                    "padding: 0; height: 400px;",
                    " border: 1px solid lightgray;"
                  )
                ),
                column(
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
          ),
          # Tab 4: Cross-drug (static ComplexHeatmap)
          tabPanel(
            "Cross-drug",
            tagList(
              fluidRow(
                style = "padding: 10px;",
                column(
                  width = 12,
                  align = "right",
                  downloadButton(
                    outputId = "cross_drug_heatmap_download",
                    label = "Export (PDF)",
                    class = "btn btn-export"
                  )
                )
              ),
              fluidRow(
                style = "padding: 10px;",
                column(
                  width = 12,
                  shinycssloaders::withSpinner(
                    plotly::plotlyOutput("cross_drug_heatmap", height = "550px")
                  ),
                  style = "padding: 5px; border: 1px solid lightgray;"
                )
              )
            )
          )
        )
      )
    )
  )
}
