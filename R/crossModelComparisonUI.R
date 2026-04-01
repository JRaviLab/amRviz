#' Cross Model Comparison UI Tab
#' @return A tabPanel for cross model comparison visualization
#' @keywords internal
crossModelComparisonUI <- function() {
  tabPanel(
    title = "Model holdouts",
    icon = icon("clock"),
    fluidRow(
      h3("Model Holdouts Comparison", style = "margin-top: 15px; margin-bottom: 15px; font-weight: bold;"),
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
            label = tags$label("Cross-train models across", style = "font-size: 15px;"),
            choices = c("Countries" = "country", "Time (5 yr intervals)" = "time"),
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
          tabPanel(
            "Model performance",
            tagList(
              column(
                width = 12,
                plotOutput("cross_model_perf_plot", height = "100%"),
                style = "padding: 0; border: 1px solid lightgray; height: 700px;"
              )
            )
          ),
          tabPanel(
            "Top features",
            tagList(
              fluidRow(
                column(
                  width = 4,
                  div(
                    style = "padding: 10px;",
                    sliderInput(inputId = "cross_model_top_n_features", label = "Top features", min = 0, max = 100, value = 10)
                  )
                )
              ),
              fluidRow(
                column(
                  width = 6,
                  plotOutput("cross_model_feature_importance_plot", height = "100%"),
                  style = "padding: 0; height: 600px; border: 1px solid lightgray;"
                ),
                column(
                  width = 6,
                  DT::dataTableOutput("cross_model_feature_importance_table"),
                  style = "padding: 0; height: 600px; border: 1px solid lightgray;"
                )
              )
            )
          )
        )
      )
    )
  )
}
