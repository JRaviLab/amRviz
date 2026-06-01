#' Model Performance UI Tab
#' @return A tabPanel for model performance visualization
#' @keywords internal
modelPerfUI <- function() {
  tabPanel(
    title = "Model performance",
    value = "model_perf_tab",
    icon = icon("chart-line"),
    fluidRow(
      column(
        width = 4,
        div(
          h3("Model Performance", style = "margin-top: 15px; margin-bottom: 15px; font-weight: bold;"),
          amr_select(
            "bug_ml_perf_id", "Bug",
            character(0),
            multiple = TRUE,
            selected = NULL
          )
        )
      ),
      column(
        width = 4,
        style = "padding-top: 52px;",
        amr_select("drug_class_ml_perf_id", "Drug class", NULL, FALSE)
      ),
      column(
        width = 4,
        style = "padding-top: 52px;",
        amr_select("drug_ml_perf_id", "Drug", NULL, FALSE)
      )
    ),
    column(
      width = 12,
      mainPanel(
        width = 12,
        style = "padding: 0;",
        tabsetPanel(
          tabPanel(
            "Model performance",
            tagList(
              fluidRow(
                column(
                  width = 4,
                  selectInput(
                    "model_scale",
                    label = tags$label("Model scale", style = "font-size: 15px;"),
                    choices = c("genes", "domains", "proteins"),
                    multiple = T,
                    selectize = TRUE,
                    selected = c("genes", "proteins", "domains"),
                    width = "80%"
                  )
                ),
                column(
                  width = 4,
                  selectInput(
                    "data_type",
                    label = tags$label("Data type", style = "font-size: 15px;"),
                    choices = c("count" = "counts", "binary" = "binary"),
                    multiple = F,
                    selectize = TRUE,
                    width = "80%",
                    selected = c("binary")
                  )
                ),
                column(
                  width = 4,
                  div(
                    style = "display:none; padding: 10px;",
                    selectInput(
                      "model_metrics",
                      label = tags$label("Performance metric", style = "font-size: 15px;"),
                      choices = c("Matthews Correlation Coefficient" = "nmcc"),
                      multiple = FALSE,
                      selectize = TRUE,
                      width = "100%"
                    )
                  )
                )
              ),
              column(
                width = 12,
                div(
                  style = "height: 600px;",
                  plotly::plotlyOutput("model_perfomance_plot", height = "100%")
                )
              )
            )
          ),
          tabPanel(
            "Performance overview",
            fluidPage(
              tags$p(
                style = "color: #555; font-size: 10px; padding-top: 10px;",
                "Baseline models only (non-stratified, non-cross-test). ",
                "Left: nMCC distribution per species and molecular scale. ",
                "Right: median nMCC by drug class across species, ",
                "molecular scale, and data encoding."
              ),
              fluidRow(
                column(4, plotly::plotlyOutput(
                  "nmcc_strip_plot",
                  height = "420px"
                )),
                column(8, plotly::plotlyOutput(
                  "nmcc_heatmap",
                  height = "420px"
                ))
              )
            )
          )
        )
      )
    )
  )
}
