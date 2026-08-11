#' Model Performance UI Tab
#' @return A tabPanel for model performance visualization
#' @keywords internal
modelPerfUI <- function() {
  shiny::tabPanel(
    title = "Model performance",
    value = "model_perf_tab",
    icon = shiny::icon("chart-line"),
    shiny::fluidRow(
      shiny::column(
        width = 4,
        shiny::div(
          shiny::h3("Model Performance", style = "margin-top: 15px; margin-bottom: 15px; font-weight: bold;"),
          amr_select(
            "bug_perf_id", "Bug",
            character(0),
            multiple = TRUE,
            selected = NULL
          )
        )
      ),
      shiny::column(
        width = 4,
        style = "padding-top: 52px;",
        amr_select("drug_class_perf_id", "Drug class", NULL, FALSE)
      ),
      shiny::column(
        width = 4,
        style = "padding-top: 52px;",
        amr_select("drug_perf_id", "Drug", NULL, FALSE)
      )
    ),
    shiny::column(
      width = 12,
      shiny::mainPanel(
        width = 12,
        style = "padding: 0;",
        shiny::tabsetPanel(
          shiny::tabPanel(
            "Performance overview",
            shiny::fluidPage(
              shiny::tags$p(
                style = "color: #555; font-size: 10px; padding-top: 10px;",
                "Only baseline models from all the available species. ",
                "Left: MCC distribution per species and molecular scale. ",
                "Right: median MCC by drug class across species, ",
                "molecular scale, and data encoding."
              ),
              shiny::fluidRow(
                shiny::column(4, plotly::plotlyOutput(
                  "mcc_strip_plot",
                  height = "420px"
                )),
                shiny::column(8, plotly::plotlyOutput(
                  "mcc_heatmap",
                  height = "420px"
                ))
              )
            )
          ),
          shiny::tabPanel(
            "Model performance",
             shiny::fluidPage(
              shiny::tags$p(
                style = "color: #555; font-size: 10px; padding-top: 10px;","The filter options above allow to select a subset of species, drug classes, and drugs. ",
                "The filter options below allow to select the model scale and data type")),
            shiny::tagList(
              shiny::fluidRow(
                shiny::column(
                  width = 4,
                  shiny::selectInput(
                    "model_scale",
                    label = shiny::tags$label("Model scale", style = "font-size: 15px;"),
                    choices = c("genes", "domains", "proteins", "cogs", "args"),
                    multiple = TRUE,
                    selectize = TRUE,
                    selected = c("genes", "proteins", "domains", "cogs", "args"),
                    width = "80%"
                  )
                ),
                shiny::column(
                  width = 4,
                  shiny::selectInput(
                    "data_type",
                    label = shiny::tags$label("Data type", style = "font-size: 15px;"),
                    choices = c("count" = "counts", "binary" = "binary"),
                    multiple = FALSE,
                    selectize = TRUE,
                    width = "80%",
                    selected = c("binary")
                  )
                ),
                shiny::column(
                  width = 4,
                  shiny::div(
                    style = "display:none; padding: 10px;",
                    shiny::selectInput(
                      "model_metrics",
                      label = shiny::tags$label("Performance metric", style = "font-size: 15px;"),
                      choices = c("Matthews Correlation Coefficient" = "mcc"),
                      multiple = FALSE,
                      selectize = TRUE,
                      width = "100%"
                    )
                  )
                )
              ),
              shiny::column(
                width = 12,
                shiny::div(
                  style = "height: 600px;",
                  plotly::plotlyOutput("model_perfomance_plot", height = "100%")
                )
              )
            )
          )
        )
      )
    )
  )
}
