#' Network UI Tab
#' @return A tabPanel for drug-feature network visualization
#' @keywords internal
networkUI <- function() {
  tabPanel(
    title = "Network",
    icon = icon("circle-nodes"),
    fluidPage(
      h3("Drug-Feature Network",
        style = paste0(
          "font-weight: bold; margin-top: 15px; margin-bottom: 15px;"
        )
      ),
      p(paste(
        "Interactive force-directed graph linking drugs or drug classes",
        "(orange) to their top predictive features (blue).",
        "Drag nodes, zoom, and hover for labels."
      )),
      fluidRow(
        column(
          width = 3,
          style = "padding: 10px;",
          selectInput(
            inputId = "network_bug_id",
            label = tags$label("Bug", style = "font-size: 15px;"),
            choices = character(0),
            multiple = FALSE,
            selectize = TRUE,
            width = "100%"
          )
        ),
        column(
          width = 3,
          style = "padding: 10px;",
          sliderInput(
            inputId = "network_top_n",
            label = tags$label(
              "Top features per drug",
              style = "font-size: 15px;"
            ),
            min = 1, max = 30, value = 5
          )
        ),
        column(
          width = 3,
          style = "padding: 10px;",
          checkboxInput(
            inputId = "network_include_clusters",
            label = "Show clusters",
            value = FALSE
          ),
          tags$span(
            "Adds cluster (fig) nodes from annotations.",
            style = "font-size: 11px; color: #666;"
          )
        ),
        column(
          width = 3,
          style = "padding: 10px;",
          checkboxInput(
            inputId = "network_include_cogs",
            label = "Show COGs",
            value = FALSE
          ),
          tags$span(
            "Adds COG nodes from annotations.",
            style = "font-size: 11px; color: #666;"
          )
        )
      ),
      fluidRow(
        column(
          width = 12,
          shinycssloaders::withSpinner(
            networkD3::forceNetworkOutput(
              "drug_feature_network",
              height = "600px"
            )
          ),
          style = "padding: 10px; border: 1px solid lightgray;"
        )
      )
    )
  )
}
