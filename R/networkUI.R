#' Network UI Tab
#' @return A tabPanel for drug-feature network visualization
#' @keywords internal
networkUI <- function() {
  shiny::tabPanel(
    title = "Network",
    icon = shiny::icon("circle-nodes"),
    shiny::fluidPage(
      shiny::h3("Drug-Feature Network",
        style = paste0(
          "font-weight: bold; margin-top: 15px; margin-bottom: 15px;"
        )
      ),
      shiny::p(paste(
        "Interactive force-directed graph linking drugs or drug classes",
        "(orange) to their top predictive features (blue).",
        "Drag nodes, zoom, and hover for labels."
      )),
      shiny::fluidRow(
        shiny::column(
          width = 3,
          style = "padding: 10px;",
          shiny::selectInput(
            inputId = "network_bug_id",
            label = shiny::tags$label("Bug", style = "font-size: 15px;"),
            choices = character(0),
            multiple = FALSE,
            selectize = TRUE,
            width = "100%"
          )
        ),
        shiny::column(
          width = 3,
          style = "padding: 10px;",
          shiny::sliderInput(
            inputId = "network_top_n",
            label = shiny::tags$label(
              "Top features per drug",
              style = "font-size: 15px;"
            ),
            min = 1, max = 30, value = 5
          )
        ),
        shiny::column(
          width = 3,
          style = "padding: 10px;",
          shiny::checkboxInput(
            inputId = "network_include_clusters",
            label = "Show clusters",
            value = FALSE
          ),
          shiny::tags$span(
            "Adds cluster (fig) nodes from annotations.",
            style = "font-size: 11px; color: #666;"
          )
        )
        # ,
        # column(
        #   width = 3,
        #   style = "padding: 10px;",
        #   checkboxInput(
        #     inputId = "network_include_cogs",
        #     label = "Show COGs",
        #     value = FALSE
        #   ),
        #   tags$span(
        #     "Adds COG nodes from annotations.",
        #     style = "font-size: 11px; color: #666;"
        #   )
        # )
      ),
      shiny::fluidRow(
        shiny::column(
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
