#' Network UI Tab
#' @param id Module id; must match the id passed to `networkServer()`.
#' @return A tabPanel for drug-feature network visualization
#' @keywords internal
networkUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::tabPanel(
    title = "Network",
    icon = shiny::icon("circle-nodes"),
    shiny::fluidPage(
      shiny::tags$style(shiny::HTML(.network_slider_css())),
      shiny::h3("Drug-Feature Network",
        style = paste0(
          "font-weight: bold; margin-top: 15px; margin-bottom: 15px;"
        )
      ),
      shiny::p(paste(
        "Interactive force-directed graph linking drugs or drug classes",
        "(dark diamonds) to their top predictive features of each type.",
        "Links are amber where a feature points toward resistance and blue",
        "toward susceptibility. Drag nodes, zoom, and hover for labels."
      )),
      shiny::fluidRow(
        shiny::column(
          width = 3,
          style = "padding: 10px;",
          shiny::selectInput(
            inputId = ns("network_bug_id"),
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
            inputId = ns("network_top_n"),
            label = shiny::tags$label(
              "Top features per type, per drug",
              style = "font-size: 15px;"
            ),
            min = 1, max = 30, value = 2
          )
        ),
        shiny::column(
          width = 3,
          style = "padding: 10px;",
          shiny::checkboxInput(
            inputId = ns("network_include_dyads"),
            label = "Show dyads",
            value = FALSE
          ),
          shiny::tags$span(
            paste(
              "A dyad pairs a protein with its gene.",
              "Adds the dyads each feature belongs to, from the ORB."
            ),
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
          width = 3,
          style = "padding: 10px;",
          shiny::radioButtons(
            inputId = ns("network_rescored"),
            label = shiny::tags$label(
              "Rescored features and dyads",
              style = "font-size: 15px;"
            ),
            choices = c("No" = "no", "Yes" = "yes"),
            selected = "no",
            inline = TRUE
          ),
          shiny::uiOutput(ns("network_rescoring_note")),
          shiny::tags$span(
            "Uses amRml's rescoring across seeds, from the run's results/rescoring.",
            style = "font-size: 11px; color: #666;"
          )
        ),
        shiny::conditionalPanel(
          condition = "input.network_rescored == 'yes'",
          ns = ns,
          # One slider per category. Model performance is grey, features keep
          # Shiny's blue and dyads are green, the colour of dyad nodes. The
          # other rescoring thresholds keep everything (see rescoringAll()).
          shiny::column(
            width = 3,
            class = "network-sliders-grey",
            style = "padding: 10px;",
            shiny::tags$h5("Model performance", style = "font-weight: bold;"),
            shiny::sliderInput(
              inputId = ns("network_resc_delta"),
              label = shiny::tags$label(
                "MCC above shuffled",
                style = "font-size: 15px;"
              ),
              min = 0, max = 1.2, value = rescoringAll()$delta_mcc, step = 0.05
            ),
            shiny::tags$span(
              "Model fits must beat their shuffled-label fit by at least this MCC.",
              style = "font-size: 11px; color: #666;"
            )
          ),
          shiny::column(
            width = 3,
            style = "padding: 10px;",
            shiny::tags$h5("Features", style = "font-weight: bold;"),
            shiny::sliderInput(
              inputId = ns("network_resc_quantile"),
              label = shiny::tags$label(
                "Rank-score quantile",
                style = "font-size: 15px;"
              ),
              min = 0, max = 1, value = rescoringAll()$quantile, step = 0.01
            ),
            shiny::tags$span(
              paste(
                "Rank score is a feature's importance rank within a model, from",
                "0 (least) to 1 (most), as the median across seeds. A feature",
                "must reach this quantile of its drug's rank scores: 0.75 keeps",
                "the top quarter."
              ),
              style = "font-size: 11px; color: #666;"
            ),
            shiny::checkboxInput(
              inputId = ns("network_resc_both"),
              label = "Found in both subtypes",
              value = rescoringAll()$both_subtypes
            ),
            shiny::tags$span(
              "Keep a feature only if it passes in binary and counts.",
              style = "font-size: 11px; color: #666;"
            )
          ),
          shiny::column(
            width = 3,
            class = "network-sliders-green",
            style = "padding: 10px;",
            shiny::tags$h5("Dyads", style = "font-weight: bold;"),
            shiny::sliderInput(
              inputId = ns("network_resc_dyad"),
              label = shiny::tags$label(
                "Dyad rank score",
                style = "font-size: 15px;"
              ),
              min = 0, max = 1, value = rescoringAll()$dyad_rank_score,
              step = 0.01
            ),
            shiny::tags$span(
              "Only applies when Show dyads is on.",
              style = "font-size: 11px; color: #666;"
            )
          )
        )
      ),
      shiny::fluidRow(
        shiny::column(
          width = 12,
          shinycssloaders::withSpinner(
            networkD3::forceNetworkOutput(
              ns("drug_feature_network"),
              height = "600px"
            )
          ),
          style = "padding: 10px; border: 1px solid lightgray;"
        )
      )
    )
  )
}


#' CSS that colours the grey and green groups of sliders
#'
#' Shiny draws sliders with ionRangeSlider; the skin's bar, handle and value
#' bubbles are recoloured inside each group.
#'
#' @return A CSS string.
#' @keywords internal
#' @noRd
.network_slider_css <- function() {
  one <- function(group, bar, dark) {
    sprintf(paste0(
      ".%1$s .irs--shiny .irs-bar { background: %2$s; border-color: %2$s; }",
      ".%1$s .irs--shiny .irs-handle { border-color: %3$s; }",
      ".%1$s .irs--shiny .irs-handle:hover, .%1$s .irs--shiny .irs-handle.state_hover",
      " { background-color: #f2f2f2; }",
      ".%1$s .irs--shiny .irs-from, .%1$s .irs--shiny .irs-to,",
      " .%1$s .irs--shiny .irs-single { background-color: %3$s; }",
      ".%1$s .irs--shiny .irs-from:before, .%1$s .irs--shiny .irs-to:before,",
      " .%1$s .irs--shiny .irs-single:before { border-top-color: %3$s; }"
    ), group, bar, dark)
  }
  paste(
    one("network-sliders-grey", "#a8a8a8", "#7a7a7a"),
    one("network-sliders-green", "#8ccf4a", .NETWORK_DYAD_COLOR)
  )
}
