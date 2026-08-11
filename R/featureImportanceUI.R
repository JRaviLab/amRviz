#' Feature Importance UI Tab
#' @return A tabPanel for feature importance visualization
#' @keywords internal
featureImportanceUI <- function() {
  shiny::tabPanel(
    title = "Bug/Drug feature comparison",
    icon = shiny::icon("bug"),
    shiny::fluidRow(
      shiny::h3("Feature Importance Comparison", style = "margin-top: 15px; margin-bottom: 15px; font-weight: bold;"),
      style = "padding: 10px;",
      shiny::column(
        width = 4,
        style = "height: 80px; display: flex; align-items: center;",
        shiny::selectInput(
          "bug_drug_comp_model_scale",
          label = shiny::tags$label("Model scale", style = "font-size: 15px;"),
          choices = c("genes", "domains", "proteins", "cogs", "args"),
          multiple = FALSE,
          selectize = TRUE,
          selected = "genes"
        )
      ),
      shiny::column(
        width = 4,
        style = "height: 80px; display: flex; align-items: center;",
        shiny::selectInput(
          "feature_data_type",
          label = shiny::tags$label("Data type", style = "font-size: 15px;"),
          choices = c("count" = "counts", "binary" = "binary"),
          multiple = FALSE,
          selectize = TRUE,
          selected = c("counts", "binary")
        )
      ),
      shiny::column(
        width = 4,
        style = "height: 80px; display: flex; align-items: center;",
        shiny::div(
          shiny::sliderInput(inputId = "top_n_features", label = "Top features (per model)", min = 0, max = 100, value = 10)
        )
      )
    ),
    shiny::column(
      width = 12, offset = 0,
      shiny::mainPanel(
        width = 12,
        style = "padding: 0;",
        shiny::tabsetPanel(
          id = "feature_importance_tabset",
          type = "tabs",
          shiny::tabPanel(
            "Across bugs",
            value = "across_bug",
            shiny::tagList(
              shiny::fluidRow(
                shiny::column(
                  width = 6,
                  style = "height: 50px; display: flex; align-items: center;",
                  shiny::h4("Toggle the button to select drug or drug class")
                ),
                shiny::column(
                  width = 4,
                  style = "height: 50px; display: flex; align-items: center;",
                  shiny::radioButtons(
                    inputId = "across_bug_id",
                    label = "",
                    choices = c("Drug class" = "drug_class", "Drug" = "drug"),
                    selected = "drug",
                    inline = TRUE
                  )
                )
              ),
              shiny::fluidRow(
                # Bug selector is the same in both modes, so define it once
                # outside the conditional panels (avoids a duplicate input id).
                shiny::column(
                  width = 6, style = "padding: 0;",
                  amr_select("bug_search_amr_across_bug", "Bug (multi-select)", character(0), selected = NULL)
                ),
                shiny::conditionalPanel(
                  condition = "input.across_bug_id == 'drug'",
                  shiny::column(
                    width = 6, style = "padding: 0;",
                    amr_select("amr_drug_across_bug", "Drug", NULL, FALSE)
                  )
                ),
                shiny::conditionalPanel(
                  condition = "input.across_bug_id == 'drug_class'",
                  shiny::column(
                    width = 6, style = "padding: 0;",
                    amr_select("amr_drug_class_across_bug", "Drug class", NULL, FALSE)
                  )
                )
              ),
              shiny::column(
                width = 6,
                style = "padding: 0; height: 600px; border: 1px solid lightgray;",
                plotly::plotlyOutput("across_bug_feature_importance_plot", height = "100%", width = "100%")
              ),
              shiny::column(
                width = 6,
                style = "padding: 0; height: 600px; border: 1px solid lightgray;",
                DT::dataTableOutput("across_bug_feature_importance_table", height = "100%")
              ),
              shiny::fluidRow(
                style = "padding: 10px 0;",
                shiny::column(
                  width = 6,
                  style = "padding: 0; border: 1px solid lightgray; height: 350px;",
                  plotly::plotlyOutput("across_bug_cog_barplot", height = "100%")
                ),
                shiny::column(
                  width = 6,
                  style = "padding: 0; border: 1px solid lightgray; height: 350px;",
                  networkD3::forceNetworkOutput("across_bug_ego_network", height = "100%")
                )
              )
            )
          ),
          shiny::tabPanel(
            "Across drugs",
            value = "across_drug",
            shiny::tagList(
              shiny::fluidRow(
                shiny::column(
                  width = 6,
                  style = "height: 50px; display: flex; align-items: center;",
                  shiny::h4("Toggle the button to select drug or drug class")
                ),
                shiny::column(
                  width = 4,
                  style = "height: 50px; display: flex; align-items: center;",
                  shiny::radioButtons(
                    inputId = "across_drug_id",
                    label = "",
                    choices = c("Drug class" = "drug_class", "Drug" = "drug"),
                    selected = "drug",
                    inline = TRUE
                  )
                )
              ),
              shiny::fluidRow(
                # Bug selector is the same in both modes, so define it once
                # outside the conditional panels (avoids a duplicate input id).
                shiny::column(
                  width = 6, style = "padding: 0;",
                  amr_select("bug_search_amr_across_drug", "Bug", character(0), selected = NULL, multiple = FALSE)
                ),
                shiny::conditionalPanel(
                  condition = "input.across_drug_id == 'drug'",
                  shiny::column(
                    width = 6, style = "padding: 0;",
                    amr_select("amr_drug_across_drug", "Drug (multi-select)", NULL)
                  )
                ),
                shiny::conditionalPanel(
                  condition = "input.across_drug_id == 'drug_class'",
                  shiny::column(
                    width = 6, style = "padding: 0;",
                    amr_select("amr_drug_class_across_drug", "Drug class (multi-select)", NULL)
                  )
                )
              ),
              shiny::column(
                width = 6, style = "padding: 0;",
                plotly::plotlyOutput("across_drug_feature_importance_plot", height = "100%"),
                style = "padding: 0; height: 600px; border: 1px solid lightgray;"
              ),
              shiny::column(
                width = 6,
                style = "padding: 0; height: 600px; border: 1px solid lightgray;",
                DT::dataTableOutput("across_drug_feature_importance_table", height = "100%")
              ),
              shiny::fluidRow(
                style = "padding: 10px 0;",
                shiny::column(
                  width = 6,
                  style = "padding: 0; border: 1px solid lightgray; height: 350px;",
                  plotly::plotlyOutput("across_drug_cog_barplot", height = "100%")
                ),
                shiny::column(
                  width = 6,
                  style = "padding: 0; border: 1px solid lightgray; height: 350px;",
                  networkD3::forceNetworkOutput("across_drug_ego_network", height = "100%")
                )
              )
            )
          )
        )
      )
    )
  )
}
