#' Feature Importance UI Tab
#' @return A tabPanel for feature importance visualization
#' @keywords internal
featureImportanceUI <- function() {
  tabPanel(
    title = "Bug/Drug feature comparison",
    icon = icon("bug"),
    fluidRow(
      h3("Feature Importance Comparison", style = "margin-top: 15px; margin-bottom: 15px; font-weight: bold;"),
      style = "padding: 10px;",
      column(
        width = 4,
        style = "height: 80px; display: flex; align-items: center;",
        selectInput(
          "bug_drug_comp_model_scale",
          label = tags$label("Model scale", style = "font-size: 15px;"),
          choices = c("genes", "domains", "proteins"),
          multiple = F,
          selectize = TRUE,
          selected = "genes"
        )
      ),
      column(
        width = 4,
        style = "height: 80px; display: flex; align-items: center;",
        selectInput(
          "data_type",
          label = tags$label("Data type", style = "font-size: 15px;"),
          choices = c("count" = "counts", "binary" = "binary"),
          multiple = F,
          selectize = TRUE,
          selected = c("counts", "binary")
        )
      ),
      column(
        width = 4,
        style = "height: 80px; display: flex; align-items: center;",
        div(
          sliderInput(inputId = "top_n_features", label = "Top features (per model)", min = 0, max = 100, value = 10)
        )
      )
    ),
    column(
      width = 12, offset = 0,
      mainPanel(
        width = 12,
        style = "padding: 0;",
        tabsetPanel(
          id = "feature_importance_tabset",
          type = "tabs",
          tabPanel(
            "Across bugs",
            value = "across_bug",
            tagList(
              fluidRow(
                column(
                  width = 6,
                  style = "height: 50px; display: flex; align-items: center;",
                  h4("Toggle the button to select drug or drug class")
                ),
                column(
                  width = 4,
                  style = "height: 50px; display: flex; align-items: center;",
                  radioButtons(
                    inputId = "across_bug_id",
                    label = "",
                    choices = c("Drug class" = "drug_class", "Drug" = "drug"),
                    selected = "drug",
                    inline = TRUE
                  )
                )
              ),
              fluidRow(
                conditionalPanel(
                  condition = "input.across_bug_id == 'drug'",
                  column(
                    width = 6, style = "padding: 0;",
                    amr_select("bug_search_amr_across_bug", "Bug (multi-select)", character(0), selected = NULL)
                  ),
                  column(
                    width = 6, style = "padding: 0;",
                    amr_select("amr_drug_ml_across_bug", "Drug", NULL, FALSE)
                  )
                ),
                conditionalPanel(
                  condition = "input.across_bug_id == 'drug_class'",
                  column(
                    width = 6, style = "padding: 0;",
                    amr_select("bug_search_amr_across_bug", "Bug (multi-select)", character(0), selected = NULL)
                  ),
                  column(
                    width = 6, style = "padding: 0;",
                    amr_select("amr_drug_class_ml_across_bug", "Drug class", NULL, FALSE)
                  )
                )
              ),
              column(
                width = 6,
                style = "padding: 0; height: 600px; border: 1px solid lightgray;",
                plotOutput("across_bug_feature_importance_plot", height = "100%", width = "100%")
              ),
              column(
                width = 6,
                style = "padding: 0; height: 600px; border: 1px solid lightgray;",
                DT::dataTableOutput("across_bug_feature_importance_table", height = "100%")
              )
            )
          ),
          tabPanel(
            "Across drugs",
            value = "across_drug",
            tagList(
              fluidRow(
                column(
                  width = 6,
                  style = "height: 50px; display: flex; align-items: center;",
                  h4("Toggle the button to select drug or drug class")
                ),
                column(
                  width = 4,
                  style = "height: 50px; display: flex; align-items: center;",
                  radioButtons(
                    inputId = "across_drug_id",
                    label = "",
                    choices = c("Drug class" = "drug_class", "Drug" = "drug"),
                    selected = "drug",
                    inline = TRUE
                  )
                )
              ),
              fluidRow(
                conditionalPanel(
                  condition = "input.across_drug_id == 'drug'",
                  column(
                    width = 6, style = "padding: 0;",
                    amr_select("bug_search_amr_across_drug", "Bug", character(0), selected = NULL, multiple = FALSE)
                  ),
                  column(
                    width = 6, style = "padding: 0;",
                    amr_select("amr_drug_ml_across_drug", "Drug (multi-select)", NULL)
                  )
                ),
                conditionalPanel(
                  condition = "input.across_drug_id == 'drug_class'",
                  column(
                    width = 6, style = "padding: 0;",
                    amr_select("bug_search_amr_across_drug", "Bug", character(0), selected = NULL, multiple = FALSE)
                  ),
                  column(
                    width = 6, style = "padding: 0;",
                    amr_select("amr_drug_class_ml_across_drug", "Drug class (multi-select)", NULL)
                  )
                )
              ),
              column(
                width = 6, style = "padding: 0;",
                plotOutput("across_drug_feature_importance_plot", height = "100%"),
                style = "padding: 0; height: 600px; border: 1px solid lightgray;"
              ),
              column(
                width = 6,
                style = "padding: 0; height: 600px; border: 1px solid lightgray;",
                DT::dataTableOutput("across_drug_feature_importance_table", height = "100%")
              )
            )
          )
        )
      )
    )
  )
}
