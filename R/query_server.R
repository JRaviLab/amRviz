# Query-data tab: table views and CSV downloads.
#
# Both tables, and their downloads, are built by queryTable().


# Columns never shown: where files sit on the computer running the app, and
# plot styling columns from older amRviz.
.QUERY_HIDDEN <- c("metadata_path", "model_shapes", "model_colors")

# Columns shown before the user picks any, when present.
.QUERY_DEFAULT <- list(
  performance = c(
    "species_label", "drug_label", "drug_or_class", "feature_type",
    "feature_subtype", "strat_label", "strat_value", "strat_value_test",
    "seed", "model", "mcc", "bal_acc", "f1", "auprc"
  ),
  features = c(
    "species_label", "drug_label", "drug_or_class", "feature_type",
    "feature_subtype", "strat_label", "strat_value", "seed",
    "Variable", "Importance", "Sign"
  )
)


#' Columns a Query Data table can show
#'
#' @param data Performance or top-features tibble.
#' @return The column names, without the hidden ones.
#' @keywords internal
#' @noRd
queryColumns <- function(data) {
  setdiff(names(data), .QUERY_HIDDEN)
}


#' A Query Data table: the chosen columns of the loaded results
#'
#' @param data Performance or top-features tibble.
#' @param columns Columns to keep, in the order chosen; NULL or empty keeps
#'   every column that can be shown. Columns not in `data` are ignored.
#' @return `data` with those columns.
#' @keywords internal
#' @noRd
queryTable <- function(data, columns = NULL) {
  available <- queryColumns(data)
  keep <- if (length(columns)) intersect(columns, available) else available
  data[, keep, drop = FALSE]
}


#' Query Data tab server
#'
#' @param id Module id; must match the id passed to `queryDataUI()`.
#' @param core Shared data from setupServerCore().
#' @return NULL, invisibly.
#' @keywords internal
#' @noRd
queryDataServer <- function(id, core) {
  shiny::moduleServer(id, function(input, output, session) {
    # Column pickers: every column of the loaded results, starting with the
    # most useful ones selected.
    fill_columns <- function(input_id, data, defaults) {
      choices <- queryColumns(data)
      shiny::updateSelectizeInput(session, input_id,
        choices = choices, selected = intersect(defaults, choices)
      )
    }
    shiny::observe({
      fill_columns("query_data_columns", core$queryData(), .QUERY_DEFAULT$performance)
    })
    shiny::observe({
      fill_columns("top_features_columns", core$topFeatures(), .QUERY_DEFAULT$features)
    })

    output$queryDataTable <- DT::renderDataTable({
      queryTable(core$queryData(), input$query_data_columns)
    })

    output$query_data_download <- shiny::downloadHandler(
      filename = function() paste0("query_data_", Sys.Date(), ".csv"),
      content = function(file) {
        utils::write.csv(
          queryTable(core$queryData(), input$query_data_columns),
          file,
          row.names = FALSE
        )
      }
    )

    output$topFeaturesTable <- DT::renderDataTable({
      queryTable(core$topFeatures(), input$top_features_columns)
    })

    output$top_features_download <- shiny::downloadHandler(
      filename = function() paste0("top_features_", Sys.Date(), ".csv"),
      content = function(file) {
        utils::write.csv(
          queryTable(core$topFeatures(), input$top_features_columns),
          file,
          row.names = FALSE
        )
      }
    )

    invisible(NULL)
  })
}
