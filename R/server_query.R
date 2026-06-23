# Query-data tab: table views and CSV downloads.

serverQuery <- function(input, output, session, core) {
  internal_cols <- c("model_shapes", "model_colors")

  output$queryDataTable <- DT::renderDataTable({
    data <- core$queryData()
    valid_cols <- setdiff(names(data), internal_cols)
    cols <- if (is.null(input$query_data_columns)) {
      valid_cols
    } else {
      intersect(valid_cols, input$query_data_columns)
    }
    data[, cols, drop = FALSE]
  })

  output$query_data_download <- shiny::downloadHandler(
    filename = function() paste0("query_data_", Sys.Date(), ".csv"),
    content = function(file) {
      data <- core$queryData()
      valid_cols <- setdiff(names(data), internal_cols)
      cols <- if (is.null(input$query_data_columns)) {
        valid_cols
      } else {
        intersect(valid_cols, input$query_data_columns)
      }
      utils::write.csv(data[, cols, drop = FALSE], file, row.names = FALSE)
    }
  )

  output$topFeaturesTable <- DT::renderDataTable({
    data <- core$topFeatures()
    cols <- input$top_features_columns %||% names(data)
    data[, cols, drop = FALSE]
  })

  output$top_features_download <- shiny::downloadHandler(
    filename = function() paste0("top_features_", Sys.Date(), ".csv"),
    content = function(file) {
      data <- core$topFeatures()
      cols <- input$top_features_columns %||% names(data)
      utils::write.csv(data[, cols, drop = FALSE], file, row.names = FALSE)
    }
  )
}
