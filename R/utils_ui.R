# Small reusable UI input/control builders.



amr_button <- function(id, label, icon_name, class_name) {
  div(
    style = "padding: 5px; text-align: center;",
    actionButton(
      inputId = id,
      label = tags$label(label, style = "font-size: 15px;"),
      icon = icon(icon_name),
      class = class_name,
      style = "width: 80%; font-weight: bold;"
    )
  )
}


styledBox <- function(outputId) {
  uiOutput(outputId, container = function(...) {
    div(style = "display:inline-block")
    div(style = "padding-top: 0px; padding-bottom: 10px; height: 80%", ...)
  })
}


# Select input helper
amr_select <- function(id, label, choices, multiple = TRUE, selected = NULL) {
  div(
    style = "padding: 10px;",
    selectInput(
      inputId = id,
      label = tags$label(label, style = "font-size: 15px;"),
      choices = choices,
      multiple = multiple,
      selectize = TRUE,
      width = "100%",
      selected = selected
    )
  )
}

