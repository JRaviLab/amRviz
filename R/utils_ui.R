# Small reusable UI input/control builders.


#' Padded action button with a styled label and icon
#'
#' @param id Input ID for the button.
#' @param label Button label text.
#' @param icon_name Name of the Font Awesome icon to show.
#' @param class_name CSS class applied to the button.
#' @return A shiny `div` wrapping the action button.
#' @keywords internal
#' @noRd
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


#' UI-output container with consistent box padding
#'
#' @param outputId Output ID to render into.
#' @return A shiny `uiOutput` with a padded container.
#' @keywords internal
#' @noRd
styledBox <- function(outputId) {
  uiOutput(outputId, container = function(...) {
    div(style = "display:inline-block")
    div(style = "padding-top: 0px; padding-bottom: 10px; height: 80%", ...)
  })
}


#' Padded selectize input with a styled label
#'
#' @param id Input ID for the select control.
#' @param label Label text.
#' @param choices Choices passed to `selectInput()`.
#' @param multiple Whether multiple selections are allowed.
#' @param selected Initially selected value(s), or NULL.
#' @return A shiny `div` wrapping the select input.
#' @keywords internal
#' @noRd
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
