# Shared data for the tab modules: the results-folder selector and the loaded
# ML results every model tab reads.


#' Results-folder selector UI
#'
#' Placeholder for the species result-folder selector, which is only shown when
#' the app is given a `results_root`. Placed inside the Query Data tab by
#' launchAMRDashboard().
#'
#' @param id Module id; must match the id passed to `resultsSourceServer()`.
#' @return A shiny `uiOutput`.
#' @keywords internal
#' @noRd
resultsSourceUI <- function(id) {
  shiny::uiOutput(shiny::NS(id, "results_species_selector_ui"))
}


#' Results-folder selector server
#'
#' Lists the amRml species folders under `results_root` and lets the user pick
#' which to load.
#'
#' @param id Module id; must match the id passed to `resultsSourceUI()`.
#' @param results_root User results root, or NULL for the packaged demo data.
#' @return A reactive returning the selected species-folder paths (NULL until
#'   the user's selection is available, and always NULL in demo mode).
#' @keywords internal
#' @noRd
resultsSourceServer <- function(id, results_root) {
  shiny::moduleServer(id, function(input, output, session) {
    output$results_species_selector_ui <- shiny::renderUI({
      if (is.null(results_root) || !nzchar(results_root)) {
        return(NULL)
      }

      choices <- listAmRmlSpeciesFolders(results_root)
      if (!length(choices)) {
        return(shiny::div(
          class = "note-box",
          paste("No amRml species folders found under:", results_root)
        ))
      }

      shiny::selectizeInput(
        inputId = session$ns("results_species_dirs"),
        label = "Select species result folders",
        choices = choices,
        selected = unname(choices),
        multiple = TRUE,
        options = list(
          plugins = list("remove_button"),
          placeholder = "Select species..."
        )
      )
    })

    shiny::reactive(input$results_species_dirs)
  })
}


#' Load the ML results shared by the model tabs
#'
#' @param results_root User results root, or NULL for the packaged demo data.
#' @param species_dirs Reactive returning the selected species-folder paths,
#'   from resultsSourceServer().
#' @return A named list of reactives: `queryData`, `topFeatures`,
#'   `available_species`, `filtered_top_features`, `loadDrugClassMapRec`.
#' @keywords internal
#' @noRd
setupServerCore <- function(results_root, species_dirs) {
  # Performance metrics + top features, refreshed when the species selection changes.
  queryData <- shiny::reactiveVal(tibble::tibble())
  topFeatures <- shiny::reactiveVal(tibble::tibble())

  shiny::observeEvent(
    list(results_root, species_dirs()),
    {
      if (!is.null(results_root) && nzchar(results_root)) {
        shiny::req(species_dirs())
      }

      queryData(loadMLResults(
        results_root = results_root,
        species_dirs = species_dirs()
      ))
      topFeatures(loadTopFeat(
        results_root = results_root,
        species_dirs = species_dirs()
      ))
    },
    ignoreInit = FALSE
  )

  # Species choices from ML perf data; excludes "cross"/"MDR" pseudo-species.
  available_species <- shiny::reactive({
    df <- queryData()
    if (is.null(df) || !nrow(df)) {
      return(character(0))
    }
    if (!all(c("species", "species_label") %in% names(df))) {
      return(character(0))
    }
    pairs <- df |>
      dplyr::filter(!(.data$species %in% c("cross", "MDR"))) |>
      dplyr::distinct(.data$species, .data$species_label)
    choices <- pairs$species
    names(choices) <- gsub("_", " ", pairs$species_label)
    sort(choices)
  })

  # Baseline (non-stratified, non-cross-test) slice consumed by the FI plots and tables.
  filtered_top_features <- shiny::reactive({
    tf <- topFeatures()
    if (is.null(tf) || !nrow(tf)) {
      return(tibble::tibble())
    }
    tf |>
      dplyr::filter(is.na(.data$strat_label) | !nzchar(.data$strat_label)) |>
      dplyr::filter(!isTRUE(.data$cross_test))
  })

  loadDrugClassMapRec <- shiny::reactive({
    loadDrugClassMap() |>
      dplyr::select(drug.antibiotic_name, drug_class) |>
      dplyr::distinct()
  })

  list(
    queryData = queryData,
    topFeatures = topFeatures,
    available_species = available_species,
    filtered_top_features = filtered_top_features,
    loadDrugClassMapRec = loadDrugClassMapRec
  )
}
