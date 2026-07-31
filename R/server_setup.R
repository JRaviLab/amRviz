# Shared reactives and cross-tab observers wired in once by setupServerCore().

setupServerCore <- function(input, output, session, results_root) {
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
      inputId = "results_species_dirs",
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

  # Performance metrics + top features, refreshed when the species selection changes.
  queryData <- shiny::reactiveVal(tibble::tibble())
  topFeatures <- shiny::reactiveVal(tibble::tibble())

  shiny::observeEvent(
    list(results_root, input$results_species_dirs),
    {
      if (!is.null(results_root) && nzchar(results_root)) {
        shiny::req(input$results_species_dirs)
      }

      queryData(loadMLResults(
        results_root = results_root,
        species_dirs = input$results_species_dirs
      ))
      topFeatures(loadTopFeat(
        results_root = results_root,
        species_dirs = input$results_species_dirs
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

  # Independent species list from metadata parquets so the Metadata tab works without ML results.
  available_metadata_species <- shiny::reactive({
    rr <- results_root
    scan_dirs <- if (!is.null(rr) && nzchar(rr)) {
      list.dirs(rr, full.names = TRUE, recursive = FALSE)
    } else {
      extdata <- system.file("extdata", package = "amRviz")
      if (nzchar(extdata)) {
        list.dirs(extdata, full.names = TRUE, recursive = FALSE)
      } else {
        character(0)
      }
    }
    # Each species subdirectory holds a single metadata.parquet; identify the
    # species by its subdirectory name.
    choices <- character(0)
    for (d in scan_dirs) {
      if (!file.exists(file.path(d, "metadata.parquet"))) next
      label <- gsub("_", " ", basename(d))
      choices <- c(choices, stats::setNames(basename(d), label))
    }
    sort(choices)
  })

  shiny::observe({
    choices <- available_species()
    sel <- if (length(choices)) choices[[1]] else NULL

    shiny::updateSelectizeInput(session, "bug_perf_id",
      choices = choices, selected = choices
    )
    shiny::updateSelectizeInput(session, "bug_search_amr_across_bug",
      choices = choices, selected = choices
    )
    shiny::updateSelectizeInput(session, "bug_search_amr_across_drug",
      choices = choices, selected = sel
    )
    shiny::updateSelectizeInput(session, "bug_cross_model_comparison_id",
      choices = choices, selected = sel
    )
    shiny::updateSelectizeInput(session, "network_bug_id",
      choices = choices, selected = sel
    )
  })

  shiny::observe({
    choices <- available_metadata_species()
    sel <- if (length(choices)) choices[[1]] else NULL
    shiny::updateSelectizeInput(session, "bug_metadata_id",
      choices = choices, selected = sel
    )
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

  # Normalised species code from the across-bug selector.
  bug_norm_input <- shiny::reactive({
    shiny::req(input$bug_search_amr_across_bug)
    normalize_species(input$bug_search_amr_across_bug)
  })

  list(
    queryData = queryData,
    topFeatures = topFeatures,
    available_species = available_species,
    available_metadata_species = available_metadata_species,
    filtered_top_features = filtered_top_features,
    loadDrugClassMapRec = loadDrugClassMapRec,
    bug_norm_input = bug_norm_input
  )
}
