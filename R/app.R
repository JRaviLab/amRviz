#' AMR Shiny dashboard application
#'
#' Creates and returns a Shiny application for exploring antimicrobial
#' resistance data and machine learning model results.
#'
#' @param results_root File path to the root directory containing amRml model
#'        output results. If `NULL` (default), the application will attempt to
#'        use example data bundled with the package, where available.
#'
#' @return A Shiny application object
#' @export
#' @import shiny
#' @importFrom magrittr %>%
#' @importFrom utils head write.csv
#' @importFrom shinyjs useShinyjs
#' @examples
#' if (interactive()) {
#'   app <- launchAMRDashboard()
#'   shiny::runApp(app)
#' }
launchAMRDashboard <- function(results_root = NULL) {
  # Bug choices
  bug_choices <- c(
    "Enterococcus faecium" = "Efa",
    "Staphylococcus aureus" = "Sau",
    "Klebsiella pneumoniae" = "Kpn",
    "Acinetobacter baumannii" = "Aba",
    "Pseudomonas aeruginosa" = "Pae",
    "Enterobacter spp." = "Esp.",
    "Escherichia coli" = "Eco",
    "Campylobacter jejuni" = "Cje",
    "Staphylococcus epidermidis" = "Sep"
  )

  # ESKAPE bugs (sorted in ESKAPE order)
  eskape_bugs <- c(
    "Enterococcus faecium" = "Efa",
    "Staphylococcus aureus" = "Sau",
    "Klebsiella pneumoniae" = "Kpn",
    "Acinetobacter baumannii" = "Aba",
    "Pseudomonas aeruginosa" = "Pae",
    "Enterobacter spp." = "Esp.",
    "Escherichia coli" = "Eco",
    "Campylobacter jejuni" = "Cje",
    "Staphylococcus epidermidis" = "Sep",
    "Streptococcus pneumoniae" = "Spn"
  )

  # UI
  ui <- tagList(
    shinyjs::useShinyjs(),
    tags$head(includeCSS(system.file("app/www/style.css", package = "amRshiny"))),
    tags$head(
      tags$style(HTML("
                      .innerbox {
                        /*border: 2px solid black;*/
                        box-shadow: 2px 2px 3px 3px #ccc;
                        margin: auto;
                        padding: 20px;
                      }

                      .bord {
                        margin: auto;
                        padding: 20px;
                      }

                      .lightblue-link{
                        color:#11aad9;
                      }

                      .iMargin{
                        margin: 5px;
                      }

                     a{
                        color: #141414;
                        text-decoration:none;
                      }


                      a:visited{
                        color:none;
                      }

                      .noDec{
                        color:white;
                        text-decoration:none;
                      }

                      .zoom:hover {
                        /*color:white;*/
                        transform: scale(1.3);
                      }

                      .note-box {
                        margin: auto;
                        margin-top: 10px;
                        margin-bottom: 10px;
                        width: 80%;
                        border: 1px solid #78adff;
                        padding: 10px;
                        background-color: #b7dceb;
                        border-radius: 5px;
                        text-align: center;
                      }
                    "))
    ),
    navbarPage(
      id = "tabselected",
      selected = "home",
      title = div("AMR"), # not clickable
      # 2. Home icon tab (right of AMR dashboard)
      tabPanel(
        title = icon("home", class = "home-tab-icon"),
        value = "home",
        fluidPage(
          h2("amR: an R package suite to predict antimicrobial resistance in bacterial pathogens"),
          tags$p(
            tags$strong("Authors: "),
            "Evan P Brenner,#, Abhirupa Ghosh,#, Ethan P Wolfe, Emily A Boyer, Charmie K Vang, Raymond L Lesiyon, David Mayer, Janani Ravi*."
          ),
          tags$p(
            "Department of Biomedical Informatics, Center for Health Artificial Intelligence, University of Colorado Anschutz, Aurora, CO, USA",
            tags$span(style = "font-style: italic", "#Co-primary, contributed equally. *Corresponding author: janani.ravi@cuanschutz.edu")
          ),
          br(),
          h4("Abstract"),
          tags$strong("Motivation: "),
          tags$p("Identifying antimicrobial resistance (AMR) in bacterial pathogens is critical for diagnostics and treatment, but resistance is a complex trait arising from diverse mechanisms spanning multiple molecular scales. Existing computational approaches often function as black boxes and rarely explore cross-species or multi-drug patterns. We developed amR, an integrated R package suite providing a complete framework from bacterial genome sequences to interpretable AMR predictions, enabling identification of resistance mechanisms across species and drugs."),
          tags$strong("Results: "),
          tags$p("The amR suite consists of three modular packages. amRdata interfaces with BV-BRC to download and process bacterial genomes with paired antimicrobial susceptibility testing data, constructs species-specific pangenomes, and extracts features at four molecular scales: gene/protein clusters, protein domains, and structural variants. All data are stored in efficient Parquet and DuckDB formats. amRml trains interpretable logistic regression machine learning models per species-drug combination, generating ranked features by importance and comprehensive performance metrics (balanced accuracy, F1, MCC). Models identify known resistance determinants (e.g., gyrA mutations for fluoroquinolones, mecA for beta-lactams) alongside poorly characterized features representing potential novel mechanisms. amRshiny provides an interactive Shiny dashboard to explore isolate metadata distributions, compare model performance across species and drugs, visualize top predictive AMR features, and analyze cross-model patterns (including features specific to geographic/temporal strata). The suite has been applied to ESKAPE pathogens, achieving balanced accuracies above 0.80. With thousands of genomes, multi-scale features, and interpretable models, amR provides the first comprehensive programmatic framework and R package for AMR research."),
          tags$strong("Availability and implementation: "),
          tags$p("https://github.com/JRaviLab/amR"),
          tags$hr(),
        )
      ),
      # other tabs
      metadataUI(),
      modelPerfUI(),
      featureImportanceUI(),
      crossModelComparisonUI(),
      queryDataUI()
    ),
    tags$footer(
      class = "footer",
      tags$div(
        "(c) JRaviLab 2026 | ",
        tags$a(href = "https://jravilab.github.io", "jravilab.github.io"),
        " | ",
        tags$a(href = "https://twitter.com/jravilab", "@jravilab"),
        " | janani.ravi@cuanschutz.edu |",
        tags$a(
          href = "https://docs.google.com/forms/d/e/1FAIpQLSdwFo5Wwt_t4WGthDGgc1EYhvvKagUEb3RiNLdsbnpDlYTk7Q/viewform?usp=dialog",
          icon("question-circle"),
          " Help Doc"
        )
      )
    ),
  )
  ## Server
  server <- function(input, output, session) {
    # Initialize reactive values
    genome_data <- reactiveVal(NULL)
    drug_class_map <- reactiveVal(NULL)
    amr_drugs <- reactiveVal(NULL)
    ml_metrics_results <- reactiveVal(NULL)
    ml_ui_trigger <- reactiveVal(0)

    ## Render species selector ui
    output$results_species_selector_ui <- renderUI({
      if (is.null(results_root) || !nzchar(results_root)) {
        return(NULL)
      }

      choices <- listAmRmlSpeciesFolders(results_root)
      if (!length(choices)) {
        return(div(
          class = "note-box",
          paste("No amRml species folders found under:", results_root)
        ))
      }

      selectizeInput(
        inputId = "results_species_dirs",
        label = "Select species result folders",
        choices = choices, # label=foldername, value=full path
        selected = unname(choices), # default select all
        multiple = TRUE,
        options = list(plugins = list("remove_button"), placeholder = "Select species...")
      )
    })

    ## Query Data Sources
    # Load performance metrics + top features based on selected species folders
    queryData <- reactiveVal(tibble::tibble())
    topFeatures <- reactiveVal(tibble::tibble())

    observeEvent(
      list(results_root, input$results_species_dirs),
      {
        # If using real results_root, wait until the user has a selection
        if (!is.null(results_root) && nzchar(results_root)) {
          req(input$results_species_dirs)
        }

        perf <- loadMLResults(
          results_root = results_root,
          species_dirs = input$results_species_dirs
        )
        top <- loadTopFeat(
          results_root = results_root,
          species_dirs = input$results_species_dirs
        )

        queryData(perf)
        topFeatures(top)
      },
      ignoreInit = FALSE
    )

    # Derive species choices from loaded ML perf data (excludes pseudo-species like "cross").
    # names = human-readable label (underscores → spaces, e.g. "Shigella flexneri")
    # values = 3-letter species code (e.g. "Sfl")
    available_species <- reactive({
      df <- queryData()
      if (is.null(df) || !nrow(df)) {
        return(character(0))
      }
      if (!all(c("species", "species_label") %in% names(df))) {
        return(character(0))
      }
      pairs <- df %>%
        dplyr::filter(!(.data$species %in% c("cross", "MDR"))) %>%
        dplyr::distinct(.data$species, .data$species_label)
      choices <- pairs$species
      names(choices) <- gsub("_", " ", pairs$species_label)
      sort(choices)
    })

    # Derive species choices from available metadata parquets (independent of ML perf data).
    # Scans for *_metadata.parquet files; species code = filename prefix, label = directory name.
    available_metadata_species <- reactive({
      rr <- results_root
      scan_dirs <- if (!is.null(rr) && nzchar(rr)) {
        list.dirs(rr, full.names = TRUE, recursive = FALSE)
      } else {
        extdata <- system.file("extdata", package = "amRshiny")
        if (nzchar(extdata)) list.dirs(extdata, full.names = TRUE, recursive = FALSE) else character(0)
      }
      choices <- character(0)
      for (d in scan_dirs) {
        fps <- list.files(d, pattern = "_metadata\\.parquet$", full.names = FALSE)
        if (!length(fps)) next
        for (f in fps) {
          code <- sub("_metadata\\.parquet$", "", f)
          label <- gsub("_", " ", basename(d))
          choices <- c(choices, setNames(code, label))
        }
      }
      sort(choices)
    })

    # Update every bug selector whenever the available species change
    observe({
      choices <- available_species()
      sel <- if (length(choices)) choices[[1]] else NULL

      updateSelectizeInput(session, "bug_ml_perf_id",
        choices = choices, selected = choices
      ) # multi-select: all by default
      updateSelectizeInput(session, "bug_search_amr_across_bug",
        choices = choices, selected = choices
      )
      updateSelectizeInput(session, "bug_search_amr_across_drug",
        choices = choices, selected = sel
      )
    })

    # Update metadata bug selector from metadata parquets (not ML perf data)
    observe({
      choices <- available_metadata_species()
      sel <- if (length(choices)) choices[[1]] else NULL
      updateSelectizeInput(session, "bug_metadata_id",
        choices = choices, selected = sel
      )
    })

    # Reactive: filtered top features for the feature importance plots/tables
    # (replaces the reactiveFileReader anti-pattern that wrote/read temp TSVs)
    filtered_top_features <- reactive({
      tf <- topFeatures()
      if (is.null(tf) || !nrow(tf)) {
        return(tibble::tibble())
      }
      tf %>%
        dplyr::filter(is.na(.data$strat_label) | !nzchar(.data$strat_label)) %>%
        dplyr::filter(!isTRUE(.data$cross_test))
    })

    loadDrugClassMapRec <- reactive({
      loadDrugClassMap() %>%
        dplyr::select(drug.antibiotic_name, drug_class) %>%
        dplyr::distinct()
    })

    # reactive normalized UI selection for bug_search_amr_across_bug
    bug_norm_input <- reactive({
      req(input$bug_search_amr_across_bug) # ensure value exists
      normalize_species(input$bug_search_amr_across_bug)
    })

    # dynamic ML ui
    output$ml_drug_toggle_ui <- renderUI({
      # For example: only allow selection among the *other two* categories
      choices <- switch(input$ml_drug_toggle_top,
        "bug" = c("Drug Class" = "class", "Drug" = "drug"),
        "class" = c("Bug" = "bug", "Drug" = "drug"),
        "drug" = c("Bug" = "bug", "Drug Class" = "class")
      )

      radioButtons(
        inputId = "ml_drug_toggle",
        label = tags$label("Select by", style = "font-size: 15px;"),
        choices = choices,
        selected = names(choices)[1],
        inline = TRUE
      )
    })

    # observe to load all the data
    observe({
      drug_class_map(loadDrugClassMapRec())
    })

    # clear the selected bug when the user changes the drug/drug class
    observeEvent(input$across_bug_id, {
      updateSelectInput(
        session,
        inputId = "bug_search_amr_across_bug",
        choices = bug_choices,
        selected = unname(bug_choices) # Default to ESKAPE bugs
      )
    })

    observeEvent(c(input$bug_search_amr_across_bug, input$across_bug_id, input$bug_drug_comp_model_scale, input$data_type), {
      req(input$bug_search_amr_across_bug)
      message("Across bug feature comparison: ")

      bn <- bug_norm_input()

      base_tf <- topFeatures() |>
        dplyr::filter(normalize_species(.data$species) %in% bn) |>
        dplyr::filter(.data$feature_type %in% input$bug_drug_comp_model_scale) |>
        dplyr::filter(.data$feature_subtype %in% input$data_type) |>
        dplyr::filter(is.na(.data$strat_label) | !nzchar(.data$strat_label))

      if (identical(input$across_bug_id, "drug")) {
        drugs_vec <- base_tf |>
          dplyr::filter(.data$drug_label == "drug") |>
          dplyr::pull(.data$drug_or_class) |>
          unique() |>
          sort()
        sel <- intersect("GEN", drugs_vec)
        if (!length(sel)) sel <- head(drugs_vec, 1)
        updateSelectInput(session, "amr_drug_ml_across_bug", choices = drugs_vec, selected = sel)
      } else {
        drugs_class_vec <- base_tf |>
          dplyr::filter(.data$drug_label == "drug_class") |>
          dplyr::pull(.data$drug_or_class) |>
          unique() |>
          sort()
        sel <- head(drugs_class_vec, 1)
        updateSelectInput(session, "amr_drug_class_ml_across_bug", choices = drugs_class_vec, selected = sel)
      }
    })


    # For drug/drug class (across drug) | filter by scale/data type
    observeEvent(c(input$bug_search_amr_across_drug, input$across_drug_id, input$bug_drug_comp_model_scale, input$data_type), {
      req(input$bug_search_amr_across_drug)
      message("Across drug feature comparison: ")

      bn <- normalize_species(input$bug_search_amr_across_drug)

      base_tf <- topFeatures() |>
        dplyr::filter(normalize_species(.data$species) %in% bn) |>
        dplyr::filter(.data$feature_type %in% input$bug_drug_comp_model_scale) |>
        dplyr::filter(.data$feature_subtype %in% input$data_type) |>
        dplyr::filter(is.na(.data$strat_label) | !nzchar(.data$strat_label))

      if (identical(input$across_drug_id, "drug")) {
        drugs_vec <- base_tf |>
          dplyr::filter(.data$drug_label == "drug") |>
          dplyr::pull(.data$drug_or_class) |>
          unique() |>
          sort()
        prev <- isolate(input$amr_drug_ml_across_drug)
        sel <- prev[prev %in% drugs_vec]
        if (!length(sel)) {
          pref <- c("OXA", "PEN", "MET")
          sel <- intersect(pref, drugs_vec)
          if (!length(sel)) sel <- head(drugs_vec, min(3, length(drugs_vec)))
        }
        updateSelectInput(session, "amr_drug_ml_across_drug", choices = drugs_vec, selected = sel)
      } else {
        drugs_class_vec <- base_tf |>
          dplyr::filter(.data$drug_label == "drug_class") |>
          dplyr::pull(.data$drug_or_class) |>
          unique() |>
          sort()
        prev <- isolate(input$amr_drug_class_ml_across_drug)
        sel <- prev[prev %in% drugs_class_vec]
        if (!length(sel)) {
          pref <- c("CEP", "LIN", "MAC", "PEN")
          sel <- intersect(pref, drugs_class_vec)
          if (!length(sel)) sel <- head(drugs_class_vec, min(4, length(drugs_class_vec)))
        }
        updateSelectInput(session, "amr_drug_class_ml_across_drug", choices = drugs_class_vec, selected = sel)
      }
    })

    # Update the model performance drug class dropdown when bug selection changes
    observeEvent(input$bug_ml_perf_id, {
      data <- queryData() %>%
        dplyr::filter(normalize_species(.data$species) %in% normalize_species(input$bug_ml_perf_id)) %>%
        dplyr::filter(is.na(.data$strat_label) | !nzchar(.data$strat_label))
      drug_class_vec <- data %>%
        dplyr::filter(.data$drug_label == "drug_class") %>%
        dplyr::pull(.data$drug_or_class) %>%
        unique() %>%
        sort()

      sel <- if ("AMG" %in% drug_class_vec) "AMG" else if (length(drug_class_vec)) drug_class_vec[1] else "all"
      updateSelectInput(
        session,
        inputId = "drug_class_ml_perf_id",
        choices = c("all", drug_class_vec),
        selected = sel
      )
    })
    # model holdouts filtering

    observeEvent(input$bug_holdouts_id,
      {
        req(input$bug_holdouts_id)

        # Build choices filtered to the selected 3-letter species code
        choices <- getHoldoutsDrugChoices(perf_data = queryData(), bug = input$bug_holdouts_id)

        # Keep the user's current selection if still valid; otherwise pick first
        prev <- isolate(input$holdouts_drug)
        sel <- if (!is.null(prev) && prev %in% choices) prev else if (length(choices)) choices[[1]] else NULL

        # Update the dropdown (use updateSelectizeInput if your UI uses selectize=TRUE)
        updateSelectizeInput(
          session,
          inputId  = "holdouts_drug",
          choices  = choices,
          selected = sel,
          server   = TRUE
        )
      },
      ignoreInit = FALSE
    ) # run once on app load so it populates immediately


    observeEvent(input$drug_class_ml_perf_id, {
      req(input$drug_class_ml_perf_id)
      base_data <- queryData() %>%
        dplyr::filter(normalize_species(.data$species) %in% normalize_species(input$bug_ml_perf_id)) %>%
        dplyr::filter(is.na(.data$strat_label) | !nzchar(.data$strat_label))
      drug_vec <- base_data %>%
        dplyr::filter(.data$drug_label == "drug") %>%
        dplyr::pull(.data$drug_or_class) %>%
        unique() %>%
        sort()

      if (input$drug_class_ml_perf_id != "all") {
        # Use metadata to map class abbreviation -> drug abbreviations
        sp_codes <- normalize_species(input$bug_ml_perf_id)
        meta <- dplyr::bind_rows(lapply(sp_codes, function(sp) {
          fp <- get_metadata_path(sp, results_root)
          if (!is.null(fp)) .read_parquet_safe(fp, verbose = FALSE) else tibble::tibble()
        }))
        if (nrow(meta) && all(c("class_abbr", "drug_abbr") %in% names(meta))) {
          drugs_in_class <- meta %>%
            dplyr::filter(.data$class_abbr %in% input$drug_class_ml_perf_id) %>%
            dplyr::pull(.data$drug_abbr) %>%
            unique()
          drug_vec <- intersect(drug_vec, drugs_in_class)
        }
        updateSelectInput(
          session,
          inputId = "drug_ml_perf_id",
          choices = drug_vec,
          selected = if (length(drug_vec)) drug_vec[1] else NULL
        )
      } else {
        updateSelectInput(
          session,
          inputId = "drug_ml_perf_id",
          choices = drug_vec,
          selected = if ("GEN" %in% drug_vec) "GEN" else if (length(drug_vec)) drug_vec[1] else NULL
        )
      }
    })


    ## plot UI
    output$geo_isolate_plot_ui <- renderUI({
      fluidRow(
        div(
          box(
            title = "",
            width = 12,
            plotly::plotlyOutput("geo_isolate_plot")
          )
        )
      )
    })

    output$r_s_across_time_ui <- renderUI({
      fluidRow(
        box(
          title = "", # "Resistance vs Susceptible across time",
          width = 12,
          plotOutput("r_s_across_time_plot")
        )
      )
    })

    output$host_isolate_plot_ui <- renderUI({
      fluidRow(
        box(
          title = "", # "Host",
          width = 12,
          plotOutput("host_isolate_plot")
        )
      )
    })

    output$isolation_source_ui <- renderUI({
      fluidRow(
        box(
          title = "", # "Isolation source",
          width = 12,
          plotOutput("isolation_source_plot")
        )
      )
    })

    output$resistance_vs_susceptible_ui <- renderUI({
      fluidRow(
        box(
          title = "", # "Data availability",
          width = 12,
          plotOutput("resistance_vs_susceptible_plot")
        )
      )
    })

    ## get a quick summary plot;
    observeEvent(input$bug_metadata_id, {
      output$quick_metadata_stats <- renderUI({
        metadata <- purrr::map_dfr(
          .x = input$bug_metadata_id,
          .f = function(x) {
            fp <- get_metadata_path(x, results_root)
            if (!is.null(fp) && file.exists(fp)) {
              arrow::read_parquet(fp) |> dplyr::mutate(species = input$bug_metadata_id)
            } else {
              return(tibble())
            }
          }
        )
        makeQuickStats(metadata)
      })
    })

    ## making the data availability plot
    observe({
      req(input$bug_metadata_id)

      output$resistance_vs_susceptible_plot <- renderPlot({
        metadata <- purrr::map_dfr(
          .x = input$bug_metadata_id,
          .f = function(x) {
            fp <- get_metadata_path(x, results_root)
            if (!is.null(fp) && file.exists(fp)) {
              arrow::read_parquet(fp) |> dplyr::mutate(species = input$bug_metadata_id)
            } else {
              return(tibble())
            }
          }
        )
        makeDatAvailabilityPlot(metadata)
      })
    })

    output$geo_isolate_plot <- plotly::renderPlotly({
      data <- purrr::map_dfr(
        .x = input$bug_metadata_id,
        .f = function(x) {
          fp <- get_metadata_path(x, results_root)
          if (!is.null(fp) && file.exists(fp)) {
            arrow::read_parquet(fp) |> dplyr::mutate(species = input$bug_metadata_id)
          } else {
            return(tibble())
          }
        }
      )

      data <- data %>%
        dplyr::filter(genome.isolation_country != "") %>%
        dplyr::mutate(
          genome_drug.evidence = case_when(
            genome_drug.laboratory_typing_method %in%
              c(
                "Disk diffusion",
                "MIC",
                "Broth dilution",
                "Agar dilution"
              ) ~ "Laboratory Method",
            genome_drug.laboratory_typing_method ==
              "Computational Prediction" ~ "Computational Method",
            TRUE ~ genome_drug.laboratory_typing_method
          )
        ) %>%
        group_by(genome.isolation_country, genome_drug.antibiotic) %>%
        summarize(count = n(), .groups = "drop") %>%
        dplyr::group_by(genome.isolation_country) %>%
        dplyr::summarise(count = sum(count), .groups = "drop") %>%
        dplyr::collect()

      makeGeoChloroPlot(data)
    })

    output$r_s_across_time_plot <- renderPlot({
      data <- purrr::map_dfr(
        .x = input$bug_metadata_id,
        .f = function(x) {
          fp <- get_metadata_path(x, results_root)
          if (!is.null(fp) && file.exists(fp)) {
            arrow::read_parquet(fp) |> dplyr::mutate(species = input$bug_metadata_id)
          } else {
            return(tibble())
          }
        }
      )

      data <- data %>%
        dplyr::mutate(
          genome_drug.evidence = case_when(
            genome_drug.laboratory_typing_method %in%
              c(
                "Disk diffusion",
                "MIC",
                "Broth dilution",
                "Agar dilution"
              ) ~ "Laboratory Method",
            genome_drug.laboratory_typing_method ==
              "Computational Prediction" ~ "Computational Method",
            TRUE ~ genome_drug.laboratory_typing_method
          )
        ) %>%
        dplyr::filter(!is.na(genome.collection_year)) %>%
        group_by(
          genome_drug.antibiotic,
          genome_drug.resistant_phenotype,
          genome.isolation_country,
          genome.collection_year
        ) %>%
        summarize(n = n()) %>%
        dplyr::collect()
      # print(data)
      makeTimeSeriesAMRPlot(data, input$amr_drug_search)
    })

    output$host_isolate_plot <- renderPlot({
      data <- purrr::map_dfr(
        .x = input$bug_metadata_id,
        .f = function(x) {
          fp <- get_metadata_path(x, results_root)
          if (!is.null(fp) && file.exists(fp)) {
            arrow::read_parquet(fp) |> dplyr::mutate(species = input$bug_metadata_id)
          } else {
            return(tibble())
          }
        }
      )

      data <- data %>%
        dplyr::filter(genome.host_common_name != "") %>%
        dplyr::mutate(
          genome_drug.evidence = case_when(
            genome_drug.laboratory_typing_method %in%
              c(
                "Disk diffusion",
                "MIC",
                "Broth dilution",
                "Agar dilution"
              ) ~ "Laboratory Method",
            genome_drug.laboratory_typing_method ==
              "Computational Prediction" ~ "Computational Method",
            TRUE ~ genome_drug.laboratory_typing_method
          )
        ) %>%
        dplyr::collect()

      makeHostIsolatePlot(data)
    })

    output$isolation_source_plot <- renderPlot({
      data <- purrr::map_dfr(
        .x = input$bug_metadata_id,
        .f = function(x) {
          fp <- get_metadata_path(x, results_root)
          if (!is.null(fp) && file.exists(fp)) {
            arrow::read_parquet(fp) |> dplyr::mutate(species = input$bug_metadata_id)
          } else {
            return(tibble())
          }
        }
      )

      data <- data %>%
        dplyr::filter(genome.host_common_name != "") %>%
        dplyr::mutate(
          genome_drug.evidence = case_when(
            genome_drug.laboratory_typing_method %in%
              c(
                "Disk diffusion",
                "MIC",
                "Broth dilution",
                "Agar dilution"
              ) ~ "Laboratory Method",
            genome_drug.laboratory_typing_method ==
              "Computational Prediction" ~ "Computational Method",
            TRUE ~ genome_drug.laboratory_typing_method
          )
        ) %>%
        dplyr::collect()
      makeIsolationSourcesPlot(data)
    })


    ## ML Metrics
    # plotly::renderPlotly
    output$model_perfomance_plot <- renderPlot({
      makeModelPerformancePlot(
        queryData(),
        input$bug_ml_perf_id,
        input$model_scale,
        input$data_type,
        input$model_metrics,
        input$drug_class_ml_perf_id,
        input$drug_ml_perf_id
      )
    })

    observe({
      output$across_bug_feature_importance_plot <- renderPlot({
        if (is.null(input$across_bug_id)) {
          return(NULL)
        }
        amr_drug <- if (input$across_bug_id == "drug") {
          input$amr_drug_ml_across_bug
        } else {
          input$amr_drug_class_ml_across_bug
        }
        ht <- makeFeatureImportancePlot(
          topFeatures(),
          input$bug_search_amr_across_bug,
          amr_drug,
          input$bug_drug_comp_model_scale,
          input$data_type,
          input$top_n_features,
          input$feature_importance_tabset
        )
        draw(ht, heatmap_legend_side = "right")
      })
    })

    observe({
      output$across_drug_feature_importance_plot <- renderPlot({
        if (is.null(input$across_drug_id)) {
          return(NULL)
        }
        amr_drug <- if (input$across_drug_id == "drug") {
          input$amr_drug_ml_across_drug
        } else {
          input$amr_drug_class_ml_across_drug
        }
        ht <- makeFeatureImportancePlot(
          topFeatures(),
          input$bug_search_amr_across_drug,
          amr_drug,
          input$bug_drug_comp_model_scale,
          input$data_type,
          input$top_n_features,
          input$feature_importance_tabset
        )
        draw(ht, heatmap_legend_side = "right")
      })
    })

    observe({
      output$feature_importance_plot <- renderPlot({
        ht <- makeFeatureImportancePlot(
          topFeatures(),
          input$bug_ml_perf_id,
          input$amr_drug_ml_across_bug,
          input$model_scale,
          input$data_type,
          input$top_n_features,
          input$feature_importance_tabset
        )
        draw(ht, heatmap_legend_side = "right")
      })
      output$feature_importance_table <- DT::renderDataTable({
        tf <- filtered_top_features() %>%
          dplyr::filter(normalize_species(.data$species) %in% normalize_species(input$bug_ml_perf_id)) %>%
          dplyr::filter(.data$drug_or_class %in% input$amr_drug_ml_across_bug)
        if (!nrow(tf)) {
          return(NULL)
        }
        makeFeatureImportTable(tf)
      })
    })

    # Feature importance tables: across bug and across drug
    output$across_bug_feature_importance_table <- DT::renderDataTable({
      amr_drug <- if (!is.null(input$across_bug_id) && input$across_bug_id == "drug_class") {
        input$amr_drug_class_ml_across_bug
      } else {
        input$amr_drug_ml_across_bug
      }
      tf <- filtered_top_features() %>%
        dplyr::filter(normalize_species(.data$species) %in% normalize_species(input$bug_search_amr_across_bug)) %>%
        dplyr::filter(.data$drug_or_class %in% amr_drug)
      if (!nrow(tf)) {
        return(NULL)
      }
      makeFeatureImportTable(tf)
    })

    output$across_drug_feature_importance_table <- DT::renderDataTable({
      amr_drug <- if (!is.null(input$across_drug_id) && input$across_drug_id == "drug_class") {
        input$amr_drug_class_ml_across_drug
      } else {
        input$amr_drug_ml_across_drug
      }
      tf <- filtered_top_features() %>%
        dplyr::filter(normalize_species(.data$species) %in% normalize_species(input$bug_search_amr_across_drug)) %>%
        dplyr::filter(.data$drug_or_class %in% amr_drug)
      if (!nrow(tf)) {
        return(NULL)
      }
      makeFeatureImportTable(tf)
    })

    # Cross model feature importance table
    output$cross_model_feature_importance_table <- DT::renderDataTable({
      strat <- if (isTRUE(input$cross_model_comparison == "country")) "country" else "year"
      tf <- topFeatures() %>%
        dplyr::filter(normalize_species(.data$species) %in% normalize_species(input$bug_cross_model_comparison_id)) %>%
        dplyr::filter(.data$drug_or_class %in% input$drug_cross_model_comparison_id) %>%
        dplyr::filter(.data$strat_label == strat) %>%
        dplyr::filter(!isTRUE(.data$cross_test))
      if (!nrow(tf)) {
        return(NULL)
      }
      makeFeatureImportTable(tf)
    })

    # model comparisons;
    observeEvent(input$bug_cross_model_comparison_id,
      {
        # Gather Drug/Drug class options for stratified holdout models for this bug
        drugs_vec <- getHoldoutsDrugChoices(perf_data = queryData(), bug = input$bug_cross_model_comparison_id)

        sel <- if (length(drugs_vec)) drugs_vec[1] else NULL
        updateSelectInput(
          session,
          inputId = "drug_cross_model_comparison_id",
          choices = drugs_vec,
          selected = sel
        )
      },
      ignoreInit = FALSE
    )

    observe({
      output$cross_model_perf_plot <- renderPlot({
        ht <- makeCrossModelPerformancePlot(
          queryData(),
          input$bug_cross_model_comparison_id,
          input$drug_cross_model_comparison_id,
          input$cross_model_comparison
        )
        draw(ht, heatmap_legend_side = "left")
      })
      output$cross_model_feature_importance_plot <- renderPlot({
        ht <- makeCrossModelFeatureImportancePlot(
          topFeatures(),
          input$bug_cross_model_comparison_id,
          input$drug_cross_model_comparison_id,
          input$cross_model_comparison,
          input$cross_model_top_n_features
        )
        draw(ht, heatmap_legend_side = "left")
      })


      # # Load performance metrics data
      # queryData <- reactiveVal(loadMLResults(results_root = results_root))
      #
      # observe({
      #   data <- queryData() # Get data
      #
      #   # Exclude columns
      #   valid_columns <- setdiff(names(data), c("model_shapes", "model_colors"))
      #
      #   # Update dropdown options
      #   updateSelectizeInput(
      #     session,
      #     "query_data_columns",
      #     choices = valid_columns, # Get column names
      #     selected = valid_columns # Default selection: all columns
      #   )
      # })

      # Render the data table
      output$queryDataTable <- DT::renderDataTable({
        data <- queryData() # Get data

        # Exclude unwanted columns
        valid_columns <- setdiff(names(data), c("model_shapes", "model_colors"))

        selected_cols <- input$query_data_columns

        if (is.null(selected_cols)) {
          data[, valid_columns, drop = FALSE] # Return data
        } else {
          data[, intersect(valid_columns, selected_cols), drop = FALSE] # Return filtered data
        }
      })

      # Download filtered data
      output$query_data_download <- downloadHandler(
        filename = function() {
          paste("query_data_", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
          data <- queryData() # Get data
          # Exclude unwanted columns
          valid_columns <- setdiff(names(data), c("model_shapes", "model_colors")) # Filter columns

          selected_cols <- input$query_data_columns # Get user-selected columns

          if (is.null(selected_cols)) {
            write.csv(data[, valid_columns, drop = FALSE], file, row.names = FALSE) # Download only valid columns
          } else {
            write.csv(data[, intersect(valid_columns, selected_cols), drop = FALSE], file, row.names = FALSE) # Download user selection intersected with valid columns
          }
        }
      )
      # Top Features Table Logic
      # topFeatures <- reactiveVal(loadTopFeat(results_root = results_root)) # Use the loadTopFeat() function to load data
      #
      # observe({
      #   data <- topFeatures() # Fetch data
      #
      #   # Dynamically update dropdown with column names
      #   updateSelectizeInput(
      #     session,
      #     "top_features_columns",
      #     choices = names(data), # Populate dropdown with column names
      #     selected = names(data) # Default: all columns selected
      #   )
      # })

      output$topFeaturesTable <- DT::renderDataTable({
        data <- topFeatures() # Get the data

        selected_columns <- input$top_features_columns # Get user-selected columns

        if (is.null(selected_columns)) {
          data # Show full data
        } else {
          data[, selected_columns, drop = FALSE] # Filter columns
        }
      })

      output$top_features_download <- downloadHandler(
        filename = function() {
          paste("top_features_", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
          data <- topFeatures() # Fetch data

          selected_columns <- input$top_features_columns # Get user-selected columns

          if (is.null(selected_columns)) {
            write.csv(data, file, row.names = FALSE) # Write complete data
          } else {
            write.csv(data[, selected_columns, drop = FALSE], file, row.names = FALSE) # Filter columns in downloaded file
          }
        }
      )
    })
  }

  # Return the Shiny application object
  shiny::shinyApp(ui = ui, server = server)
}
