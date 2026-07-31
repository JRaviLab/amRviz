#' AMR Shiny dashboard application
#'
#' Creates and returns a Shiny application for exploring antimicrobial
#' resistance data and machine learning model results.
#'
#' @param results_root File path to the root directory containing amRml model
#'        output results. If `NULL` (default), the application will attempt to
#'        use example data bundled with the package, where available.
#' @param amrdata_root File path to the root directory containing amRdata
#'        annotation parquets (e.g. `cluster_feature.parquet`,
#'        `gene_names.parquet`). If `NULL` (default), `~/amRdata/data` is used
#'        when present; otherwise annotation-based features are disabled.
#'
#' @return A Shiny application object
#' @export
#' @import shiny
#' @importFrom utils head write.csv
#' @importFrom shinyjs useShinyjs
#' @examples
#' if (interactive()) {
#'   app <- launchAMRDashboard()
#'   shiny::runApp(app)
#' }
launchAMRDashboard <- function(results_root = NULL,
                               amrdata_root = NULL) {
  # Default amrdata_root: ~/amRdata/data if it exists
  if (is.null(amrdata_root)) {
    default_amrdata <- file.path(path.expand("~"), "amRdata", "data")
    if (dir.exists(default_amrdata)) amrdata_root <- default_amrdata
  }
  # UI
  ui <- tagList(
    shinyjs::useShinyjs(),
    tags$head(includeCSS(system.file("app/www/style.css", package = "amRviz"))),
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

                      /* App header (title row, above nav tabs) */
                      .amr-app-header {
                        background-color: #2b2b2b;
                        padding: 14px 24px;
                        display: flex;
                        justify-content: space-between;
                        align-items: center;
                        border-bottom: 1px solid #3a3a3a;
                      }
                      .amr-app-header-title {
                        color: #ffffff;
                        font-weight: bold;
                        font-size: 26px;
                        letter-spacing: 1.5px;
                        font-family: sans-serif;
                      }
                      .amr-app-header-logo img {
                        height: 60px;
                        border-radius: 50%;
                      }

                      /* Navbar (tabs row) - dark theme matching app header */
                      .navbar {
                        background-color: #1a1a1a !important;
                        border: none !important;
                        margin-bottom: 0 !important;
                        min-height: 44px !important;
                        border-radius: 0 !important;
                      }
                      .navbar-default {
                        background-color: #1a1a1a !important;
                        border-color: #1a1a1a !important;
                      }
                      /* hide the empty brand spacer */
                      .navbar-brand {
                        display: none !important;
                      }
                      .navbar-nav > li > a {
                        color: #cccccc !important;
                        font-size: 14px !important;
                        padding: 13px 18px !important;
                        transition: background-color 0.2s, color 0.2s;
                      }
                      .navbar-nav > li > a:hover {
                        color: #ffffff !important;
                        background-color: #333333 !important;
                      }
                      .navbar-nav > li.active > a,
                      .navbar-nav > li.active > a:focus,
                      .navbar-nav > li.active > a:hover {
                        color: #ffffff !important;
                        background-color: #2b2b2b !important;
                        border-bottom: 3px solid #5b9bd5 !important;
                      }
                      .home-tab-icon {
                        font-size: 16px;
                      }
                    "))
    ),
    # App header: title row separate from the nav tabs
    tags$div(
      class = "amr-app-header",
      tags$div(
        class = "amr-app-header-title",
        "amRviz"
      ),
      tags$div(
        class = "amr-app-header-logo",
        tags$img(src = "www/logo.png", onerror = "this.style.display='none'")
      )
    ),
    navbarPage(
      id = "tabselected",
      selected = "home",
      title = "",
      # 2. Home icon tab (right of AMR dashboard)
      tabPanel(
        title = icon("home", class = "home-tab-icon"),
        value = "home",
        fluidPage(
          style = "max-width: 960px; margin: 0 auto; padding: 24px 16px;",

          # Overview & Features
          h3("Overview and Features",
            style = "font-weight: bold; margin-bottom: 12px;"
          ),
          tags$p("amRviz allows users to:"),
          tags$ul(
            style = "line-height: 2; margin-bottom: 16px;",
            tags$li("Explore AMR model performance across species, drugs, and drug classes."),
            tags$li("Compare predictive performance across molecular feature scales (gene, protein, domain, structure)."),
            tags$li("Identify key genomic features driving resistance predictions."),
            tags$li("Analyze model generalization across geography and time."),
            tags$li("Visualize and filter isolate metadata and model results interactively.")
          ),
          tags$p(
            style = "margin-bottom: 28px;",
            tags$em("amRviz is interactive, modular, and scalable for exploring AMR data and machine learning outputs.")
          ),

          # Workflow figure
          div(
            style = "text-align: center; margin-bottom: 32px;",
            tags$img(
              src = "www/amr_overview.png",
              style = "max-width: 70%; border-radius: 6px; box-shadow: 0 2px 10px rgba(0,0,0,0.12);",
              onerror = "this.style.display='none'"
            )
          ),
          tags$hr(),
          h2("amR: an R package suite to predict antimicrobial resistance in bacterial pathogens"),
          tags$p(
            tags$strong("Authors: "),
            "Evan P Brenner^, Abhirupa Ghosh^, Emily A Boyer, Charmie K Vang, Ethan P Wolfe, Alexander P McKim, Raymond L Lesiyon, David Mayer, Janani Ravi*."
          ),
          tags$p(
            "Department of Biomedical Informatics, Center for Health Artificial Intelligence, University of Colorado Anschutz, Aurora, CO 80045.",
            tags$span(style = "font-style: italic", " ^Co-primary authors contributing equally. *Corresponding author: janani.ravi@cuanschutz.edu")
          ),
          tags$p(
            tags$strong("Keywords: "),
            "antimicrobial resistance, machine learning, bacterial genomics, pangenomics, interpretable models, drug resistance prediction, multiscale features, R package"
          ),
          br(),
          h4("Abstract"),
          tags$strong("Motivation: "),
          tags$p("Identifying antimicrobial resistance (AMR) in bacterial pathogens is critical for diagnostics and treatment, but resistance is a complex trait arising from diverse mechanisms spanning multiple molecular scales. Existing computational approaches often function as black boxes and rarely explore cross-species or multi-drug patterns. We developed amR, an integrated R package suite providing a complete framework from bacterial genome sequences to interpretable AMR predictions, enabling identification of resistance mechanisms across species and drugs."),
          tags$strong("Results: "),
          tags$p("The amR suite contains three modular packages. amRdata interfaces with BV-BRC to download and process bacterial genomes with paired antimicrobial susceptibility testing data, constructs pangenomes, and extracts features at gene/protein cluster, protein domain, and structural variant scales, along with annotating protein clusters with Clusters of Orthologous Genes and ResFinder AMR-associated features. Data are stored in memory-efficient Parquet and DuckDB formats. amRml trains interpretable machine learning models per species-drug combination, calculates feature importance and comprehensive performance metrics, and provides rich ground for mechanism discovery. amRviz provides an interactive Shiny dashboard to explore metadata distributions, model performance across species and drugs, visualize top predictive AMR features, and analyze cross-model patterns (including across geographic/temporal strata). The suite has been applied to ESKAPE pathogens, achieving balanced accuracies above 0.80. With thousands of genomes, multi-scale features, and interpretable models, amR provides an accessible, comprehensive framework for AMR research."),
          tags$strong("Availability and implementation: "),
          tags$p(
            "amR is developed in R. We use Dockerized software for data curation and feature extraction, perform modeling with the tidymodels framework, and visualize results using Shiny. The suite is available at: ",
            tags$a(href = "https://github.com/JRaviLab/amR", "https://github.com/JRaviLab/amR"), "."
          ),
          tags$p(
            tags$strong("Contact: "),
            tags$a(href = "mailto:janani.ravi@cuanschutz.edu", "janani.ravi@cuanschutz.edu")
          ),
          tags$hr(),
        )
      ),
      # other tabs
      metadataUI(),
      modelPerfUI(),
      featureImportanceUI(),
      crossModelComparisonUI(),
      networkUI(),
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
    # Shared reactives, data loaders, and cross-tab UI-update observers used
    # by every tab module below. `core` is a named list of reactive accessors.
    core <- setupServerCore(input, output, session, results_root)
    queryData <- core$queryData
    topFeatures <- core$topFeatures
    available_species <- core$available_species
    available_metadata_species <- core$available_metadata_species
    filtered_top_features <- core$filtered_top_features
    loadDrugClassMapRec <- core$loadDrugClassMapRec
    bug_norm_input <- core$bug_norm_input

    # Metadata tab
    serverMetadata(input, output, session, results_root)

    # Model Performance + Performance overview tabs (includes the bug/class/drug
    # cascading selector observers).
    serverModelPerf(input, output, session, core, results_root)

    # Bug/Drug feature comparison tab
    serverFeatureImportance(input, output, session, core, results_root, amrdata_root)

    # Model-holdouts tab
    serverCrossModel(input, output, session, core)


    # Drug-feature network tab
    serverNetwork(input, output, session, core, results_root)

    # Query data + top-features tables and CSV downloads
    serverQuery(input, output, session, core)
  }

  # Register inst/app/www/ so images and CSS are served from the package
  shiny::addResourcePath(
    "www",
    system.file("app/www", package = "amRviz")
  )

  # Return the Shiny application object
  shiny::shinyApp(ui = ui, server = server)
}
