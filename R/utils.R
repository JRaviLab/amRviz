#' @importFrom shinydashboard box tabBox
#' @importFrom dplyr filter mutate select group_by summarize ungroup arrange
#' @importFrom dplyr pull distinct left_join bind_rows slice_max slice_head
#' @importFrom dplyr collect desc n c_across all_of rowwise case_when join_by
#' @importFrom dplyr count summarise where
#' @importFrom ggplot2 ggplot aes geom_col geom_line geom_point geom_boxplot
#' @importFrom ggplot2 theme_bw theme_minimal theme element_text labs
#' @importFrom ggplot2 scale_fill_brewer scale_color_brewer coord_cartesian
#' @importFrom ggplot2 position_jitterdodge ggtitle
#' @importFrom here here
#' @importFrom DBI dbConnect dbGetQuery dbWriteTable
#' @importFrom duckdb duckdb
#' @importFrom tibble tibble column_to_rownames
#' @importFrom tidyr pivot_wider separate
#' @importFrom purrr map_dfr
#' @importFrom stringr str_glue str_extract str_split_i str_replace_all
#' @importFrom stringr str_to_lower str_trunc str_flatten str_match str_detect
#' @importFrom stringr str_remove str_to_sentence str_extract_all
#' @importFrom readr read_tsv read_csv
#' @importFrom rlang sym .data
#' @importFrom arrow read_parquet write_parquet
#' @importFrom plotly plot_ly layout colorbar renderPlotly plotlyOutput
#' @importFrom ComplexHeatmap Heatmap draw
#' @importFrom circlize colorRamp2
#' @importFrom grid gpar unit
#' @importFrom DT datatable
#' @importFrom glue glue
#' @importFrom stats setNames
#' @importFrom utils packageVersion
NULL

# Species code regex (note: Esp. includes escaped period)
SPECIES_PATTERN <- "(Efa|Sau|Kpn|Aba|Pae|Esp\\.?)"

# normalize_species helper: make "Esp." and "Esp" equivalent by removing
# a single trailing dot for comparisons (preserves NA).
normalize_species <- function(x) {
  x_chr <- as.character(x)
  x_chr[is.na(x_chr)] <- NA_character_
  stringr::str_replace_all(x_chr, "\\.$", "")
}

# getHoldoutsDrugChoices: derive drug/class choices for the holdouts tab
# perf_data: combined performance tibble from loadMLResults()
# bug: optional 3-letter species code to filter to
getHoldoutsDrugChoices <- function(perf_data, bug = NULL) {
  if (is.null(perf_data) || !is.data.frame(perf_data) || !nrow(perf_data)) {
    return(character(0))
  }

  # Filter to stratified models only (country or year, not baseline)
  df <- perf_data %>%
    dplyr::filter(!is.na(.data$strat_label) & nzchar(.data$strat_label))

  if ("species" %in% names(df)) {
    df <- df %>% dplyr::mutate(species = normalize_species(.data$species))
  }

  if (!is.null(bug) && "species" %in% names(df)) {
    bug_norm <- normalize_species(bug)
    df <- dplyr::filter(df, .data$species %in% bug_norm)
  }

  if (!nrow(df)) {
    return(character(0))
  }

  x <- as.character(df$drug_or_class)
  x[!is.na(x) & nzchar(trimws(x))] |>
    unique() |>
    sort()
}


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

# Bug choices
bug_choices <- c(
  "Staphylococcus aureus" = "Sau",
  "Klebsiella pneumoniae" = "Kpn",
  "Acinetobacter baumannii" = "Aba",
  "Pseudomonas aeruginosa" = "Pae",
  "Enterobacter spp." = "Esp.",
  "Enterococcus faecium" = "Efa"
)


getTopFeatures <- function(top_features_dir, pattern) {
  # List files in the directory
  file_list <- list.files(top_features_dir, pattern, full.names = TRUE)

  # Exclude files with double underscores
  valid_files <- file_list[!grepl("country__", basename(file_list)) &
    !(grepl("trimethoprim", basename(file_list)) &
      grepl("sulfonamides", basename(file_list)) &
      grepl("derivatives", basename(file_list)))]

  # load drugs names;
  df <- lapply(X = valid_files, FUN = function(x) {
    transformed_x <- stringr::str_replace_all(x, pattern = "United_Kingdom", "UK")

    df <- readr::read_tsv(file.path(x), show_col_types = F) |>
      dplyr::mutate(species = stringr::str_extract(transformed_x, pattern = SPECIES_PATTERN)) |>
      dplyr::mutate(drug_type = stringr::str_extract(transformed_x, pattern = "drug_class|drug")) |>
      dplyr::mutate(feature_type = stringr::str_extract(transformed_x, pattern = "proteins|domains|genes|struct")) |>
      dplyr::mutate(data_type = stringr::str_extract(transformed_x, pattern = "counts|binary|struct")) |>
      dplyr::mutate(configurations = ifelse(!is.na(stringr::str_extract(transformed_x, pattern = "PCA|shuffled")), stringr::str_extract(transformed_x, pattern = "PCA|shuffled"), "full")) |>
      dplyr::mutate(configurations = stringr::str_to_lower(configurations)) |>
      dplyr::mutate(source_file = basename(x))
    if (stringr::str_detect(transformed_x, pattern = "country_models")) {
      df <- dplyr::mutate(
        df,
        drug_or_class_name = stringr::str_match(transformed_x, pattern = stringr::str_glue("{drug_type}_country_\\s*(.*?)\\s*_{feature_type}"))[, 2] |>
          str_remove("_[^_]+$")
      ) |>
        dplyr::mutate(trained_country = stringr::str_match(transformed_x, pattern = stringr::str_glue("{drug_or_class_name}_\\s*(.*?)\\s*_{feature_type}"))[, 2]) |>
        dplyr::mutate(tested_country = trained_country)
    } else if (stringr::str_detect(transformed_x, pattern = "cross_models")) {
      df <- dplyr::mutate(
        df,
        drug_or_class_name = stringr::str_match(transformed_x, pattern = stringr::str_glue("{drug_type}_country_\\s*(.*?)\\s*_{pattern}"))[, 2] |>
          stringr::str_remove("(_[^_]+_[^_]+)$")
      ) |>
        dplyr::mutate(trained_country = stringr::str_match(transformed_x, pattern = stringr::str_glue("{drug_or_class_name}_\\s*(.*?)\\s*_{pattern}"))[, 2]) |>
        tidyr::separate(
          col = trained_country,
          into = c("trained_country", "tested_country"),
          sep = "_",
          remove = FALSE
        )
    } else if (stringr::str_detect(transformed_x, pattern = "year_models")) {
      df <- dplyr::mutate(
        df,
        drug_or_class_name = stringr::str_match(transformed_x, pattern = stringr::str_glue("{drug_type}_year_\\s*(.*?)\\s*_(?=\\d)"))[, 2],
        trained_year = stringr::str_extract(transformed_x, pattern = stringr::str_glue("(\\d+)_(\\d+)")),
        tested_year = trained_year
      )
    } else if (stringr::str_detect(transformed_x, pattern = "cross_year")) {
      if (stringr::str_detect(x, "(\\d+)_(\\d+)_(\\d+)_(\\d+)")) {
        years <- unlist(stringr::str_extract_all(transformed_x, pattern = stringr::str_glue("(\\d+)_(\\d+)")))
        second_year <- as.integer(stringr::str_extract(years[2], pattern = "(\\d+)$"))
        df <- dplyr::mutate(
          df,
          drug_or_class_name = stringr::str_match(transformed_x, pattern = stringr::str_glue("{drug_type}_year_\\s*(.*?)\\s*_(?=\\d)"))[, 2],
          trained_year = stringr::str_extract(transformed_x, pattern = stringr::str_glue("(\\d+)_(\\d+)")),
          trained_year = years[1],
          tested_year = stringr::str_glue("{second_year-4}_{second_year}")
        )
      }
    } else {
      df <- dplyr::mutate(
        df,
        drug_or_class_name = stringr::str_match(transformed_x, pattern = stringr::str_glue("{drug_type}_\\s*(.*?)\\s*_{feature_type}"))[, 2]
      )
    }
  }) |>
    bind_rows()
  # convert the feature_type with struct to genes
  df <- df |>
    dplyr::mutate(feature_type = ifelse(feature_type == "struct", "genes", feature_type))

  ## write to parquet file
  if (stringr::str_detect(top_features_dir, pattern = "country_models") | stringr::str_detect(top_features_dir, pattern = "cross_models")) {
    results_folder <- ifelse(stringr::str_detect(top_features_dir, pattern = "country_models"), "country_models", "cross_models")
    data_file <- ifelse(pattern == "top_features", "top_features.parquet", "performance_metrics.parquet")
    arrow::write_parquet(df, sink = here::here("shinyapp", "data", results_folder, data_file))
  } else if (stringr::str_detect(top_features_dir, pattern = "cross_year") | stringr::str_detect(top_features_dir, pattern = "year_models")) {
    results_folder <- ifelse(stringr::str_detect(top_features_dir, pattern = "cross_year"), "cross_year", "year_models")
    data_file <- ifelse(pattern == "top_features", "top_features.parquet", "performance_metrics.parquet")
    arrow::write_parquet(df, sink = here::here("shinyapp", "data", results_folder, data_file))
  } else {
    results_folder <- ifelse(pattern == "top_features", "top_features", "performance_metrics")
    data_file <- ifelse(pattern == "top_features", "top_features.parquet", "performance_metrics.parquet")
    arrow::write_parquet(df, sink = here::here("shinyapp", "data", results_folder, data_file))
  }
  return(df)
}

combinePerformanceMetrics <- function() {
  # Define the paths to the individual performance_metrics.parquet files
  country_metrics_path <- here::here("shinyapp", "data", "country_models", "performance_metrics.parquet")
  cross_country_metrics_path <- here::here("shinyapp", "data", "cross_models", "performance_metrics.parquet")
  year_metrics_path <- here::here("shinyapp", "data", "year_models", "performance_metrics.parquet")
  cross_year_metrics_path <- here::here("shinyapp", "data", "cross_year", "performance_metrics.parquet")

  # Read the individual files and filter out rows with "country__" files
  country_metrics <- read_parquet(country_metrics_path) %>%
    dplyr::filter(!grepl("country__", source_file))

  cross_country_metrics <- read_parquet(cross_country_metrics_path) %>%
    dplyr::filter(!grepl("country__", source_file))

  year_metrics <- read_parquet(year_metrics_path) %>%
    dplyr::filter(!grepl("country__", source_file))

  cross_year_metrics <- read_parquet(cross_year_metrics_path) %>%
    dplyr::filter(!grepl("country__", source_file))

  # Combine all the data
  combined_metrics <- bind_rows(
    country_metrics,
    cross_country_metrics,
    year_metrics,
    cross_year_metrics
  )

  # Write the combined data to the main performance_metrics.parquet file
  combined_metrics_path <- here::here("shinyapp", "data", "performance_metrics", "performance_metrics.parquet")
  arrow::write_parquet(combined_metrics, combined_metrics_path)

  message("Combined performance_metrics.parquet file has been updated.")
}

getAllFilteredData <- function(amr_db_path, filtered_table_path) {
  amr_eskape_files <- list.files(amr_db_path)
  amr_filtered_tbl_db <- DBI::dbConnect(
    duckdb::duckdb(),
    dbdir = filtered_table_path
  )
  for (i in seq_along(amr_eskape_files)) {
    eskape_path <- file.path(amr_db_path, amr_eskape_files[i])
    eskape_spp_name <- stringr::str_split_i(amr_eskape_files[i], "\\.", 1)

    eskape_db_conn <- DBI::dbConnect(
      duckdb::duckdb(),
      dbdir = eskape_path
    )
    # filter the filtered data table
    eskape_filterered_data_table <- dbGetQuery(
      eskape_db_conn,
      stringr::str_glue(
        "SELECT * FROM {eskape_spp_name}_filtered"
      )
    )
    # store this table in a new connections;
    dbWriteTable(
      amr_filtered_tbl_db,
      name = stringr::str_glue("{eskape_spp_name}_filtered"),
      value = eskape_filterered_data_table,
      overwrite = TRUE
    )
  }
}


loadAMRMetadata <- function(spp_name) {
  message(stringr::str_glue("Duckdb version {packageVersion('duckdb')}"))
  message(stringr::str_glue("DBI version {packageVersion('DBI')}"))
  cwd <- getwd()
  db_fp <- file.path(cwd, "data", "amr_filtered_tbls.db")
  cleaned_iso_src_fp <- file.path(cwd, "data", "llm_gpt_amr_isolation_source_cleaned.csv")
  cleaned_iso_src_df <- readr::read_csv(
    cleaned_iso_src_fp,
    show_col_types = FALSE
  )
  message(stringr::str_glue("Looking for DB at: {db_fp}"))
  if (!file.exists(db_fp)) {
    stop(glue::glue("Database file not found at {db_fp}"))
  }
  amr_filtered_tbl_db <- DBI::dbConnect(
    duckdb::duckdb(),
    dbdir = db_fp
  )

  dplyr::tbl(
    amr_filtered_tbl_db,
    stringr::str_glue("{spp_name}_filtered")
  ) %>%
    dplyr::filter(!is.na(genome_drug.antibiotic)) %>%
    dplyr::filter(genome_drug.resistant_phenotype %in% c("Resistant", "Susceptible")) %>%
    dplyr::collect() %>%
    dplyr::left_join(
      cleaned_iso_src_df,
      by = "genome.isolation_source"
    ) %>%
    dplyr::mutate(genome.isolation_source = group)
}

loadDrugClassMap <- function() {
  cwd <- getwd()
  # drug_class_map_fp <- file.path(cwd, "data", "drug_class_map.tsv")
  drug_class_map_fp <- system.file("extdata", "drug_class_map.tsv", package = "amRshiny")
  message(stringr::str_glue("loadDrugClassMap(): Looking for TSV at: {drug_class_map_fp}"))
  drug_class_map_df <- readr::read_tsv(
    here(drug_class_map_fp),
    show_col_types = FALSE
  ) %>%
    dplyr::select(drug.antibiotic_name, drug_class) %>%
    dplyr::distinct()
  return(drug_class_map_df)
}

## NOTE: loadMLResults() and loadTopFeat() with results_root/species_dirs args
## are defined below (lines ~1391+). These single-arg stubs are superseded.

getAmrDrugs <- function(amr_drug, amr_drug_class) {
  # Load the drug class mapping data
  amr_drug_class_df <- loadDrugClassMap()
  # Initialize amr_drugs as NULL
  amr_drugs <- NULL

  # If a specific drug is selected, use it
  if (length(amr_drug) > 0) {
    return(amr_drug)
  }

  # If a drug class is selected, get all drugs for this class
  if (length(amr_drug_class) > 0) {
    amr_drugs <- amr_drug_class_df %>%
      dplyr::filter(drug_class %in% amr_drug_class) %>%
      pull(drug.antibiotic_name) %>%
      unique()
  }

  # Return the vector of selected drugs
  return(amr_drugs)
}

quickStatBox <- function(title, value, icon_name, bg_color, text_color = "white") {
  div(
    style = glue::glue("
      background-color: {bg_color};
      color: {text_color};
      padding: 10px 14px;
      margin: 4px 3px;
      border-radius: 10px;
      box-shadow: 0 2px 6px rgba(0,0,0,0.12);
      display: flex;
      justify-content: space-between;
      align-items: center;
      min-height: 55px;
      font-family: sans-serif;
    "),
    div(
      div(
        style = "font-size: 18px; font-weight: bold; line-height: 1.2;",
        value
      ),
      div(
        style = "font-size: 11px; margin-top: 3px; opacity: 0.9;",
        title
      )
    ),
    div(
      style = "font-size: 24px; opacity: 0.4;",
      icon(icon_name)
    )
  )
}

makeQuickStats <- function(data) {
  df <- data

  total_records <- nrow(df)
  total_genomes <- length(unique(df$genome_drug.genome_id))
  n_drugs <- length(unique(df$genome_drug.antibiotic))
  n_drug_classes <- length(unique(df$drug_class[!is.na(df$drug_class)]))

  pheno_tbl <- df[df$genome_drug.resistant_phenotype %in% c("Resistant", "Susceptible"), ] |>
    (\(x) table(x$genome_drug.resistant_phenotype))() |>
    as.data.frame()
  n_resistant <- if ("Resistant" %in% pheno_tbl$Var1) pheno_tbl[pheno_tbl$Var1 == "Resistant", "Freq"] else 0
  n_susceptible <- if ("Susceptible" %in% pheno_tbl$Var1) pheno_tbl[pheno_tbl$Var1 == "Susceptible", "Freq"] else 0

  top_drugs <- df %>%
    dplyr::count(.data$genome_drug.antibiotic, name = "n") %>%
    dplyr::arrange(dplyr::desc(.data$n)) %>%
    dplyr::slice_head(n = 5) %>%
    dplyr::pull(.data$genome_drug.antibiotic)

  top_classes <- df %>%
    dplyr::filter(!is.na(.data$drug_class)) %>%
    dplyr::count(.data$drug_class, name = "n") %>%
    dplyr::arrange(dplyr::desc(.data$n)) %>%
    dplyr::slice_head(n = 5) %>%
    dplyr::pull(.data$drug_class)

  top_countries <- df %>%
    dplyr::filter(
      !is.na(.data$genome.isolation_country),
      .data$genome.isolation_country != "",
      !grepl("^NA$", .data$genome.isolation_country)
    ) %>%
    dplyr::count(.data$genome.isolation_country, name = "n") %>%
    dplyr::arrange(dplyr::desc(.data$n)) %>%
    dplyr::slice_head(n = 5) %>%
    dplyr::pull(.data$genome.isolation_country)

  spp_label <- unique(df$species)
  tagList(
    tags$p(
      style = "font-weight: bold; padding: 4px 4px 4px; font-family: sans-serif; font-size: 14px;",
      paste("Data summary for", stringr::str_to_sentence(spp_label[1]))
    ),
    fluidRow(
      column(4, quickStatBox("Isolate-drug records", total_records, "database", "#3c5a6f")),
      column(4, quickStatBox("Unique genomes", total_genomes, "dna", "#7aab6e")),
      column(4, quickStatBox("Drugs tested", n_drugs, "pills", "#9b7fba"))
    ),
    fluidRow(
      column(4, quickStatBox("Resistant isolates", n_resistant, "virus", "#d4872a")),
      column(4, quickStatBox("Susceptible isolates", n_susceptible, "shield-halved", "#5b8db8")),
      column(4, quickStatBox("Drug classes", n_drug_classes, "layer-group", "#8b6b7a"))
    ),
    fluidRow(
      column(4, quickStatBox(
        "Top 5 drugs",
        tags$span(
          style = "font-size:11px; line-height:1.5;",
          HTML(paste(top_drugs, collapse = "<br/>"))
        ),
        "star", "#4e9a9a"
      )),
      column(4, quickStatBox(
        "Top 5 drug classes",
        tags$span(
          style = "font-size:11px; line-height:1.5;",
          HTML(paste(top_classes, collapse = "<br/>"))
        ),
        "list", "#b5b5b5"
      )),
      column(4, quickStatBox(
        "Top 5 countries",
        tags$span(
          style = "font-size:11px; line-height:1.5;",
          HTML(paste(top_countries, collapse = "<br/>"))
        ),
        "globe", "#d4735e"
      ))
    )
  )
}

makeDatAvailabilityPlot <- function(data) {
  # Palette: resistant=warm amber, susceptible=muted slate blue, intermediate=peach
  phenotype_pal <- c(
    R = "#d4872a", Resistant = "#d4872a", resistant = "#d4872a",
    S = "#5b8db8", Susceptible = "#5b8db8", susceptible = "#5b8db8",
    I = "#e6ab80", Intermediate = "#e6ab80", intermediate = "#e6ab80"
  )
  data <- data %>%
    dplyr::group_by(.data$genome_drug.antibiotic, .data$genome_drug.resistant_phenotype) %>%
    dplyr::count() %>%
    dplyr::ungroup()
  g <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data$genome_drug.antibiotic,
      y = .data$n,
      fill = .data$genome_drug.resistant_phenotype,
      text = paste0(
        "Drug: ", .data$genome_drug.antibiotic,
        "<br>Phenotype: ", .data$genome_drug.resistant_phenotype,
        "<br>Count: ", .data$n
      )
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = phenotype_pal, na.value = "gray70") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(angle = 90, hjust = 1, size = 10),
      legend.title = ggplot2::element_text(size = 12),
      legend.text  = ggplot2::element_text(size = 10),
      axis.text    = ggplot2::element_text(size = 10),
      axis.title   = ggplot2::element_text(size = 10)
    ) +
    ggplot2::labs(x = "Drug", y = "No. of isolates", fill = "AMR phenotype")
  plotly::ggplotly(g, tooltip = "text")
}

makeGeoChloroPlot <- function(data) {
  if (!requireNamespace("countrycode", quietly = TRUE)) {
    stop("Package 'countrycode' is required for this function. Install it with install.packages('countrycode').")
  }
  data$iso3 <- countrycode::countrycode(data$genome.isolation_country, origin = "country.name", destination = "iso3c")
  plot_ly(
    data = data,
    type = "choropleth",
    locations = ~iso3,
    z = ~count,
    text = ~ paste(iso3, "<br>Genomes:", count),
    hoverinfo = "text",
    colorscale = list(c(0, "#fdf3e3"), c(0.5, "#d4a050"), c(1, "#d4872a")),
    marker = list(line = list(width = 0.5, color = "white"))
  ) %>%
    colorbar(title = list(text = "No. of isolates", font = list(size = 14))) %>%
    layout(
      # title = list(text = "Geographic distribution", font = list(size = 14)),
      geo = list(
        projection = list(type = "natural earth"),
        showcoastlines = TRUE,
        showland = TRUE,
        landcolor = "lightgrey"
      ),
      font = list(size = 14)
    )
}

makeTimeSeriesAMRPlot <- function(data, amr_drug) {
  whole_data_title <- stringr::str_glue(
    "AMR resistance trend"
  )

  title_amr_drug <- stringr::str_glue(
    "AMR resistance trend for {amr_drug}"
  )

  title_amr_drug <- ifelse(
    amr_drug == "all",
    whole_data_title,
    title_amr_drug
  )

  phenotype_pal <- c(
    R = "#d4872a", Resistant = "#d4872a", resistant = "#d4872a",
    S = "#5b8db8", Susceptible = "#5b8db8", susceptible = "#5b8db8",
    I = "#e6ab80", Intermediate = "#e6ab80", intermediate = "#e6ab80"
  )
  g <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data$genome.collection_year,
      y = .data$n,
      colour = .data$genome_drug.resistant_phenotype,
      group = .data$genome_drug.resistant_phenotype,
      text = paste0(
        "Year: ", .data$genome.collection_year,
        "<br>Phenotype: ", .data$genome_drug.resistant_phenotype,
        "<br>Count: ", .data$n
      )
    )
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::scale_color_manual(values = phenotype_pal, na.value = "gray70") +
    ggplot2::labs(
      title  = title_amr_drug,
      x      = "Year",
      y      = "No. of isolates",
      colour = "AMR phenotype"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.title = ggplot2::element_text(size = 12),
      legend.text  = ggplot2::element_text(size = 10),
      axis.text    = ggplot2::element_text(size = 10),
      axis.title   = ggplot2::element_text(size = 10)
    )
  plotly::ggplotly(g, tooltip = "text")
}

makeHostIsolatePlot <- function(data) {
  data <- data %>%
    dplyr::mutate(
      genome.host_common_name = stringr::str_to_lower(genome.host_common_name),
      genome.isolation_source = stringr::str_to_lower(genome.isolation_source)
    )

  host_df <- data %>%
    dplyr::group_by(genome_drug.antibiotic, genome.host_common_name) %>%
    dplyr::summarize(n = n()) %>%
    dplyr::ungroup()

  meta_pal <- c(
    "#4a6b8a", "#87ceeb", "#e6ab80", "#8b6b7a",
    "#4e9a9a", "#c4a35a", "#9b7fba", "#5b8db8",
    "#d4735e", "#d4872a", "#b5b5b5"
  )
  g <- ggplot2::ggplot(
    host_df,
    ggplot2::aes(
      x = .data$genome_drug.antibiotic,
      y = .data$n,
      fill = .data$genome.host_common_name,
      text = paste0(
        "Drug: ", .data$genome_drug.antibiotic,
        "<br>Host: ", .data$genome.host_common_name,
        "<br>Count: ", .data$n
      )
    )
  ) +
    ggplot2::geom_col(position = "stack") +
    ggplot2::scale_fill_manual(values = meta_pal, na.value = "gray70") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 10),
      axis.text.y = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 10),
      axis.title  = ggplot2::element_text(size = 10)
    ) +
    ggplot2::labs(x = "Drug", y = "No. of isolates", fill = "Host")
  plotly::ggplotly(g, tooltip = "text")
}

makeIsolationSourcesPlot <- function(data) {
  isolation_source <- data %>%
    dplyr::mutate(genome.isolation_source = stringr::str_to_lower(genome.isolation_source)) |>
    dplyr::group_by(genome_drug.antibiotic, genome.isolation_source) %>%
    dplyr::summarize(n = n()) %>%
    dplyr::ungroup()

  top_isolate_source <- isolation_source |>
    dplyr::group_by(genome.isolation_source) |>
    dplyr::summarize(n = sum(n)) |>
    dplyr::filter(genome.isolation_source != "") |>
    dplyr::slice_max(order_by = n, n = 10) |>
    dplyr::pull(genome.isolation_source)

  # Filter the isolation source data to include only the top five sources
  isolation_source <- isolation_source |>
    dplyr::mutate(
      genome.isolation_source_ = ifelse(
        genome.isolation_source %in% top_isolate_source,
        genome.isolation_source,
        "Other"
      )
    )
  isolation_source <- isolation_source %>%
    mutate(genome.isolation_source_ = stringr::str_trunc(genome.isolation_source_, width = 20))

  meta_pal <- c(
    "#4a6b8a", "#87ceeb", "#e6ab80", "#8b6b7a",
    "#4e9a9a", "#c4a35a", "#9b7fba", "#5b8db8",
    "#d4735e", "#d4872a", "#b5b5b5"
  )
  isolation_source_plot <- ggplot(
    isolation_source,
    aes(
      x = genome_drug.antibiotic,
      y = n,
      fill = genome.isolation_source_,
      text = paste0(
        "Drug: ", genome_drug.antibiotic,
        "<br>Source: ", genome.isolation_source_,
        "<br>Count: ", n
      )
    )
  ) +
    geom_col(position = "stack") +
    ggplot2::scale_fill_manual(values = meta_pal, na.value = "gray70") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      legend.text = element_text(size = 10),
      axis.title = element_text(size = 10)
    ) +
    labs(
      x = "",
      y = "No. of isolates",
      fill = "Isolation source"
    )
  plotly::ggplotly(isolation_source_plot, tooltip = "text")
}

# makeModelPerformancePlot: plot baseline ML model performance metrics.
# data: pre-loaded performance tibble from loadMLResults() / queryData()
# Columns used from amRml parquet schema:
#   species, feature_type (scale), feature_subtype (data type),
#   drug_or_class (drug/class abbrev), drug_label ("drug"/"drug_class"), nmcc/bal_acc/f1
makeModelPerformancePlot <- function(
  data, bug, model_scale, data_type, metrics,
  amr_drug_class, amr_drug
) {
  empty_plot <- function(msg) {
    plotly::plot_ly() %>%
      plotly::layout(title = list(text = msg, x = 0))
  }

  if (is.null(data) || !is.data.frame(data) || !nrow(data)) {
    return(empty_plot("No data available"))
  }

  # Filter to baseline models (strat_label is NA = no country/year stratification)
  df <- data %>%
    dplyr::filter(normalize_species(.data$species) %in% normalize_species(bug)) %>%
    dplyr::filter(.data$feature_type %in% model_scale) %>%
    dplyr::filter(.data$feature_subtype %in% data_type) %>%
    dplyr::filter(is.na(.data$strat_label) | !nzchar(.data$strat_label))

  # Filter by drug class or drug
  has_class <- !is.null(amr_drug_class) && length(amr_drug_class) > 0 &&
    !identical(amr_drug_class, "all")
  has_drug <- !is.null(amr_drug) && length(amr_drug) > 0 &&
    nzchar(amr_drug[1]) && !identical(amr_drug, "all")

  if (has_class && has_drug) {
    df <- df %>%
      dplyr::filter(
        (.data$drug_label == "drug_class" & .data$drug_or_class %in% amr_drug_class) |
          (.data$drug_label == "drug" & .data$drug_or_class %in% amr_drug)
      )
  } else if (has_class) {
    df <- df %>%
      dplyr::filter(
        .data$drug_label == "drug_class" & .data$drug_or_class %in% amr_drug_class
      )
  } else if (has_drug) {
    df <- df %>%
      dplyr::filter(
        .data$drug_label == "drug" & .data$drug_or_class %in% amr_drug
      )
  }

  if (!nrow(df)) {
    return(empty_plot("No data for current selection"))
  }

  # Use species_label for display if available
  if ("species_label" %in% names(df)) {
    df <- df %>%
      dplyr::mutate(species_display = gsub("_", " ", .data$species_label))
  } else {
    df <- df %>%
      dplyr::mutate(species_display = normalize_species(.data$species))
  }

  df <- df %>%
    dplyr::mutate(
      feature_type = factor(
        .data$feature_type,
        levels = c("domains", "genes", "proteins", "struct")
      )
    )

  # Match Performance Overview palette
  scale_colors <- c(
    domains  = "#66a61e",
    genes    = "#e6ab80",
    proteins = "#87ceeb",
    struct   = "#a52a2a"
  )

  # Build native plotly grouped boxplot — ggplotly doesn't handle dodge reliably
  p <- plotly::plot_ly()

  for (sc in levels(df$feature_type)) {
    sub <- df[df$feature_type == sc, ]
    if (!nrow(sub)) next
    col <- if (sc %in% names(scale_colors)) scale_colors[[sc]] else "#999999"
    y_vals <- sub[[metrics]]
    tooltip_vals <- paste0(
      "Species: ", sub$species_display,
      "<br>Scale: ", sc,
      "<br>Drug/class: ", sub$drug_or_class,
      "<br>", metrics, ": ", round(y_vals, 3),
      "<br>Encoding: ", sub$feature_subtype
    )
    p <- p %>%
      plotly::add_trace(
        type       = "box",
        x          = sub$species_display,
        y          = y_vals,
        name       = sc,
        text       = tooltip_vals,
        hoverinfo  = "text",
        boxpoints  = "all",
        jitter     = 0.3,
        pointpos   = 0,
        marker     = list(color = col, opacity = 0.7, size = 5),
        line       = list(color = col),
        fillcolor  = paste0(col, "66")
      )
  }

  p %>% plotly::layout(
    boxmode = "group",
    title   = list(text = "Performance metrics", x = 0, font = list(size = 13, color = "#333333", family = "Arial, sans-serif")),
    xaxis   = list(title = "Species"),
    yaxis   = list(title = metrics, range = c(0, 1.05)),
    legend  = list(title = list(text = "Scale")),
    margin  = list(t = 50)
  )
}

# makeFeatureImportancePlot: heatmap of top features across bugs or drugs.
# data: pre-loaded top-features tibble from loadTopFeat() / topFeatures()
# amRml column mapping (new → expected here):
#   drug_or_class  → drug/class abbreviation identifier
#   feature_subtype → data encoding (binary/counts)
#   feature_type    → molecular scale for baseline (genes/domains/proteins/struct)
#   strat_label     → NA for baseline models
# Annotation join is attempted from results_root/Annotated/ or extdata/Annotated/;
# if no annotated files found, Variable name is used directly as feature label.
makeFeatureImportancePlot <- function(
  data, bug, amr_drug, model_scale, data_type_,
  top_n_features, feature_importance_tabset,
  annotated_dir = NULL
) {
  if (is.null(data) || !is.data.frame(data) || !nrow(data)) {
    return(NULL)
  }

  # Map UI scale label to singular for annotated file matching
  scale <- dplyr::case_when(
    model_scale == "proteins" ~ "protein",
    model_scale == "domains" ~ "domain",
    model_scale == "genes" ~ "gene",
    TRUE ~ model_scale
  )

  top_n_features <- ifelse(
    identical(top_n_features, "all"), "all",
    as.numeric(top_n_features)
  )

  bug_norm <- normalize_species(bug)

  # Filter top features using new amRml column names:
  #   feature_type  = molecular scale (genes/domains/proteins/struct) for baseline
  #   feature_subtype = encoding (binary/counts)
  #   drug_or_class = drug/class abbreviation
  #   strat_label   = NA for baseline (no country/year stratification)
  top_features_df <- data %>%
    dplyr::mutate(species = normalize_species(.data$species)) %>%
    dplyr::filter(.data$species %in% bug_norm) %>%
    dplyr::filter(.data$drug_or_class %in% amr_drug) %>%
    dplyr::filter(.data$feature_type %in% c(model_scale, "struct")) %>%
    dplyr::filter(.data$feature_subtype %in% data_type_) %>%
    dplyr::filter(is.na(.data$strat_label) | !nzchar(.data$strat_label))

  if (!nrow(top_features_df)) {
    return(NULL)
  }

  # group column: across_bug groups by species; across_drug by drug_or_class
  group_column <- dplyr::case_when(
    feature_importance_tabset == "across_drug" ~ "drug_or_class",
    feature_importance_tabset == "across_bug" ~ "species"
  )

  # Attempt to load annotated files for COG name lookup
  ann_dirs <- c(
    annotated_dir,
    system.file("extdata", "Annotated", package = "amRshiny")
  )
  ann_dirs <- ann_dirs[!is.null(ann_dirs) & nzchar(ann_dirs) & dir.exists(ann_dirs)]

  annotated_files <- character(0)
  if (length(ann_dirs) > 0) {
    annotated_files <- unlist(lapply(ann_dirs, function(d) {
      fls <- list.files(
        d,
        pattern = stringr::str_flatten(bug_norm, collapse = "|"),
        full.names = TRUE
      )
      Filter(function(x) grepl(scale, x, fixed = TRUE), fls)
    }))
  }

  has_annotation <- length(annotated_files) > 0

  if (has_annotation) {
    annotated_table <- purrr::map_dfr(annotated_files, function(x) {
      sp <- stringr::str_extract(basename(x), SPECIES_PATTERN)
      arrow::read_parquet(x) |>
        dplyr::mutate(species = normalize_species(sp))
    })

    join_by_expr <- switch(paste(group_column, scale, sep = "_"),
      "species_protein" = join_by(Variable == "proteinID", "species" == "species"),
      "species_domain" = join_by(Variable == "PfamID", "species" == "species"),
      "species_gene" = join_by(Variable == "Gene", "species" == "species"),
      "drug_or_class_protein" = join_by(Variable == "proteinID"),
      "drug_or_class_domain" = join_by(Variable == "PfamID"),
      "drug_or_class_gene" = join_by(Variable == "Gene"),
      NULL
    )

    if (scale == "protein") {
      annotated_table <- annotated_table |>
        dplyr::mutate(proteinID = stringr::str_replace(.data$proteinID, "\\|", "."))
    }
    if (scale == "domain") {
      top_features_df <- top_features_df |>
        dplyr::mutate(Variable = stringr::str_split_i(.data$Variable, "_", 1))
    }

    if (!is.null(join_by_expr)) {
      top_features_df <- tryCatch(
        top_features_df %>%
          dplyr::inner_join(annotated_table, by = join_by_expr) %>%
          dplyr::filter(!is.na(.data$COG_name)),
        error = function(e) {
          message("Annotation join failed: ", conditionMessage(e))
          top_features_df
        }
      )
    }
  }

  # If no annotation join produced COG_name, fall back to Variable
  if (!"COG_name" %in% names(top_features_df)) {
    top_features_df <- top_features_df %>%
      dplyr::mutate(COG_name = .data$Variable)
  }
  if (!nrow(top_features_df)) {
    return(NULL)
  }

  # Aggregate: max importance per group x COG
  top_features_df <- top_features_df %>%
    dplyr::group_by(!!rlang::sym(group_column), .data$COG_name) %>%
    dplyr::summarize(Importance = max(.data$Importance, na.rm = TRUE), .groups = "drop")

  # Min-max normalise within each group
  top_features_df <- top_features_df %>%
    dplyr::group_by(!!rlang::sym(group_column)) %>%
    dplyr::mutate(
      Importance = (.data$Importance - min(.data$Importance, na.rm = TRUE)) /
        (max(.data$Importance, na.rm = TRUE) - min(.data$Importance, na.rm = TRUE))
    ) %>%
    dplyr::ungroup()

  # Slice top N features per group
  if (!identical(top_n_features, "all")) {
    top_features_df <- top_features_df %>%
      dplyr::group_by(!!rlang::sym(group_column)) %>%
      dplyr::slice_max(order_by = .data$Importance, n = top_n_features) %>%
      dplyr::ungroup()
  }

  if (!nrow(top_features_df)) {
    return(NULL)
  }

  # Build wide matrix
  if (feature_importance_tabset == "across_bug") {
    vi_wider <- top_features_df %>%
      dplyr::select(.data$COG_name, .data$Importance, .data$species) %>%
      dplyr::distinct() %>%
      tidyr::pivot_wider(names_from = "species", values_from = "Importance")

    group_cols <- setdiff(colnames(vi_wider), "COG_name")
    if (!length(group_cols)) {
      return(NULL)
    }

    vi_wider <- vi_wider %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        .n = sum(!is.na(c_across(all_of(group_cols)))),
        .mx = max(c_across(all_of(group_cols)), na.rm = TRUE)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(dplyr::desc(.data$.n), dplyr::desc(.data$.mx)) %>%
      dplyr::select(-.data$.n, -.data$.mx)

    vi_mat <- vi_wider %>%
      tibble::column_to_rownames("COG_name") %>%
      as.matrix()

    eskape_order <- c("Efa", "Sau", "Kpn", "Aba", "Pae", "Esp")
    col_ord <- intersect(eskape_order, colnames(vi_mat))
    if (length(col_ord)) vi_mat <- vi_mat[, col_ord, drop = FALSE]
  }

  if (feature_importance_tabset == "across_drug") {
    top_features_df <- top_features_df %>%
      dplyr::mutate(drug_or_class = stringr::str_trim(as.character(.data$drug_or_class)))

    vi_wider <- top_features_df %>%
      dplyr::select(.data$COG_name, .data$Importance, .data$drug_or_class) %>%
      dplyr::group_by(.data$COG_name, .data$drug_or_class) %>%
      dplyr::summarise(Importance = max(.data$Importance, na.rm = TRUE), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = "drug_or_class", values_from = "Importance")

    group_cols <- setdiff(colnames(vi_wider), "COG_name")
    if (!length(group_cols)) {
      return(NULL)
    }

    vi_wider <- vi_wider %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        .n = sum(!is.na(c_across(all_of(group_cols)))),
        .mx = max(c_across(all_of(group_cols)), na.rm = TRUE)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(dplyr::desc(.data$.n), dplyr::desc(.data$.mx)) %>%
      dplyr::select(-.data$.n, -.data$.mx)

    vi_mat <- vi_wider %>%
      tibble::column_to_rownames("COG_name") %>%
      as.matrix()
  }

  if (!exists("vi_mat") || !length(vi_mat)) {
    return(NULL)
  }

  # Build plotly heatmap (features = rows, groups = columns)
  plotly::plot_ly(
    x = colnames(vi_mat),
    y = rownames(vi_mat),
    z = vi_mat,
    type = "heatmap",
    colorscale = list(c(0, "#c8e8e8"), c(1, "#4e9a9a")),
    colorbar = list(title = "Importance"),
    hovertemplate = paste0(
      "<b>Feature:</b> %{y}<br>",
      "<b>Group:</b> %{x}<br>",
      "<b>Importance:</b> %{z:.3f}<extra></extra>"
    )
  ) %>%
    plotly::layout(
      xaxis = list(title = "", tickangle = -45, side = "top"),
      yaxis = list(title = "", autorange = "reversed"),
      margin = list(l = 200, b = 20, t = 80)
    )
}


# makeCrossModelFeatureImportancePlot: heatmap of top features for holdout models.
# top_data: pre-loaded top-features tibble (country or year stratified rows).
# amRml column mapping:
#   drug_or_class → drug/class abbreviation
#   strat_label   → "country" or "year"
#   strat_value   → trained-on country/year
makeCrossModelFeatureImportancePlot <- function(
  top_data, bug, drug, cross_model, top_n_features,
  annotated_dir = NULL
) {
  if (is.null(top_data) || !is.data.frame(top_data) || !nrow(top_data)) {
    return(NULL)
  }

  top_n_features <- ifelse(
    identical(top_n_features, "all"), "all",
    as.numeric(top_n_features)
  )

  strat <- if (cross_model == "country") "country" else "year"
  strat_col <- "strat_value" # trained-on value in new schema

  features_df <- top_data %>%
    dplyr::filter(.data$species %in% bug) %>%
    dplyr::filter(.data$drug_or_class %in% drug) %>%
    dplyr::filter(.data$strat_label == strat) %>%
    dplyr::filter(!.data$cross_test)

  if (!nrow(features_df)) {
    return(NULL)
  }

  if (!identical(top_n_features, "all")) {
    features_df <- features_df %>%
      dplyr::group_by(!!rlang::sym(strat_col)) %>%
      dplyr::slice_max(order_by = .data$Importance, n = top_n_features) %>%
      dplyr::ungroup()
  }

  vi_wider <- features_df %>%
    dplyr::select(.data$Variable, .data$Importance, !!rlang::sym(strat_col)) %>%
    tidyr::pivot_wider(
      names_from = strat_col,
      values_from = "Importance",
      values_fn = mean
    )

  if (!nrow(vi_wider)) {
    return(NULL)
  }

  vi_mat <- vi_wider %>%
    tibble::column_to_rownames("Variable") %>%
    as.matrix()

  # column-wise min-max normalisation
  vi_mat <- apply(vi_mat, 2, function(x) {
    rng <- range(x, na.rm = TRUE)
    if (diff(rng) == 0) {
      return(x)
    }
    (x - rng[1]) / diff(rng)
  })

  max_val <- max(vi_mat, na.rm = TRUE)
  min_val <- min(vi_mat, na.rm = TRUE)
  if (min_val == max_val) {
    min_val <- max(0, min_val - 0.1)
    max_val <- min(1, max_val + 0.1)
  }

  plotly::plot_ly(
    x = colnames(vi_mat),
    y = rownames(vi_mat),
    z = vi_mat,
    type = "heatmap",
    colorscale = list(c(0, "#c8e8e8"), c(1, "#4e9a9a")),
    colorbar = list(title = "Importance"),
    hovertemplate = paste0(
      "<b>Feature:</b> %{y}<br>",
      "<b>Group:</b> %{x}<br>",
      "<b>Importance:</b> %{z:.3f}<extra></extra>"
    )
  ) %>%
    plotly::layout(
      xaxis = list(title = "", tickangle = -45, side = "top"),
      yaxis = list(title = "", autorange = "reversed"),
      margin = list(l = 200, b = 20, t = 80)
    )
}

# makeCrossModelPerformancePlot: heatmap of balanced accuracy for holdout models.
# perf_data: pre-loaded performance tibble (country or year stratified rows).
# amRml column mapping:
#   drug_or_class   → drug/class abbreviation
#   strat_label     → "country" or "year"
#   strat_value     → trained-on country/year
#   strat_value_test → tested-on country/year (NA for self-evaluation)
makeCrossModelPerformancePlot <- function(perf_data, bug, drug, cross_model) {
  if (is.null(perf_data) || !is.data.frame(perf_data) || !nrow(perf_data)) {
    return(NULL)
  }

  strat <- if (cross_model == "country") "country" else "year"

  df <- perf_data %>%
    dplyr::filter(normalize_species(.data$species) %in% normalize_species(bug)) %>%
    dplyr::filter(.data$drug_or_class %in% drug) %>%
    dplyr::filter(.data$strat_label == strat)

  if (!nrow(df)) {
    return(NULL)
  }

  # For self-evaluation rows (strat_value_test is NA), set tested = trained
  df <- df %>%
    dplyr::mutate(
      strat_value_test = dplyr::if_else(
        is.na(.data$strat_value_test),
        .data$strat_value,
        .data$strat_value_test
      )
    )

  # Aggregate bal_acc (mean across scales/encodings for same train/test pair)
  models_performance <- df %>%
    dplyr::group_by(.data$strat_value, .data$strat_value_test) %>%
    dplyr::summarise(bal_acc = mean(.data$bal_acc, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(
      names_from = "strat_value",
      values_from = "bal_acc"
    ) %>%
    tibble::column_to_rownames("strat_value_test") %>%
    as.matrix()

  if (!length(models_performance)) {
    return(NULL)
  }

  min_val <- min(models_performance, na.rm = TRUE)
  max_val <- max(models_performance, na.rm = TRUE)
  if (min_val == max_val) {
    min_val <- max(0, min_val - 0.1)
    max_val <- min(1, max_val + 0.1)
  }

  row_title <- if (cross_model == "country") "Tested Country" else "Tested Year"
  col_title <- if (cross_model == "country") "Trained Country" else "Trained Year"

  plotly::plot_ly(
    x = colnames(models_performance),
    y = rownames(models_performance),
    z = models_performance,
    type = "heatmap",
    colorscale = list(c(0, "lightblue"), c(1, "darkblue")),
    colorbar = list(title = "Balanced\nAccuracy"),
    hovertemplate = paste0(
      "<b>", col_title, ":</b> %{x}<br>",
      "<b>", row_title, ":</b> %{y}<br>",
      "<b>Bal. Accuracy:</b> %{z:.3f}<extra></extra>"
    )
  ) %>%
    plotly::layout(
      xaxis = list(title = col_title, side = "top"),
      yaxis = list(title = row_title, autorange = "reversed"),
      margin = list(l = 120, b = 20, t = 80)
    )
}

makeFeatureImportTable <- function(feature_import_table) {
  # Early exit for empty/zero-column data
  if (is.null(feature_import_table) || ncol(feature_import_table) == 0) {
    return(DT::datatable(tibble::tibble(), options = list(dom = "t"), rownames = FALSE))
  }

  # Preferred display order (only those that exist will be used)
  cols_priority <- c(
    "species", "drug_or_class",
    "COG_name", "COG_description", "ARG_name", "ARG_description",
    "Gene", "Annotation", "accession",
    "feature_type", "feature_subtype", "Importance"
  )
  existing <- intersect(cols_priority, names(feature_import_table))

  feature_import_table <- feature_import_table |>
    dplyr::mutate(
      dplyr::across(where(is.numeric), ~ formatC(.x, format = "e", digits = 3))
    ) |>
    dplyr::mutate(dplyr::across(
      dplyr::any_of("ARG_name"),
      ~ stringr::str_replace_all(.x, "non-ARG", "-")
    )) |>
    # Reorder for display: preferred columns first, then everything else
    dplyr::select(dplyr::all_of(existing), dplyr::everything())

  tbl <- feature_import_table
  if ("accession" %in% names(tbl)) {
    tbl <- tbl |>
      dplyr::mutate(accession = stringr::str_split_i(.data$accession, "\\.", 1)) |>
      dplyr::mutate(accession = dplyr::case_when(
        !is.na(.data$accession) ~ stringr::str_glue(
          "<a href='https://www.ncbi.nlm.nih.gov/genome/annotation_prok/evidence/{accession}'",
          " target='_blank' style='color:#1a73e8; text-decoration: underline;'>{accession}</a>"
        ),
        TRUE ~ .data$accession
      ))
  }
  if ("COG_name" %in% names(tbl)) {
    tbl <- tbl |>
      dplyr::mutate(COG_name = stringr::str_glue(
        "<a href='https://www.ncbi.nlm.nih.gov/research/cog/cog/{COG_name}'",
        " target='_blank' style='color:#1a73e8; text-decoration: underline;'>{COG_name}</a>"
      ))
  }

  DT::datatable(
    tbl,
    options = list(scrollX = TRUE, autoWidth = FALSE, orderClasses = TRUE),
    class = "display nowrap stripe",
    rownames = FALSE,
    width = "100%",
    escape = FALSE
  )
}

# helpers to discover species folders and load to combine results generated from amRml
.normalize_results_root <- function(results_root) {
  if (is.null(results_root) || length(results_root) != 1 || is.na(results_root) || !nzchar(results_root)) {
    return(NULL)
  }
  normalizePath(results_root, winslash = "/", mustWork = FALSE)
}

.read_parquet_safe <- function(fp, verbose = TRUE) {
  if (is.null(fp) || !file.exists(fp)) {
    if (isTRUE(verbose)) message("File not found: ", fp)
    return(tibble::tibble())
  }
  tryCatch(
    arrow::read_parquet(fp),
    error = function(e) {
      if (isTRUE(verbose)) message("Failed to read parquet: ", fp, " (", conditionMessage(e), ")")
      tibble::tibble()
    }
  )
}

# Discover species subdirectories under results_root that contain amRml output.
# Files live inside per-species subdirectories: {root}/{SpeciesDir}/{code}_ML_perf.parquet
# Returns a named character vector: names = directory basename (display label),
# values = full path to the species subdirectory.
listAmRmlSpeciesFolders <- function(results_root, verbose = TRUE) {
  rr <- .normalize_results_root(results_root)
  if (is.null(rr) || !dir.exists(rr)) {
    return(character(0))
  }

  subdirs <- list.dirs(rr, full.names = TRUE, recursive = FALSE)
  if (!length(subdirs)) {
    return(character(0))
  }

  # Keep subdirs that contain at least one baseline *_ML_perf.parquet file
  has_perf <- vapply(subdirs, function(d) {
    fps <- list.files(d, pattern = "_ML_perf\\.parquet$", full.names = FALSE)
    any(!grepl("_(country|year|MDR|cross)_ML_perf\\.parquet$", fps))
  }, logical(1))

  ok <- subdirs[has_perf]
  if (!length(ok)) {
    return(character(0))
  }
  setNames(ok, basename(ok))
}

# Load all performance parquets (baseline + country + year + cross) from a species directory.
# species_dir: full path to the species subdirectory (e.g. "/results/Shigella_flexneri")
# Attaches species_label = basename(species_dir) so the full name is available for display.
.load_one_species_perf <- function(species_dir, verbose = TRUE) {
  fps <- list.files(species_dir, pattern = "_ML_perf\\.parquet$", full.names = TRUE)
  fps <- fps[!grepl("_MDR_ML_perf\\.parquet$", fps)]
  if (!length(fps)) {
    return(tibble::tibble())
  }
  df <- dplyr::bind_rows(lapply(fps, .read_parquet_safe, verbose = verbose))
  if (nrow(df)) df$species_label <- basename(species_dir)
  df
}

# Load all top-feature parquets (baseline + country + year) from a species directory.
.load_one_species_top <- function(species_dir, verbose = TRUE) {
  fps <- list.files(species_dir, pattern = "_ML_top_features\\.parquet$", full.names = TRUE)
  fps <- fps[!grepl("_MDR_ML_top_features\\.parquet$", fps)]
  if (!length(fps)) {
    return(tibble::tibble())
  }
  df <- dplyr::bind_rows(lapply(fps, .read_parquet_safe, verbose = verbose))
  if (nrow(df)) df$species_label <- basename(species_dir)
  df
}

# Public loaders used by app.R (multi-species aware).
# species_dirs: character vector of full paths to species subdirectories selected by the user.
loadMLResults <- function(results_root = NULL, species_dirs = NULL, verbose = TRUE) {
  rr <- .normalize_results_root(results_root)

  # User mode: results_root + species selected → load from selected subdirectories
  if (!is.null(rr) && !is.null(species_dirs) && length(species_dirs) > 0) {
    dfs <- lapply(species_dirs, .load_one_species_perf, verbose = verbose)
    return(dplyr::bind_rows(dfs))
  }

  # User mode: results_root provided but nothing selected yet
  if (!is.null(rr) && is.null(species_dirs)) {
    return(tibble::tibble())
  }

  # Demo fallback: scan extdata subdirectories recursively for *_ML_perf.parquet
  extdata <- system.file("extdata", package = "amRshiny")
  if (!nzchar(extdata)) {
    return(tibble::tibble())
  }
  subdirs <- list.dirs(extdata, full.names = TRUE, recursive = FALSE)
  if (isTRUE(verbose)) message("loadMLResults(): using packaged demo parquets")
  dplyr::bind_rows(lapply(subdirs, .load_one_species_perf, verbose = verbose))
}

loadTopFeat <- function(results_root = NULL, species_dirs = NULL, verbose = TRUE) {
  rr <- .normalize_results_root(results_root)

  if (!is.null(rr) && !is.null(species_dirs) && length(species_dirs) > 0) {
    dfs <- lapply(species_dirs, .load_one_species_top, verbose = verbose)
    return(dplyr::bind_rows(dfs))
  }

  if (!is.null(rr) && is.null(species_dirs)) {
    return(tibble::tibble())
  }

  # Demo fallback: scan extdata subdirectories for *_ML_top_features.parquet
  extdata <- system.file("extdata", package = "amRshiny")
  if (!nzchar(extdata)) {
    return(tibble::tibble())
  }
  subdirs <- list.dirs(extdata, full.names = TRUE, recursive = FALSE)
  if (isTRUE(verbose)) message("loadTopFeat(): using packaged demo parquets")
  dplyr::bind_rows(lapply(subdirs, .load_one_species_top, verbose = verbose))
}

# Return path to a species metadata parquet.
# Searches inside the species subdirectory under results_root (user mode) or
# extdata (demo mode). The metadata file follows the pattern {code}_metadata.parquet
# and lives alongside the other parquets for that species.
get_metadata_path <- function(species_code, results_root = NULL) {
  fname <- paste0(species_code, "_metadata.parquet")

  # User mode: search species subdirectories under results_root
  if (!is.null(results_root) && nzchar(results_root)) {
    subdirs <- list.dirs(results_root, full.names = TRUE, recursive = FALSE)
    for (d in subdirs) {
      fp <- file.path(d, fname)
      if (file.exists(fp)) {
        return(fp)
      }
    }
  }

  # Demo mode: search species subdirectories under extdata
  extdata <- system.file("extdata", package = "amRshiny")
  if (nzchar(extdata)) {
    subdirs <- list.dirs(extdata, full.names = TRUE, recursive = FALSE)
    for (d in subdirs) {
      fp <- file.path(d, fname)
      if (file.exists(fp)) {
        return(fp)
      }
    }
  }
  NULL
}

# ── nMCC overview plots ──────────────────────────────────────────────────────

# Shared helper: filter perf data to baseline (non-stratified, non-cross) rows
# and attach a species_display column from species_label (or species code).
.prep_nmcc_data <- function(data) {
  if (is.null(data) || !is.data.frame(data) || !nrow(data)) {
    return(NULL)
  }
  df <- data %>%
    dplyr::filter(is.na(.data$strat_label) | !nzchar(.data$strat_label)) %>%
    dplyr::filter(!.data$cross_test, !is.na(.data$nmcc))
  if (!nrow(df)) {
    return(NULL)
  }
  if ("species_label" %in% names(df)) {
    df <- df %>%
      dplyr::mutate(species_display = gsub("_", " ", .data$species_label))
  } else {
    df <- df %>%
      dplyr::mutate(species_display = .data$species)
  }
  df
}

# makeNmccStripPlot: Panel A — nMCC distribution by species, faceted by
# molecular scale.  Uses ggplot2 box + jitter converted to plotly.
# data: performance tibble from queryData() / loadMLResults()
makeNmccStripPlot <- function(data, selected_drug_class = NULL, selected_drug = NULL) {
  df <- .prep_nmcc_data(data)
  if (is.null(df)) {
    return(plotly::plot_ly() %>% plotly::layout(title = "No data available"))
  }

  scale_colors <- c(
    "Domain"  = "#66a61e",
    "Gene"    = "#e6ab80",
    "Protein" = "#87ceeb",
    "Struct"  = "#a52a2a"
  )
  scale_order <- c("domains", "genes", "proteins", "struct")
  scale_labels <- c(domains = "Domain", genes = "Gene", proteins = "Protein", struct = "Struct")

  df <- df %>%
    dplyr::filter(.data$feature_type %in% scale_order) %>%
    dplyr::mutate(
      scale_label = factor(
        scale_labels[.data$feature_type],
        levels = unname(scale_labels)
      )
    )

  if (!nrow(df)) {
    return(plotly::plot_ly() %>% plotly::layout(title = "No data for selected filters"))
  }

  # ── Highlight selected drug or drug class ──
  # Priority: specific drug > drug class. Points matching the selection get
  # full alpha/larger size; all others are dimmed when a selection is active.
  # When drug_class == "all", suppress both highlights — no specific class
  # is selected so the accompanying drug selection is not meaningful.
  class_active <- !is.null(selected_drug_class) &&
    nzchar(selected_drug_class) &&
    selected_drug_class != "all"
  use_drug <- class_active && !is.null(selected_drug) && nzchar(selected_drug)
  use_class <- class_active && !use_drug

  df <- df %>%
    dplyr::mutate(
      highlighted = dplyr::case_when(
        use_drug ~ .data$drug_or_class == selected_drug & .data$drug_label == "drug",
        use_class ~ .data$drug_or_class == selected_drug_class & .data$drug_label == "drug_class",
        TRUE ~ FALSE
      )
    )

  any_highlighted <- use_drug || use_class

  df <- df %>%
    dplyr::mutate(
      pt_alpha = dplyr::if_else(.data$highlighted & any_highlighted, 0.9,
        dplyr::if_else(any_highlighted, 0.15, 0.7)
      ),
      pt_size = dplyr::if_else(.data$highlighted & any_highlighted, 3.5,
        dplyr::if_else(any_highlighted, 1.5, 2.5)
      )
    )

  g <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data$species_display,
      y = .data$nmcc,
      color = .data$scale_label,
      fill = .data$scale_label,
      alpha = .data$pt_alpha,
      size = .data$pt_size,
      text = paste0(
        "Drug/class: ", .data$drug_or_class,
        "\nnMCC: ", round(.data$nmcc, 3),
        "\nEncoding: ", .data$feature_subtype
      )
    )
  ) +
    ggplot2::geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.5, linewidth = 0.4) +
    ggplot2::geom_jitter(width = 0.15) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50", linewidth = 0.4) +
    ggplot2::facet_grid(scale_label ~ ., switch = "y") +
    ggplot2::scale_color_manual(values = scale_colors) +
    ggplot2::scale_fill_manual(values = scale_colors) +
    ggplot2::scale_alpha_identity() +
    ggplot2::scale_size_identity() +
    ggplot2::coord_cartesian(ylim = c(0.35, 1.05)) +
    ggplot2::labs(x = NULL, y = "nMCC") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position   = "none",
      strip.placement   = "outside",
      strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1),
      panel.spacing     = ggplot2::unit(0.3, "lines")
    )

  plotly::ggplotly(g, tooltip = "text") %>%
    plotly::layout(
      title = list(
        text = "nMCC per species and molecular scale", x = 0,
        font = list(size = 13, color = "#333333", family = "Arial, sans-serif")
      ),
      margin = list(t = 50)
    )
}

# makeNmccHeatmap: Panel B — compound heatmap showing median nMCC by drug class
# across three groupings: species (grayscale), molecular scale (categorical
# colors), and data encoding (categorical colors).
# Only drug_class rows are used (drug_label == "drug_class").
makeNmccHeatmap <- function(data, selected_drug_class = NULL) {
  df <- .prep_nmcc_data(data)
  if (is.null(df)) {
    return(plotly::plot_ly() %>% plotly::layout(title = "No data available"))
  }
  df <- df %>% dplyr::filter(.data$drug_label == "drug_class")
  if (!nrow(df)) {
    return(plotly::plot_ly() %>% plotly::layout(title = "No drug-class data available"))
  }

  # Drug class order: most-represented first (bottom of y-axis = first row)
  drug_order <- df %>%
    dplyr::count(.data$drug_or_class) %>%
    dplyr::arrange(dplyr::desc(.data$n)) %>%
    dplyr::pull(.data$drug_or_class)
  drug_order_rev <- rev(drug_order) # ggplot y reversed so top row = first

  df <- df %>%
    dplyr::mutate(drug_or_class = factor(.data$drug_or_class, levels = drug_order_rev))

  # ── Section 1: species × drug_class (grayscale median nMCC) ──
  spp_summ <- df %>%
    dplyr::group_by(.data$drug_or_class, .data$species_display) %>%
    dplyr::summarise(med_nmcc = median(.data$nmcc, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(drug_or_class = factor(.data$drug_or_class, levels = drug_order_rev))

  g1 <- ggplot2::ggplot(
    spp_summ,
    ggplot2::aes(
      x = .data$species_display,
      y = .data$drug_or_class,
      fill = .data$med_nmcc,
      text = paste0(
        "Species: ", .data$species_display,
        "\nDrug class: ", .data$drug_or_class,
        "\nMedian nMCC: ", round(.data$med_nmcc, 3)
      )
    )
  ) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient(
      low = "#f7fbff", high = "#08306b",
      limits = c(0.5, 1.0), name = "nMCC", na.value = "white"
    ) +
    ggplot2::labs(x = NULL, y = "Drug class", title = "Species") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 12, color = "#333333", family = "sans"),
      axis.text.x = ggplot2::element_text(angle = 40, hjust = 1),
      legend.position = "bottom"
    )

  # ── Section 2: molecular scale × drug_class (nMCC gradient per scale color) ──
  scale_colors <- c(
    Domain  = "#66a61e",
    Gene    = "#e6ab80",
    Protein = "#87ceeb",
    Struct  = "#a52a2a"
  )
  scale_order <- c("domains", "genes", "proteins", "struct")
  scale_labels <- c(domains = "Domain", genes = "Gene", proteins = "Protein", struct = "Struct")

  sc_all <- df %>%
    dplyr::filter(.data$feature_type %in% scale_order) %>%
    dplyr::group_by(.data$drug_or_class, .data$feature_type) %>%
    dplyr::summarise(med_nmcc = median(.data$nmcc, na.rm = TRUE), .groups = "drop")

  # ## Top-ranked only (commented out — revert by uncommenting and swapping sc_summ below)
  # sc_max <- sc_all %>%
  #   dplyr::group_by(.data$drug_or_class) %>%
  #   dplyr::summarise(max_nmcc = max(.data$med_nmcc, na.rm = TRUE), .groups = "drop")
  # sc_summ_top <- sc_all %>%
  #   dplyr::left_join(sc_max, by = "drug_or_class") %>%
  #   dplyr::mutate(
  #     scale_label   = factor(scale_labels[.data$feature_type], levels = unname(scale_labels)),
  #     drug_or_class = factor(.data$drug_or_class, levels = drug_order_rev),
  #     fill_label    = dplyr::if_else(
  #       .data$med_nmcc >= .data$max_nmcc - 1e-9,
  #       as.character(.data$scale_label), NA_character_
  #     )
  #   )

  # All tiles, alpha encodes nMCC intensity (0.5 → light, 1.0 → full color)
  sc_summ <- sc_all %>%
    dplyr::mutate(
      scale_label   = factor(scale_labels[.data$feature_type], levels = unname(scale_labels)),
      drug_or_class = factor(.data$drug_or_class, levels = drug_order_rev),
      alpha_val     = (.data$med_nmcc - 0.5) / 0.5 # normalize 0.5–1.0 → 0–1
    )

  g2 <- ggplot2::ggplot(
    sc_summ,
    ggplot2::aes(
      x = .data$scale_label,
      y = .data$drug_or_class,
      fill = .data$scale_label,
      alpha = .data$alpha_val,
      text = paste0(
        "Scale: ", .data$scale_label,
        "\nDrug class: ", .data$drug_or_class,
        "\nMedian nMCC: ", round(.data$med_nmcc, 3)
      )
    )
  ) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_manual(values = scale_colors, guide = "none") +
    ggplot2::scale_alpha_continuous(range = c(0.15, 1.0), guide = "none") +
    ggplot2::labs(x = NULL, y = NULL, title = "Molecular scale") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title  = ggplot2::element_text(hjust = 0.5, size = 12, color = "#333333", family = "sans"),
      axis.text.x = ggplot2::element_text(angle = 40, hjust = 1),
      axis.text.y = ggplot2::element_blank()
    )

  # ── Section 3: data encoding × drug_class (nMCC gradient per encoding color) ──
  subtype_colors <- c(Binary = "#6495ED", Counts = "#BA55D3")

  st_all <- df %>%
    dplyr::group_by(.data$drug_or_class, .data$feature_subtype) %>%
    dplyr::summarise(med_nmcc = median(.data$nmcc, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(
      subtype_label = dplyr::case_when(
        .data$feature_subtype == "binary" ~ "Binary",
        .data$feature_subtype == "counts" ~ "Counts",
        TRUE ~ .data$feature_subtype
      )
    )

  # ## Top-ranked only (commented out — revert by uncommenting and swapping st_summ below)
  # st_max <- st_all %>%
  #   dplyr::group_by(.data$drug_or_class) %>%
  #   dplyr::summarise(max_nmcc = max(.data$med_nmcc, na.rm = TRUE), .groups = "drop")
  # st_summ_top <- st_all %>%
  #   dplyr::left_join(st_max, by = "drug_or_class") %>%
  #   dplyr::mutate(
  #     drug_or_class = factor(.data$drug_or_class, levels = drug_order_rev),
  #     fill_label    = dplyr::if_else(
  #       .data$med_nmcc >= .data$max_nmcc - 1e-9,
  #       .data$subtype_label, NA_character_
  #     )
  #   )

  # All tiles, alpha encodes nMCC intensity
  st_summ <- st_all %>%
    dplyr::mutate(
      drug_or_class = factor(.data$drug_or_class, levels = drug_order_rev),
      alpha_val     = (.data$med_nmcc - 0.5) / 0.5
    )

  g3 <- ggplot2::ggplot(
    st_summ,
    ggplot2::aes(
      x = .data$subtype_label,
      y = .data$drug_or_class,
      fill = .data$subtype_label,
      alpha = .data$alpha_val,
      text = paste0(
        "Encoding: ", .data$subtype_label,
        "\nDrug class: ", .data$drug_or_class,
        "\nMedian nMCC: ", round(.data$med_nmcc, 3)
      )
    )
  ) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_manual(values = subtype_colors, guide = "none") +
    ggplot2::scale_alpha_continuous(range = c(0.15, 1.0), guide = "none") +
    ggplot2::labs(x = NULL, y = NULL, title = "Data type") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title  = ggplot2::element_text(hjust = 0.5, size = 12, color = "#333333", family = "sans"),
      axis.text.x = ggplot2::element_text(angle = 40, hjust = 1),
      axis.text.y = ggplot2::element_blank()
    )

  # ── Highlight selected drug class row ──
  # Compute the 1-indexed factor level position of the selected drug class.
  # drug_order_rev[1] = bottom row (factor level 1 in ggplot2 y-axis).
  highlight_y <- if (
    !is.null(selected_drug_class) &&
      nzchar(selected_drug_class) &&
      selected_drug_class != "all" &&
      selected_drug_class %in% drug_order_rev
  ) {
    which(drug_order_rev == selected_drug_class)
  } else {
    NULL
  }
  if (!is.null(highlight_y)) {
    highlight_rect <- ggplot2::annotate(
      "rect",
      xmin = -Inf, xmax = Inf,
      ymin = highlight_y - 0.5, ymax = highlight_y + 0.5,
      fill = NA, color = "#FFD700", linewidth = 1.5
    )
    g1 <- g1 + highlight_rect
    g2 <- g2 + highlight_rect
    g3 <- g3 + highlight_rect
  }

  p1 <- plotly::ggplotly(g1, tooltip = "text")
  p2 <- plotly::ggplotly(g2, tooltip = "text")
  p3 <- plotly::ggplotly(g3, tooltip = "text")

  plotly::subplot(p1, p2, p3,
    nrows  = 1,
    shareY = TRUE,
    widths = c(0.45, 0.35, 0.20),
    margin = 0.03
  ) %>%
    plotly::layout(
      title = list(
        text = "nMCC by drug class, species, molecular scale, and data type", x = 0,
        font = list(size = 13, color = "#333333", family = "Arial, sans-serif")
      ),
      legend = list(orientation = "h", y = -0.15)
    )
}
