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
  drug_class_map_fp <- system.file("extdata", "drug_class_map.tsv", package = "amRviz")
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

quickStatBox <- function(title, value, unit = "", bg_color = "#f0f0f0", text_color = "black") {
  div(
    style = glue::glue("
      background-color: {bg_color};
      color: {text_color};
      padding: 20px;
      margin: 5px auto;
      border-radius: 10px;
      width:95%;
      box-shadow: 0 2px 4px rgba(0,0,0,0.1);
      text-align: center;
      font-family: sans-serif;
    "),
    div(style = "font-size: 16px;", title),
    div(
      style = "font-size: 12px; margin-top: 5px;",
      paste0(value, if (unit != "") paste0(" ", unit))
    )
  )
}

makeQuickStats <- function(data) { # , drug_class_df, spp_name, amr_drugs) {
  data_with_drug_class <- data
  # Sample stat calculations
  total_genomes <- nrow(data_with_drug_class)
  total_uniques_genomes <- length(unique(data_with_drug_class$genome_drug.genome_id))
  n_amr_drugs <- length(unique(data_with_drug_class$genome_drug.antibiotic))
  n_amr_drug_class <- length(unique(data_with_drug_class$drug_class[!is.na(data_with_drug_class$drug_class)]))
  total_strains <- data_with_drug_class[
    data_with_drug_class$genome_drug.resistant_phenotype %in% c("Resistant", "Susceptible"),
  ][["genome_drug.resistant_phenotype"]] |>
    table() |>
    as.data.frame()
  n_genomes_resistant <- total_strains[total_strains$Var1 == "Resistant", "Freq"]
  n_genomes_susceptible <- total_strains[total_strains$Var1 == "Susceptible", "Freq"]

  top_n_drugs <- data_with_drug_class %>%
    dplyr::group_by(genome_drug.antibiotic) %>%
    dplyr::summarize(n = n()) %>%
    dplyr::arrange(desc(n)) %>%
    dplyr::slice_head(n = 5) %>%
    dplyr::pull(genome_drug.antibiotic)

  top_n_drugs_class <- data_with_drug_class %>%
    dplyr::filter(!is.na(drug_class)) %>%
    dplyr::group_by(drug_class) %>%
    dplyr::summarize(n = n()) %>%
    dplyr::arrange(desc(n)) %>%
    dplyr::slice_head(n = 5) %>%
    dplyr::pull(drug_class)

  top_n_countries <- data_with_drug_class %>%
    dplyr::filter(!grepl(pattern = "NA", x = genome.isolation_country)) %>%
    dplyr::group_by(genome.isolation_country) %>%
    dplyr::summarize(n = n()) %>%
    dplyr::arrange(desc(n)) %>%
    dplyr::slice_head(n = 5) %>%
    dplyr::pull(genome.isolation_country)
  spp_name <- unique(data_with_drug_class$species)
  summary_paragraph <- tags$p(
    style = "text-align: center; font-weight: bold; padding: 12px 0; font-family: sans-serif; font-size: 20px;",
    stringr::str_glue("Data summary for {stringr::str_to_sentence(spp_name)}")
  )
  tagList(
    summary_paragraph,
    fluidRow(
      column(2, quickStatBox("# isolate-drug-AMR_phenotype combo", total_genomes,
        bg_color = "#000000", text_color = "white"
      )),
      column(2, quickStatBox("# unique genomes", total_uniques_genomes,
        bg_color = "#56B4E9", text_color = "black"
      )),
      column(2, quickStatBox("# resistant isolates", n_genomes_resistant,
        bg_color = "#009E73", text_color = "black"
      )),
      column(2, quickStatBox("# susceptible isolates", n_genomes_susceptible,
        bg_color = "#F0E442", text_color = "black"
      )),
      column(2, quickStatBox("# drugs", paste(n_amr_drugs, collapse = ", "),
        bg_color = "#0072B2"
      ), text_color = "white"),
      column(2, quickStatBox("# drug classes", paste(n_amr_drug_class, collapse = ", "),
        bg_color = "#D55E00", text_color = "white"
      )),
    ),
    fluidRow(
      column(4, quickStatBox("Top 5 drugs", paste(top_n_drugs, collapse = "\n"),
        bg_color = "#CC79A7", , text_color = "black"
      )),
      column(4, quickStatBox("Top 5 drug classes", paste(top_n_drugs_class, collapse = "\n"),
        bg_color = "#999999", text_color = "black"
      )),
      column(4, quickStatBox("Top countries", paste(top_n_countries, collapse = "\n"),
        bg_color = "#E69F00", text_color = "white"
      ))
    )
  )
}

makeDatAvailabilityPlot <- function(data) {
  data <- data %>%
    dplyr::group_by(genome_drug.antibiotic, genome_drug.resistant_phenotype) %>%
    count() %>%
    ungroup()
  ggplot(
    data,
    aes(
      x = genome_drug.antibiotic,
      y = n,
      color = genome_drug.resistant_phenotype,
      fill = genome_drug.resistant_phenotype
    )
  ) +
    geom_col() +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 10),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 10)
    ) +
    labs(
      x = "Drugs",
      y = "No. of isolates",
      color = "AMR phenotype",
      fill = "AMR phenotype"
    )
}

makeGeoChloroPlot <- function(data) {
  data$iso3 <- countrycode::countrycode(data$genome.isolation_country, origin = "country.name", destination = "iso3c")
  plot_ly(
    data = data,
    type = "choropleth",
    locations = ~iso3,
    z = ~count,
    text = ~ paste(iso3, "<br>Genomes:", count),
    hoverinfo = "text",
    colorscale = "Viridis",
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

  ggplot(
    data,
    aes(
      x = genome.collection_year,
      y = n,
      colour = genome_drug.resistant_phenotype
    )
  ) +
    geom_line() +
    geom_point() +
    labs(
      title = title_amr_drug,
      x = "Year",
      y = "No. of isolates",
      colour = "AMR phenotype"
    ) +
    theme_minimal() +
    theme(
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12)
    )
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

  ggplot(
    host_df,
    aes(
      x = genome_drug.antibiotic,
      y = n,
      fill = genome.host_common_name
    )
  ) +
    geom_col(position = "stack") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      legend.text = element_text(size = 10),
      axis.title = element_text(size = 10),
    ) +
    labs(
      x = "Drugs",
      y = "No. of isolates",
      fill = "Host"
    )
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

  isolation_source_plot <- ggplot(
    isolation_source,
    aes(
      x = genome_drug.antibiotic,
      y = n,
      fill = genome.isolation_source_
    )
  ) +
    geom_col(position = "stack") +
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
  return(isolation_source_plot)
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
  if (is.null(data) || !is.data.frame(data) || !nrow(data)) {
    return(ggplot() +
      labs(title = "No data available") +
      theme_bw())
  }

  # Filter to baseline models (strat_label is NA = no country/year stratification)
  df <- data %>%
    dplyr::filter(normalize_species(.data$species) %in% normalize_species(bug)) %>%
    dplyr::filter(.data$feature_type %in% model_scale) %>%
    dplyr::filter(.data$feature_subtype %in% data_type) %>%
    dplyr::filter(is.na(.data$strat_label) | !nzchar(.data$strat_label))

  # Filter by drug class or drug if not "all"
  if (!is.null(amr_drug_class) && length(amr_drug_class) > 0 &&
    !identical(amr_drug_class, "all")) {
    df <- df %>%
      dplyr::filter(
        (.data$drug_label == "drug_class" & .data$drug_or_class %in% amr_drug_class) |
          (.data$drug_label == "drug" & .data$drug_or_class %in% amr_drug)
      )
  }

  if (!nrow(df)) {
    return(ggplot() +
      labs(title = "No data for current selection") +
      theme_bw())
  }

  # Normalize species and set ESKAPE factor order
  eskape_order <- c("Efa", "Sau", "Kpn", "Aba", "Pae", "Esp")
  df <- df %>%
    dplyr::mutate(species = normalize_species(.data$species))
  present <- unique(df$species)
  lvls <- intersect(eskape_order, present)
  if (length(lvls) > 0) {
    df <- df %>% dplyr::mutate(species = factor(.data$species, levels = lvls))
  }

  metric_sym <- rlang::sym(metrics)
  g <- ggplot(
    df,
    aes(
      x = .data$species,
      y = !!metric_sym,
      fill = .data$feature_type,
      color = .data$feature_type
    )
  ) +
    geom_boxplot(alpha = 0.4) +
    labs(
      title = "Performance metrics",
      x = "Species",
      y = metrics,
      fill = "Scale",
      color = "Scale"
    ) +
    scale_fill_brewer(palette = "Set2") +
    scale_color_brewer(palette = "Dark2") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    coord_cartesian(ylim = c(0, 1))

  if (!is.null(amr_drug) && length(amr_drug) > 0 && nzchar(amr_drug[1])) {
    g <- g +
      geom_point(
        data = df %>% dplyr::filter(.data$drug_or_class %in% amr_drug),
        position = position_jitterdodge(jitter.width = 0.1),
        size = 4
      )
  }
  g
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
  annotated_dir = NULL,
  amrdata_root = NULL,
  results_root = NULL
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
    system.file("extdata", "Annotated", package = "amRviz")
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

  # Replace opaque feature IDs (e.g. "group_6367") with human-readable names
  # from {scale}_names.parquet when available.
  name_map <- load_feature_name_map(
    species_code = bug_norm[1],
    model_scale = model_scale,
    amrdata_root = amrdata_root,
    results_root = results_root
  )
  if (!is.null(name_map) && nrow(name_map)) {
    # Domain variables in top features look like "PF21279_IPR...": split
    # on "_" to extract the join key.
    join_key <- if (scale == "domain") {
      stringr::str_split_i(top_features_df$COG_name, "_", 1)
    } else {
      top_features_df$COG_name
    }
    lookup <- stats::setNames(name_map$label, name_map$Variable)
    new_label <- lookup[join_key]
    # Keep original id only when no label is available or it's blank
    replace <- !is.na(new_label) & nzchar(new_label)
    top_features_df$COG_name[replace] <- paste0(
      top_features_df$COG_name[replace],
      " (", new_label[replace], ")"
    )
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

  max_val <- max(vi_mat, na.rm = TRUE)
  min_val <- min(vi_mat, na.rm = TRUE)

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

# makeMetadataSankey: multi-tier sankey of resistance flow:
# phenotype -> drug class -> antibiotic -> country -> host -> isolation source.
# Filters to a chosen drug class (or top N classes) to keep the diagram
# legible, since unfiltered metadata has too many flows.
# data: metadata tibble (one row per genome-drug record)
# drug_classes: vector of drug class names to keep (NULL = top 3 by count)
makeMetadataSankey <- function(data, drug_classes = NULL,
                               max_classes = 3) {
  if (!requireNamespace("networkD3", quietly = TRUE)) {
    return(NULL)
  }
  if (is.null(data) || !is.data.frame(data) || !nrow(data)) {
    return(NULL)
  }

  required_cols <- c(
    "genome_drug.resistant_phenotype", "genome_drug.antibiotic",
    "genome.isolation_country", "genome.host_common_name",
    "genome.isolation_source", "drug_class"
  )
  if (!all(required_cols %in% names(data))) {
    return(NULL)
  }

  df <- data %>%
    dplyr::filter(
      !is.na(.data$genome_drug.resistant_phenotype),
      !is.na(.data$genome_drug.antibiotic),
      !is.na(.data$drug_class),
      !is.na(.data$genome.isolation_country),
      !is.na(.data$genome.host_common_name),
      !is.na(.data$genome.isolation_source),
      nzchar(.data$genome.isolation_country),
      nzchar(.data$genome.host_common_name),
      nzchar(.data$genome.isolation_source),
      .data$genome_drug.resistant_phenotype %in%
        c("Resistant", "Susceptible")
    )
  if (!nrow(df)) {
    return(NULL)
  }

  # Pick top drug classes if not explicitly given
  if (is.null(drug_classes) || !length(drug_classes)) {
    drug_classes <- df %>%
      dplyr::count(.data$drug_class, name = "n") %>%
      dplyr::arrange(dplyr::desc(.data$n)) %>%
      dplyr::slice_head(n = max_classes) %>%
      dplyr::pull(.data$drug_class)
  }
  df <- df %>% dplyr::filter(.data$drug_class %in% drug_classes)
  if (!nrow(df)) {
    return(NULL)
  }

  # Group rare isolation sources into "Other" to keep the diagram readable
  top_sources <- df %>%
    dplyr::count(.data$genome.isolation_source, name = "n") %>%
    dplyr::arrange(dplyr::desc(.data$n)) %>%
    dplyr::slice_head(n = 8) %>%
    dplyr::pull(.data$genome.isolation_source)
  df <- df %>%
    dplyr::mutate(
      genome.isolation_source = dplyr::if_else(
        .data$genome.isolation_source %in% top_sources,
        .data$genome.isolation_source,
        "Other"
      )
    )

  agg <- df %>%
    dplyr::count(
      .data$genome_drug.resistant_phenotype,
      .data$drug_class,
      .data$genome_drug.antibiotic,
      .data$genome.isolation_country,
      .data$genome.host_common_name,
      .data$genome.isolation_source,
      name = "n"
    )

  # Build node list (unique values across all tiers)
  nodes <- data.frame(
    name = unique(c(
      agg$genome_drug.resistant_phenotype,
      agg$drug_class,
      agg$genome_drug.antibiotic,
      agg$genome.isolation_country,
      agg$genome.host_common_name,
      agg$genome.isolation_source
    )),
    stringsAsFactors = FALSE
  )

  # Build links across each adjacent pair of tiers
  mk_links <- function(src_col, tgt_col, group_col) {
    sub <- agg
    sub$.src <- sub[[src_col]]
    sub$.tgt <- sub[[tgt_col]]
    sub$.grp <- sub[[group_col]]
    sub %>%
      dplyr::group_by(.data$.src, .data$.tgt, .data$.grp) %>%
      dplyr::summarise(value = sum(.data$n), .groups = "drop") %>%
      dplyr::transmute(
        source = match(.data$.src, nodes$name) - 1,
        target = match(.data$.tgt, nodes$name) - 1,
        value  = .data$value,
        group  = as.character(.data$.grp)
      )
  }

  links <- dplyr::bind_rows(
    mk_links("genome_drug.resistant_phenotype", "drug_class", "drug_class"),
    mk_links("drug_class", "genome_drug.antibiotic", "drug_class"),
    mk_links(
      "genome_drug.antibiotic", "genome.isolation_country",
      "drug_class"
    ),
    mk_links(
      "genome.isolation_country", "genome.host_common_name",
      "drug_class"
    ),
    mk_links(
      "genome.host_common_name", "genome.isolation_source",
      "drug_class"
    )
  )

  networkD3::sankeyNetwork(
    Links = as.data.frame(links),
    Nodes = nodes,
    Source = "source",
    Target = "target",
    Value = "value",
    NodeID = "name",
    LinkGroup = "group",
    fontSize = 12,
    nodeWidth = 25,
    sinksRight = FALSE
  )
}

# makeCogBarChart: horizontal bar chart of the most common COGs across the
# features in a top-features tibble (already enriched with annotations).
# top_n: number of COGs to display.
makeCogBarChart <- function(enriched_tbl, top_n = 15) {
  if (is.null(enriched_tbl) || !nrow(enriched_tbl) ||
    !"COG" %in% names(enriched_tbl)) {
    return(plotly::plot_ly() %>%
      plotly::layout(title = list(text = "No annotations available", x = 0)))
  }

  # Split the comma-separated COG cells and count occurrences per Variable.
  cog_df <- enriched_tbl %>%
    dplyr::filter(!is.na(.data$COG), nzchar(.data$COG)) %>%
    dplyr::select(.data$Variable, .data$COG, dplyr::any_of("COG_name")) %>%
    dplyr::distinct()

  if (!nrow(cog_df)) {
    return(plotly::plot_ly() %>%
      plotly::layout(title = list(text = "No COGs in selection", x = 0)))
  }

  rows <- do.call(rbind, lapply(seq_len(nrow(cog_df)), function(i) {
    cogs <- trimws(strsplit(cog_df$COG[i], ",", fixed = TRUE)[[1]])
    names <- if ("COG_name" %in% names(cog_df) &&
      !is.na(cog_df$COG_name[i])) {
      n <- trimws(strsplit(cog_df$COG_name[i], ";", fixed = TRUE)[[1]])
      rep_len(n, length(cogs))
    } else {
      rep(NA_character_, length(cogs))
    }
    data.frame(COG = cogs, COG_name = names, stringsAsFactors = FALSE)
  }))

  counts <- rows %>%
    dplyr::filter(nzchar(.data$COG)) %>%
    dplyr::count(.data$COG, .data$COG_name, name = "n") %>%
    dplyr::arrange(dplyr::desc(.data$n)) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::mutate(
      label = dplyr::if_else(
        is.na(.data$COG_name) | !nzchar(.data$COG_name),
        .data$COG,
        paste0(
          .data$COG, ": ",
          stringr::str_trunc(.data$COG_name, 40)
        )
      )
    )

  counts$label <- factor(counts$label, levels = rev(counts$label))

  plotly::plot_ly(
    data = counts,
    type = "bar",
    orientation = "h",
    x = ~n,
    y = ~label,
    marker = list(color = "#5b8db8"),
    hovertemplate = paste0(
      "<b>%{y}</b><br>Features: %{x}<extra></extra>"
    )
  ) %>%
    plotly::layout(
      title = list(
        text = "Top COGs among selected features",
        x = 0, font = list(size = 13)
      ),
      xaxis = list(title = "Features"),
      yaxis = list(title = ""),
      margin = list(l = 10, t = 50, r = 20, b = 40)
    )
}

# makeDrugFeatureNetwork: interactive force-directed graph linking drugs (or
# drug classes) to their top features (Variables). Optionally extends to
# cluster and COG tiers when an annotations parquet is available.
# top_data: pre-loaded top-features tibble from loadTopFeat().
# bug: 3-letter species code.
# top_n: number of top features per drug to include as edges.
# include_clusters / include_cogs: add annotation tiers when TRUE.
# results_root: path for annotation lookup (falls back to extdata).
makeDrugFeatureNetwork <- function(top_data, bug, top_n = 10,
                                   include_clusters = FALSE,
                                   include_cogs = FALSE,
                                   results_root = NULL) {
  if (!requireNamespace("networkD3", quietly = TRUE)) {
    stop(
      "Package 'networkD3' is required. ",
      "Install it with install.packages('networkD3')."
    )
  }
  if (is.null(top_data) || !is.data.frame(top_data) || !nrow(top_data)) {
    return(NULL)
  }

  df <- top_data %>%
    dplyr::mutate(species = normalize_species(.data$species)) %>%
    dplyr::filter(.data$species %in% normalize_species(bug)) %>%
    dplyr::filter(is.na(.data$strat_label) | !nzchar(.data$strat_label)) %>%
    dplyr::filter(!is.na(.data$Variable), !is.na(.data$drug_or_class))

  if (!nrow(df)) {
    return(NULL)
  }

  # Take top N features per drug/class by max importance
  edges_df <- df %>%
    dplyr::group_by(.data$drug_or_class, .data$Variable) %>%
    dplyr::summarise(
      Importance = max(.data$Importance, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::group_by(.data$drug_or_class) %>%
    dplyr::slice_max(order_by = .data$Importance, n = top_n) %>%
    dplyr::ungroup()

  if (!nrow(edges_df)) {
    return(NULL)
  }

  drugs <- unique(edges_df$drug_or_class)
  variables <- unique(edges_df$Variable)

  drug_var_edges <- edges_df %>%
    dplyr::transmute(
      source_name = .data$drug_or_class,
      target_name = .data$Variable,
      value = .data$Importance
    )

  extra_nodes_cluster <- character(0)
  extra_nodes_cog <- character(0)
  var_cluster_edges <- NULL
  cluster_cog_edges <- NULL

  if ((include_clusters || include_cogs)) {
    ann <- load_feature_annotations(bug, results_root)
    if (!is.null(ann) && nrow(ann)) {
      # Match Variable ("PF23840_IPR056912") against ann$feature ("PF23840")
      variables_key <- stringr::str_split_i(variables, "_", 1)
      var_key_map <- stats::setNames(variables, variables_key)
      ann_sub <- ann %>%
        dplyr::filter(.data$feature %in% variables_key)

      if (include_clusters && nrow(ann_sub)) {
        cl_edges <- ann_sub %>%
          dplyr::filter(!is.na(.data$cluster)) %>%
          dplyr::distinct(.data$feature, .data$cluster) %>%
          dplyr::transmute(
            source_name = var_key_map[.data$feature],
            target_name = .data$cluster,
            value = 0.5
          )
        if (nrow(cl_edges)) {
          var_cluster_edges <- cl_edges
          extra_nodes_cluster <- unique(cl_edges$target_name)
        }
      }

      if (include_cogs && nrow(ann_sub)) {
        if (include_clusters && !is.null(var_cluster_edges)) {
          # cluster -> COG edges
          cc_edges <- ann_sub %>%
            dplyr::filter(!is.na(.data$cluster), !is.na(.data$COG)) %>%
            dplyr::distinct(.data$cluster, .data$COG) %>%
            dplyr::transmute(
              source_name = .data$cluster,
              target_name = .data$COG,
              value = 0.3
            )
          if (nrow(cc_edges)) {
            cluster_cog_edges <- cc_edges
            extra_nodes_cog <- unique(cc_edges$target_name)
          }
        } else {
          # variable -> COG directly
          cc_edges <- ann_sub %>%
            dplyr::filter(!is.na(.data$COG)) %>%
            dplyr::distinct(.data$feature, .data$COG) %>%
            dplyr::transmute(
              source_name = var_key_map[.data$feature],
              target_name = .data$COG,
              value = 0.3
            )
          if (nrow(cc_edges)) {
            cluster_cog_edges <- cc_edges
            extra_nodes_cog <- unique(cc_edges$target_name)
          }
        }
      }
    }
  }

  nodes <- data.frame(
    name = c(drugs, variables, extra_nodes_cluster, extra_nodes_cog),
    node_type = c(
      rep("drug", length(drugs)),
      rep("variable", length(variables)),
      rep("cluster", length(extra_nodes_cluster)),
      rep("cog", length(extra_nodes_cog))
    ),
    stringsAsFactors = FALSE
  )

  all_edges <- dplyr::bind_rows(
    drug_var_edges, var_cluster_edges, cluster_cog_edges
  )
  edges <- all_edges %>%
    dplyr::transmute(
      source = match(.data$source_name, nodes$name) - 1,
      target = match(.data$target_name, nodes$name) - 1,
      value = .data$value
    )

  color_scale <- 'd3.scaleOrdinal()
    .domain(["drug","variable","cluster","cog"])
    .range(["#d4872a","#5b8db8","#66a61e","#9b7fba"])'

  networkD3::forceNetwork(
    Links = edges,
    Nodes = nodes,
    Source = "source",
    Target = "target",
    Value = "value",
    NodeID = "name",
    Group = "node_type",
    colourScale = networkD3::JS(color_scale),
    opacity = 0.9,
    zoom = TRUE,
    fontSize = 13,
    linkDistance = 100,
    charge = -80,
    legend = TRUE
  )
}

# makeCrossModelRidgePlot: balanced accuracy by drug class for holdout models,
# coloured by Same (self-eval) vs Different (cross-eval).
# cross_model: "country" or "time"
makeCrossModelRidgePlot <- function(perf_data, bug, cross_model) {
  if (is.null(perf_data) || !is.data.frame(perf_data) || !nrow(perf_data)) {
    return(plotly::plot_ly() %>%
      plotly::layout(title = list(text = "No data available", x = 0)))
  }

  strat <- if (cross_model == "country") "country" else "year"

  df <- perf_data %>%
    dplyr::filter(
      normalize_species(.data$species) %in% normalize_species(bug)
    ) %>%
    dplyr::filter(.data$strat_label == strat) %>%
    dplyr::filter(.data$drug_label == "drug_class")

  if (!nrow(df)) {
    return(plotly::plot_ly() %>%
      plotly::layout(title = list(text = "No data for selection", x = 0)))
  }

  # Label: Same = trained & tested on same stratum, Different = cross-tested
  df <- df %>%
    dplyr::mutate(
      test_type = dplyr::if_else(
        .data$cross_test, "Different", "Same"
      )
    )

  strat_label <- if (cross_model == "country") "Countries" else "Time periods"

  colors <- c(Same = "#d4872a", Different = "#5b8db8")

  p <- plotly::plot_ly()

  for (tt in c("Same", "Different")) {
    sub <- df[df$test_type == tt, ]
    if (!nrow(sub)) next
    col <- colors[[tt]]
    p <- p %>%
      plotly::add_trace(
        type = "box",
        x = sub$bal_acc,
        y = sub$drug_or_class,
        name = tt,
        orientation = "h",
        boxpoints = "all",
        jitter = 0.4,
        pointpos = 0,
        marker = list(
          color = col, opacity = 0.7, size = 6
        ),
        line = list(color = col),
        fillcolor = paste0(col, "44"),
        hovertemplate = paste0(
          "<b>Drug class:</b> %{y}<br>",
          "<b>Bal. Accuracy:</b> %{x:.3f}<br>",
          "<b>Test:</b> ", tt, "<extra></extra>"
        )
      )
  }

  p %>% plotly::layout(
    title = list(
      text = paste("Balanced accuracy by drug class -", strat_label),
      x = 0, font = list(size = 13)
    ),
    xaxis = list(title = "Balanced accuracy", range = c(0, 1.05)),
    yaxis = list(title = ""),
    legend = list(title = list(text = "Test")),
    boxmode = "group",
    margin = list(l = 100, t = 50)
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

  plotly::plot_ly(
    x = colnames(models_performance),
    y = rownames(models_performance),
    z = models_performance,
    type = "heatmap",
    colorscale = list(
      c(0, "#f2f0f7"), c(0.25, "#cbc9e2"),
      c(0.5, "#9e9ac8"), c(0.75, "#756bb1"),
      c(1, "#54278f")
    ),
    zmin = 0.5, zmax = 1.0,
    colorbar = list(title = "Balanced\naccuracy"),
    hovertemplate = paste0(
      "<b>Train data:</b> %{x}<br>",
      "<b>Test data:</b> %{y}<br>",
      "<b>Bal. Accuracy:</b> %{z:.3f}<extra></extra>"
    )
  ) %>%
    plotly::layout(
      title = list(
        text = if (cross_model == "country") {
          "Cross-country performance"
        } else {
          "Cross-time performance"
        },
        x = 0, font = list(size = 13)
      ),
      xaxis = list(title = "Train data", side = "bottom"),
      yaxis = list(title = "Test data", autorange = "reversed"),
      margin = list(l = 100, b = 60, t = 50)
    )
}

# Load a feature-id -> human-readable name mapping from the amRdata-style
# directory layout. Returns a tibble with columns `Variable` and `label`,
# or NULL if no matching parquet is found. Looks in:
#   {amrdata_root}/{species_dir}/{gene|domain|protein}_names.parquet
#   {results_root}/{species_dir}/{gene|domain|protein}_names.parquet
#   {extdata}/{species_dir}/{gene|domain|protein}_names.parquet
load_feature_name_map <- function(species_code, model_scale,
                                  amrdata_root = NULL,
                                  results_root = NULL) {
  scale <- dplyr::case_when(
    model_scale == "proteins" ~ "protein",
    model_scale == "domains" ~ "domain",
    model_scale == "genes" ~ "gene",
    TRUE ~ model_scale
  )
  fname <- paste0(scale, "_names.parquet")

  # Species-code and species-label maps use the same lookup paths as
  # cluster_feature_COG.parquet.
  roots <- c(
    amrdata_root, results_root,
    system.file("extdata", package = "amRviz")
  )
  roots <- roots[!is.null(roots) & nzchar(roots) & dir.exists(roots)]

  fp <- NULL
  for (r in roots) {
    for (d in list.dirs(r, full.names = TRUE, recursive = FALSE)) {
      cand <- file.path(d, fname)
      if (file.exists(cand)) {
        fp <- cand
        break
      }
    }
    if (!is.null(fp)) break
  }
  if (is.null(fp)) {
    return(NULL)
  }

  df <- tryCatch(arrow::read_parquet(fp), error = function(e) NULL)
  if (is.null(df) || !nrow(df)) {
    return(NULL)
  }

  # Normalise to {Variable, label}
  if (scale == "gene" && all(c("Gene", "Annotation") %in% names(df))) {
    return(tibble::tibble(
      Variable = df$Gene, label = df$Annotation
    ))
  }
  if (scale == "protein" && all(
    c("proteinID", "proteinName") %in% names(df)
  )) {
    return(tibble::tibble(
      Variable = df$proteinID, label = df$proteinName
    ))
  }
  if (scale == "domain" && all(c("DB.ID", "SignDesc") %in% names(df))) {
    # domain ids need deduping since one Pfam can occur many times
    agg <- df %>%
      dplyr::distinct(.data$DB.ID, .data$SignDesc) %>%
      dplyr::group_by(.data$DB.ID) %>%
      dplyr::summarise(
        label = dplyr::first(.data$SignDesc),
        .groups = "drop"
      )
    return(tibble::tibble(Variable = agg$DB.ID, label = agg$label))
  }
  NULL
}

# Load cluster/COG annotations for a species if the parquet exists.
# Searches in results_root/<species_dir>/cluster_feature_COG.parquet first,
# then falls back to extdata/<species_dir>/cluster_feature_COG.parquet.
load_feature_annotations <- function(species_code, results_root = NULL) {
  fname <- "cluster_feature_COG.parquet"
  if (!is.null(results_root) && nzchar(results_root)) {
    for (d in list.dirs(results_root, full.names = TRUE, recursive = FALSE)) {
      fp <- file.path(d, fname)
      if (file.exists(fp)) {
        return(arrow::read_parquet(fp))
      }
    }
  }
  extdata <- system.file("extdata", package = "amRviz")
  if (nzchar(extdata)) {
    for (d in list.dirs(extdata, full.names = TRUE, recursive = FALSE)) {
      fp <- file.path(d, fname)
      if (file.exists(fp)) {
        return(arrow::read_parquet(fp))
      }
    }
  }
  NULL
}

# Enrich a top-features tibble with cluster/COG annotations joined on
# Variable -> feature. Collapses multiple COGs per feature into one comma-
# separated cell. Returns the input unchanged if no annotations found.
enrich_with_annotations <- function(tbl, species_code, results_root = NULL) {
  if (is.null(tbl) || !nrow(tbl) || !"Variable" %in% names(tbl)) {
    return(tbl)
  }
  ann <- load_feature_annotations(species_code, results_root)
  if (is.null(ann) || !nrow(ann)) {
    return(tbl)
  }

  # Top features Variable can be "PF23840_IPR056912" (domains) or the raw
  # feature id (genes/proteins). Split on first "_" to extract the join key.
  tbl <- tbl %>%
    dplyr::mutate(
      .join_key = stringr::str_split_i(.data$Variable, "_", 1)
    )

  join_keys <- unique(tbl$.join_key)

  ann_collapsed <- ann %>%
    dplyr::filter(.data$feature %in% join_keys) %>%
    dplyr::group_by(.data$feature) %>%
    dplyr::summarise(
      cluster = dplyr::first(.data$cluster),
      cluster_name = dplyr::first(.data$cluster_name),
      COG = paste(unique(stats::na.omit(.data$COG)), collapse = ", "),
      COG_name = paste(unique(stats::na.omit(.data$COG_name)),
        collapse = "; "
      ),
      .groups = "drop"
    )

  tbl %>%
    dplyr::left_join(ann_collapsed, by = c(".join_key" = "feature")) %>%
    dplyr::select(-.data$.join_key)
}

makeFeatureImportTable <- function(feature_import_table) {
  # Early exit for empty/zero-column data
  if (is.null(feature_import_table) || ncol(feature_import_table) == 0) {
    return(DT::datatable(tibble::tibble(), options = list(dom = "t"), rownames = FALSE))
  }

  # Preferred display order (only those that exist will be used)
  cols_priority <- c(
    "species", "drug_or_class", "Variable",
    "cluster", "cluster_name", "COG", "COG_name",
    "COG_description", "ARG_name", "ARG_description",
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
  # Link cluster (fig IDs) to BVBRC, one <a> per unique id (comma-sep cells).
  if ("cluster" %in% names(tbl)) {
    link_fig <- function(ids) {
      if (is.na(ids) || !nzchar(ids)) {
        return(ids)
      }
      parts <- trimws(strsplit(ids, ",", fixed = TRUE)[[1]])
      linked <- vapply(parts, function(id) {
        url <- paste0(
          "https://www.bv-brc.org/view/Feature/",
          utils::URLencode(id, reserved = TRUE)
        )
        paste0(
          "<a href='", url, "' target='_blank' ",
          "style='color:#1a73e8; text-decoration: underline;'>",
          id, "</a>"
        )
      }, character(1))
      paste(linked, collapse = ", ")
    }
    tbl$cluster <- vapply(tbl$cluster, link_fig, character(1))
  }
  # Link COG ids (comma-separated) to NCBI COG page.
  if ("COG" %in% names(tbl)) {
    link_cog <- function(ids) {
      if (is.na(ids) || !nzchar(ids)) {
        return(ids)
      }
      parts <- trimws(strsplit(ids, ",", fixed = TRUE)[[1]])
      linked <- vapply(parts, function(id) {
        paste0(
          "<a href='https://www.ncbi.nlm.nih.gov/research/cog/cog/",
          id, "' target='_blank' ",
          "style='color:#1a73e8; text-decoration: underline;'>",
          id, "</a>"
        )
      }, character(1))
      paste(linked, collapse = ", ")
    }
    tbl$COG <- vapply(tbl$COG, link_cog, character(1))
  }

  DT::datatable(
    tbl,
    options = list(scrollX = TRUE, autoWidth = FALSE, orderClasses = TRUE),
    class = "display nowrap stripe",
    rownames = FALSE,
    width = "100%",
    escape = FALSE,
    selection = "single"
  )
}

# makeFeatureEgoNetwork: small force-directed graph for a single selected
# feature, showing feature -> cluster -> COG links.
makeFeatureEgoNetwork <- function(enriched_tbl, variable) {
  if (!requireNamespace("networkD3", quietly = TRUE)) {
    return(NULL)
  }
  if (is.null(enriched_tbl) || !nrow(enriched_tbl) ||
    is.null(variable) || !nzchar(variable)) {
    return(NULL)
  }

  row <- enriched_tbl %>%
    dplyr::filter(.data$Variable == variable) %>%
    dplyr::slice_head(n = 1)
  if (!nrow(row)) {
    return(NULL)
  }

  cogs <- character(0)
  if ("COG" %in% names(row) && !is.na(row$COG) && nzchar(row$COG)) {
    cogs <- trimws(strsplit(row$COG, ",", fixed = TRUE)[[1]])
  }
  cluster <- if ("cluster" %in% names(row) && !is.na(row$cluster) &&
    nzchar(row$cluster)) {
    row$cluster
  } else {
    NA_character_
  }

  if (is.na(cluster) && !length(cogs)) {
    return(NULL)
  }

  nodes_name <- c(variable)
  nodes_type <- c("variable")
  edges_src <- c()
  edges_tgt <- c()

  if (!is.na(cluster)) {
    nodes_name <- c(nodes_name, cluster)
    nodes_type <- c(nodes_type, "cluster")
    edges_src <- c(edges_src, variable)
    edges_tgt <- c(edges_tgt, cluster)
    for (cg in cogs) {
      nodes_name <- c(nodes_name, cg)
      nodes_type <- c(nodes_type, "cog")
      edges_src <- c(edges_src, cluster)
      edges_tgt <- c(edges_tgt, cg)
    }
  } else {
    for (cg in cogs) {
      nodes_name <- c(nodes_name, cg)
      nodes_type <- c(nodes_type, "cog")
      edges_src <- c(edges_src, variable)
      edges_tgt <- c(edges_tgt, cg)
    }
  }

  nodes <- data.frame(
    name = nodes_name, node_type = nodes_type,
    stringsAsFactors = FALSE
  )
  edges <- data.frame(
    source = match(edges_src, nodes$name) - 1,
    target = match(edges_tgt, nodes$name) - 1,
    value = 1
  )

  color_scale <- 'd3.scaleOrdinal()
    .domain(["variable","cluster","cog"])
    .range(["#5b8db8","#66a61e","#9b7fba"])'

  networkD3::forceNetwork(
    Links = edges,
    Nodes = nodes,
    Source = "source",
    Target = "target",
    Value = "value",
    NodeID = "name",
    Group = "node_type",
    colourScale = networkD3::JS(color_scale),
    opacity = 0.9,
    zoom = TRUE,
    fontSize = 12,
    linkDistance = 80,
    charge = -120,
    legend = TRUE
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
  extdata <- system.file("extdata", package = "amRviz")
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
  extdata <- system.file("extdata", package = "amRviz")
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
  extdata <- system.file("extdata", package = "amRviz")
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
