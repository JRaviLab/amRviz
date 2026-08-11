# Metadata tab visualisations.


#' Styled summary-statistic card
#'
#' @param title Caption shown beneath the value.
#' @param value Value (string or tag) shown prominently.
#' @param icon_name Font Awesome icon name.
#' @param bg_color Background color.
#' @param text_color Text color (default "white").
#' @return A shiny `div` styled as a stat card.
#' @keywords internal
#' @noRd
quickStatBox <- function(title, value, icon_name, bg_color,
                         text_color = "white") {
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


#' Summary statistic cards for the metadata tab
#'
#' Builds the data-summary header: totals (records, genomes, drugs, classes,
#' resistant/susceptible isolates) and the top-5 drugs, classes, and countries.
#'
#' @param data Metadata tibble (one row per genome-drug record).
#' @return A shiny `tagList` of stat cards.
#' @keywords internal
#' @noRd
makeQuickStats <- function(data) {
  data_with_drug_class <- data
  # Sample stat calculations
  total_genomes <- nrow(data_with_drug_class)
  total_uniques_genomes <- length(unique(data_with_drug_class$genome.genome_id))
  n_amr_drugs <- length(unique(data_with_drug_class$genome_drug.antibiotic))
  n_amr_drug_class <- length(unique(data_with_drug_class$drug_class[!is.na(data_with_drug_class$drug_class)]))
  total_strains <- data_with_drug_class[
    data_with_drug_class$genome_drug.resistant_phenotype %in% c("Resistant", "Susceptible"),
  ][["genome_drug.resistant_phenotype"]] |>
    table() |>
    as.data.frame()
  n_genomes_resistant <- total_strains[total_strains$Var1 == "Resistant", "Freq"]
  n_genomes_susceptible <- total_strains[total_strains$Var1 == "Susceptible", "Freq"]

  top_n_drugs <- data_with_drug_class |>
    dplyr::group_by(genome_drug.antibiotic) |>
    dplyr::summarize(n = n()) |>
    dplyr::arrange(desc(n)) |>
    dplyr::slice_head(n = 5) |>
    dplyr::pull(genome_drug.antibiotic)

  top_n_drugs_class <- data_with_drug_class |>
    dplyr::filter(!is.na(drug_class)) |>
    dplyr::group_by(drug_class) |>
    dplyr::summarize(n = n()) |>
    dplyr::arrange(desc(n)) |>
    dplyr::slice_head(n = 5) |>
    dplyr::pull(drug_class)

  top_n_countries <- data_with_drug_class |>
    dplyr::filter(!is.na(genome.isolation_country)) |>
    dplyr::group_by(genome.isolation_country) |>
    dplyr::summarize(n = n()) |>
    dplyr::arrange(desc(n)) |>
    dplyr::slice_head(n = 5) |>
    dplyr::pull(genome.isolation_country)
  spp_name <- unique(data_with_drug_class$species)
  summary_paragraph <- tags$p(
    style = "font-weight: bold; padding: 4px 4px 4px; font-family: sans-serif; font-size: 14px;",
    stringr::str_glue("Data summary for {stringr::str_to_sentence(spp_name[1])}")
  )
  tagList(
    summary_paragraph,
    fluidRow(
      column(4, quickStatBox("Isolate-drug records", total_genomes, "database", "#3c5a6f")),
      column(4, quickStatBox("Unique genomes", total_uniques_genomes, "dna", "#7aab6e")),
      column(4, quickStatBox("Drugs tested", n_amr_drugs, "pills", "#9b7fba"))
    ),
    fluidRow(
      column(4, quickStatBox("Resistant isolates", n_genomes_resistant, "virus", "#d4872a")),
      column(4, quickStatBox("Susceptible isolates", n_genomes_susceptible, "shield-halved", "#8a8a8a")),
      column(4, quickStatBox("Drug classes", n_amr_drug_class, "layer-group", "#8b6b7a"))
    ),
    fluidRow(
      column(4, quickStatBox(
        "Top 5 drugs",
        tags$span(
          style = "font-size:11px; line-height:1.5;",
          HTML(paste(top_n_drugs, collapse = "<br/>"))
        ),
        "star", "#4e9a9a"
      )),
      column(4, quickStatBox(
        "Top 5 drug classes",
        tags$span(
          style = "font-size:11px; line-height:1.5;",
          HTML(paste(top_n_drugs_class, collapse = "<br/>"))
        ),
        "list", "#6a6f4e"
      )),
      column(4, quickStatBox(
        "Top 5 countries",
        tags$span(
          style = "font-size:11px; line-height:1.5;",
          HTML(paste(top_n_countries, collapse = "<br/>"))
        ),
        "globe", "#d4735e"
      ))
    )
  )
}


#' Stacked bar of isolate counts by drug and phenotype
#'
#' @param data Metadata tibble.
#' @return A plotly stacked bar chart.
#' @keywords internal
#' @noRd
makeDataAvailabilityPlot <- function(data) {
  data <- data |>
    dplyr::distinct(genome.genome_id, genome_drug.antibiotic, genome_drug.resistant_phenotype) |>
    dplyr::group_by(genome_drug.antibiotic, genome_drug.resistant_phenotype) |>
    count() |>
    ungroup()
  g <- ggplot(
    data,
    aes(
      x = genome_drug.antibiotic,
      y = n,
      color = genome_drug.resistant_phenotype,
      fill = genome_drug.resistant_phenotype,
      text = paste0(
        "Drug: ", genome_drug.antibiotic,
        "<br>Phenotype: ", genome_drug.resistant_phenotype,
        "<br>Isolates: ", n
      )
    )
  ) +
    geom_col() +
    scale_fill_manual(values = PHENOTYPE_COLORS, na.value = "gray70") +
    scale_color_manual(values = PHENOTYPE_COLORS, na.value = "gray70") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 10)
    ) +
    labs(
      x = "Drug",
      y = "No. of isolates",
      color = "AMR phenotype",
      fill = "AMR phenotype"
    )
  plotly::ggplotly(g, tooltip = "text")
}


#' World choropleth of genome counts by country
#'
#' @param data Tibble with `genome.isolation_country` and `count` columns.
#' @return A plotly choropleth map.
#' @keywords internal
#' @noRd
makeGeoChloroPlot <- function(data) {
  data$iso3 <- countrycode::countrycode(data$genome.isolation_country, origin = "country.name", destination = "iso3c")
  plot_ly(
    data = data,
    type = "choropleth",
    locations = ~iso3,
    z = ~count,
    text = ~ paste(iso3, "<br>Genomes:", count),
    hoverinfo = "text",
    colorscale = list(
      c(0, "#fdf3e6"),
      c(0.25, "#f0d3a5"),
      c(0.5, "#e2ad6a"),
      c(0.75, "#d4872a"),
      c(1, "#9a5e1c")
    ),
    marker = list(line = list(width = 0.5, color = "white"))
  ) |>
    colorbar(title = list(text = "No. of isolates", font = list(size = 14))) |>
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


#' Resistance trend over collection year
#'
#' @param data Tibble with collection year, phenotype, and isolate counts.
#' @param amr_drug Drug to title the plot with, or "all".
#' @return A plotly line/point time series.
#' @keywords internal
#' @noRd
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

  g <- ggplot(
    data,
    aes(
      x = genome.collection_year,
      y = n,
      colour = genome_drug.resistant_phenotype,
      text = paste0(
        "Year: ", genome.collection_year,
        "<br>Phenotype: ", genome_drug.resistant_phenotype,
        "<br>Isolates: ", n
      )
    )
  ) +
    geom_line(aes(group = genome_drug.resistant_phenotype)) +
    geom_point() +
    scale_color_manual(values = PHENOTYPE_COLORS, na.value = "gray70") +
    labs(
      title = title_amr_drug,
      x = "Year",
      y = "No. of isolates",
      colour = "AMR phenotype"
    ) +
    theme_bw() +
    theme(
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      axis.text = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_text(size = 10)
    )
  plotly::ggplotly(g, tooltip = "text")
}


#' Stacked bar of isolate counts by drug and host
#'
#' @param data Metadata tibble.
#' @return A plotly stacked bar chart.
#' @keywords internal
#' @noRd
makeHostIsolatePlot <- function(data) {
  data <- data |>
    dplyr::mutate(
      genome.host_common_name = stringr::str_to_lower(genome.host_common_name),
      genome.isolation_source = stringr::str_to_lower(genome.isolation_source)
    )

  host_df <- data |>
    dplyr::group_by(genome_drug.antibiotic, genome.host_common_name) |>
    dplyr::summarize(n = n()) |>
    dplyr::ungroup()

  n_hosts <- length(unique(host_df$genome.host_common_name))

  g <- ggplot(
    host_df,
    aes(
      x = genome_drug.antibiotic,
      y = n,
      fill = genome.host_common_name,
      text = paste0(
        "Drug: ", genome_drug.antibiotic,
        "<br>Host: ", genome.host_common_name,
        "<br>Isolates: ", n
      )
    )
  ) +
    geom_col(position = "stack") +
    scale_fill_manual(values = meta_palette(n_hosts), na.value = "gray70") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      axis.title = element_text(size = 10)
    ) +
    labs(
      x = "Drug",
      y = "No. of isolates",
      fill = "Host"
    )
  plotly::ggplotly(g, tooltip = "text")
}


#' Stacked bar of isolate counts by drug and isolation source
#'
#' Keeps the top 10 isolation sources and groups the rest into "Other".
#'
#' @param data Metadata tibble.
#' @return A plotly stacked bar chart.
#' @keywords internal
#' @noRd
makeIsolationSourcesPlot <- function(data) {
  isolation_source <- data |>
    dplyr::mutate(genome.isolation_source = stringr::str_to_lower(genome.isolation_source)) |>
    dplyr::group_by(genome_drug.antibiotic, genome.isolation_source) |>
    dplyr::summarize(n = n()) |>
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
  isolation_source <- isolation_source |>
    mutate(genome.isolation_source_ = stringr::str_trunc(genome.isolation_source_, width = 20))

  n_sources <- length(unique(isolation_source$genome.isolation_source_))

  isolation_source_plot <- ggplot(
    isolation_source,
    aes(
      x = genome_drug.antibiotic,
      y = n,
      fill = genome.isolation_source_,
      text = paste0(
        "Drug: ", genome_drug.antibiotic,
        "<br>Source: ", genome.isolation_source_,
        "<br>Isolates: ", n
      )
    )
  ) +
    geom_col(position = "stack") +
    scale_fill_manual(values = meta_palette(n_sources), na.value = "gray70") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      axis.title = element_text(size = 10)
    ) +
    labs(
      x = "Drug",
      y = "No. of isolates",
      fill = "Isolation source"
    )
  plotly::ggplotly(isolation_source_plot, tooltip = "text")
}


#' Multi-tier resistance-flow Sankey
#'
#' Sankey of resistance flow:
#' phenotype -> drug class -> antibiotic -> country -> host -> isolation source.
#' Filters to the chosen drug classes (or the top `max_classes` by count) to
#' keep the diagram legible, since unfiltered metadata has too many flows.
#'
#' @param data Metadata tibble (one row per genome-drug record).
#' @param drug_classes Drug classes to keep (NULL = top `max_classes` by count).
#' @param max_classes Number of top drug classes to keep when none are given.
#' @return A `networkD3` sankeyNetwork widget, or NULL when prerequisites are
#'   missing.
#' @keywords internal
#' @noRd
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

  df <- data |>
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
    drug_classes <- df |>
      dplyr::count(.data$drug_class, name = "n") |>
      dplyr::arrange(dplyr::desc(.data$n)) |>
      dplyr::slice_head(n = max_classes) |>
      dplyr::pull(.data$drug_class)
  }
  df <- df |> dplyr::filter(.data$drug_class %in% drug_classes)
  if (!nrow(df)) {
    return(NULL)
  }

  # Group rare isolation sources into "Other" to keep the diagram readable
  top_sources <- df |>
    dplyr::count(.data$genome.isolation_source, name = "n") |>
    dplyr::arrange(dplyr::desc(.data$n)) |>
    dplyr::slice_head(n = 8) |>
    dplyr::pull(.data$genome.isolation_source)
  df <- df |>
    dplyr::mutate(
      genome.isolation_source = dplyr::if_else(
        .data$genome.isolation_source %in% top_sources,
        .data$genome.isolation_source,
        "Other"
      )
    )

  agg <- df |>
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
    sub |>
      dplyr::group_by(.data$.src, .data$.tgt, .data$.grp) |>
      dplyr::summarise(value = sum(.data$n), .groups = "drop") |>
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

  # Sankey-specific muted palette. Resistant -> amber and Susceptible -> grey
  # are explicitly mapped; remaining range colors are deliberately lower-
  # saturation than META_COLORS so the translucent link bands don't read as
  # bright pastels at the Sankey's default opacity.
  range_colors <- c(
    "#d4872a", # Resistant (amber)
    "#8a8a8a", # Susceptible (grey)
    "#7a8a9b", # muted slate blue
    "#c4a991", # dusty peach
    "#9ab18a", # sage green
    "#a99c8f", # warm taupe
    "#8b6b7a", # mauve
    "#4e9a9a", # teal
    "#c4a35a", # mustard
    "#9b7fba", # dusty purple
    "#d4735e", # terracotta
    "#6a6f4e" # olive
  )
  colour_scale <- networkD3::JS(paste0(
    'd3.scaleOrdinal().domain(["Resistant","Susceptible"]).range(["',
    paste(range_colors, collapse = '","'),
    '"])'
  ))

  networkD3::sankeyNetwork(
    Links = as.data.frame(links),
    Nodes = nodes,
    Source = "source",
    Target = "target",
    Value = "value",
    NodeID = "name",
    LinkGroup = "group",
    colourScale = colour_scale,
    fontSize = 12,
    nodeWidth = 25,
    sinksRight = FALSE
  )
}
