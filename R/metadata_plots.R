# Metadata tab data prep and visualizations, shared by the dashboard and
# exportAMRVisualizations().
#
# Each plot draws one table, built by a function in the "Tables" section
# below (metadataSummary(), metadataAvailability(), metadataByCountry(),
# metadataByYear(), metadataByHost(), metadataBySource(), metadataSankey()).
# One metadata row is one genome tested against one drug, so a table counts
# either genomes or genome-drug records, and its columns say which.


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
  shiny::div(
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
    shiny::div(
      shiny::div(
        style = "font-size: 18px; font-weight: bold; line-height: 1.2;",
        value
      ),
      shiny::div(
        style = "font-size: 11px; margin-top: 3px; opacity: 0.9;",
        title
      )
    ),
    shiny::div(
      style = "font-size: 24px; opacity: 0.4;",
      shiny::icon(icon_name)
    )
  )
}


#' Summary statistic cards for the metadata tab
#'
#' Draws metadataSummary(): totals (records, genomes, drugs, classes,
#' resistant and susceptible records) and the top-5 drugs, classes, and
#' countries.
#'
#' @param data Metadata tibble (one row per genome-drug record).
#' @return A shiny `tagList` of stat cards.
#' @keywords internal
#' @noRd
makeQuickStats <- function(data) {
  s <- metadataSummary(data)
  top_list <- function(x) {
    shiny::tags$span(
      style = "font-size:11px; line-height:1.5;",
      shiny::HTML(paste(x, collapse = "<br/>"))
    )
  }
  summary_paragraph <- shiny::tags$p(
    style = "font-weight: bold; padding: 4px 4px 4px; font-family: sans-serif; font-size: 14px;",
    stringr::str_glue("Data summary for {stringr::str_to_sentence(s$species)}")
  )
  shiny::tagList(
    summary_paragraph,
    shiny::fluidRow(
      shiny::column(4, quickStatBox("Genome-drug records", s$records, "database", "#3c5a6f")),
      shiny::column(4, quickStatBox("Unique genomes", s$genomes, "dna", "#7aab6e")),
      shiny::column(4, quickStatBox("Drugs tested", s$drugs, "pills", "#9b7fba"))
    ),
    shiny::fluidRow(
      shiny::column(4, quickStatBox("Resistant records", s$resistant_records, "virus", "#d4872a")),
      shiny::column(4, quickStatBox("Susceptible records", s$susceptible_records, "shield-halved", "#8a8a8a")),
      shiny::column(4, quickStatBox("Drug classes", s$drug_classes, "layer-group", "#8b6b7a"))
    ),
    shiny::fluidRow(
      shiny::column(4, quickStatBox("Top 5 drugs", top_list(s$top_drugs), "star", "#4e9a9a")),
      shiny::column(4, quickStatBox("Top 5 drug classes", top_list(s$top_drug_classes), "list", "#6a6f4e")),
      shiny::column(4, quickStatBox("Top 5 countries", top_list(s$top_countries), "globe", "#d4735e"))
    )
  )
}


#' Shared look for the metadata bar and line charts
#' @keywords internal
#' @noRd
.metadata_theme <- function(vertical_x = TRUE) {
  ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = if (vertical_x) {
        ggplot2::element_text(angle = 90, hjust = 1, size = 10)
      } else {
        ggplot2::element_text(angle = 45, hjust = 1, size = 10)
      },
      axis.text.y = ggplot2::element_text(size = 10),
      legend.title = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 10),
      axis.title = ggplot2::element_text(size = 10)
    )
}


#' Stacked bar of genomes tested per drug, by phenotype
#'
#' @param data Output of metadataAvailability().
#' @return A plotly stacked bar chart.
#' @keywords internal
#' @noRd
makeDataAvailabilityPlot <- function(data) {
  if (is.null(data) || !nrow(data)) {
    return(plotly_placeholder("No phenotype data"))
  }
  g <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data$drug,
      y = .data$genomes,
      color = .data$phenotype,
      fill = .data$phenotype,
      text = paste0(
        "Drug: ", .data$drug,
        "<br>Phenotype: ", .data$phenotype,
        "<br>Genomes: ", .data$genomes
      )
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = PHENOTYPE_COLORS, na.value = "gray70") +
    ggplot2::scale_color_manual(values = PHENOTYPE_COLORS, na.value = "gray70") +
    .metadata_theme() +
    ggplot2::labs(
      x = "Drug", y = "No. of genomes",
      color = "AMR phenotype", fill = "AMR phenotype"
    )
  plotly::ggplotly(g, tooltip = "text")
}


#' World choropleth of genomes by country
#'
#' @param data Output of metadataByCountry().
#' @return A plotly choropleth map.
#' @keywords internal
#' @noRd
makeGeoChloroPlot <- function(data) {
  data <- data[!is.na(data$iso3), , drop = FALSE]
  if (!nrow(data)) {
    return(plotly_placeholder("No country data"))
  }
  plotly::plot_ly(
    data = data,
    type = "choropleth",
    locations = ~iso3,
    z = ~genomes,
    text = ~ paste0(
      country, "<br>Genomes: ", genomes, "<br>Genome-drug records: ", records
    ),
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
    plotly::colorbar(title = list(text = "No. of genomes", font = list(size = 14))) |>
    plotly::layout(
      geo = list(
        projection = list(type = "natural earth"),
        showcoastlines = TRUE,
        showland = TRUE,
        landcolor = "lightgrey"
      ),
      font = list(size = 14)
    )
}


#' Resistant and susceptible records per collection year
#'
#' @param data Output of metadataByYear().
#' @param amr_drug Drug to title the plot with, or "all".
#' @return A plotly line/point time series.
#' @keywords internal
#' @noRd
makeTimeSeriesAMRPlot <- function(data, amr_drug) {
  if (is.null(data) || !nrow(data)) {
    return(plotly_placeholder("No collection years recorded"))
  }
  title <- if (identical(amr_drug, "all")) {
    "AMR phenotypes by collection year"
  } else {
    stringr::str_glue("AMR phenotypes by collection year for {amr_drug}")
  }
  g <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data$year,
      y = .data$records,
      colour = .data$phenotype,
      text = paste0(
        "Year: ", .data$year,
        "<br>Phenotype: ", .data$phenotype,
        "<br>Genome-drug records: ", .data$records
      )
    )
  ) +
    ggplot2::geom_line(ggplot2::aes(group = .data$phenotype)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_manual(values = PHENOTYPE_COLORS, na.value = "gray70") +
    ggplot2::labs(
      title = title, x = "Collection year", y = "Genome-drug records",
      colour = "AMR phenotype"
    ) +
    .metadata_theme(vertical_x = FALSE)
  plotly::ggplotly(g, tooltip = "text")
}


#' Stacked bar of records per drug, by host
#'
#' @param data Output of metadataByHost().
#' @return A plotly stacked bar chart.
#' @keywords internal
#' @noRd
makeHostIsolatePlot <- function(data) {
  if (is.null(data) || !nrow(data)) {
    return(plotly_placeholder("No hosts recorded"))
  }
  g <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data$drug,
      y = .data$records,
      fill = .data$host,
      text = paste0(
        "Drug: ", .data$drug,
        "<br>Host: ", .data$host,
        "<br>Genome-drug records: ", .data$records
      )
    )
  ) +
    ggplot2::geom_col(position = "stack") +
    ggplot2::scale_fill_manual(
      values = meta_palette(length(unique(data$host))), na.value = "gray70"
    ) +
    .metadata_theme() +
    ggplot2::labs(x = "Drug", y = "Genome-drug records", fill = "Host")
  plotly::ggplotly(g, tooltip = "text")
}


#' Stacked bar of records per drug, by isolation source
#'
#' @param data Output of metadataBySource().
#' @return A plotly stacked bar chart.
#' @keywords internal
#' @noRd
makeIsolationSourcesPlot <- function(data) {
  if (is.null(data) || !nrow(data)) {
    return(plotly_placeholder("No isolation sources recorded"))
  }
  sources <- levels(data$source)
  named <- setdiff(sources, c("other", "not recorded"))
  colors <- c(
    stats::setNames(meta_palette(length(named)), named),
    other = "#b5b5b5", `not recorded` = "#e0e0e0"
  )
  g <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data$drug,
      y = .data$records,
      fill = .data$source,
      text = paste0(
        "Drug: ", .data$drug,
        "<br>Source: ", .data$source,
        "<br>Genome-drug records: ", .data$records
      )
    )
  ) +
    ggplot2::geom_col(position = "stack") +
    ggplot2::scale_fill_manual(
      values = colors[sources],
      labels = stringr::str_trunc(sources, width = 24)
    ) +
    .metadata_theme() +
    ggplot2::labs(x = "Drug", y = "Genome-drug records", fill = "Isolation source")
  plotly::ggplotly(g, tooltip = "text")
}


#' Multi-tier resistance-flow Sankey
#'
#' Draws metadataSankey(): phenotype -> drug class -> antibiotic -> country ->
#' host -> isolation source.
#'
#' @param data Metadata tibble (one row per genome-drug record).
#' @param drug_classes Drug classes to keep (NULL = top `max_classes` by count).
#' @param max_classes Number of top drug classes to keep when none are given.
#' @return A `networkD3` sankeyNetwork widget, or NULL when there is nothing to
#'   draw.
#' @keywords internal
#' @noRd
makeMetadataSankey <- function(data, drug_classes = NULL,
                               max_classes = 3) {
  flow <- metadataSankey(data, drug_classes, max_classes)
  if (is.null(flow)) {
    return(NULL)
  }
  nodes <- as.data.frame(flow$nodes)
  links <- data.frame(
    source = match(flow$links$from, nodes$id) - 1,
    target = match(flow$links$to, nodes$id) - 1,
    value = flow$links$records,
    group = flow$links$drug_class
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
    Links = links,
    Nodes = nodes,
    Source = "source",
    Target = "target",
    Value = "value",
    NodeID = "label",
    LinkGroup = "group",
    colourScale = colour_scale,
    fontSize = 12,
    nodeWidth = 25,
    sinksRight = FALSE
  )
}


# ---- Reading metadata ----

#' List the species that have metadata
#'
#' The ORBs whose manifests are found under `results_root` (species and
#' metadata file from each manifest). Without `results_root`, the packaged demo
#' ORB.
#'
#' @param results_root User results root, or NULL for the packaged demo data.
#' @return A sorted named character vector (names = display labels, values =
#'   metadata parquet paths), or empty when none are found.
#' @keywords internal
#' @noRd
listMetadataSpecies <- function(results_root = NULL) {
  .orb_choices(.metadata_orbs(results_root))
}


#' The ORBs that listMetadataSpecies() lists
#'
#' @param results_root User results root, or NULL for the packaged demo data.
#' @return A `.list_available_orbs()` tibble, with each ORB's species code
#'   (`dataset_id`), label and metadata path.
#' @keywords internal
#' @noRd
.metadata_orbs <- function(results_root = NULL) {
  root <- if (!is.null(results_root) && nzchar(results_root)) {
    results_root
  } else {
    system.file("extdata", package = "amRviz")
  }
  .list_available_orbs(.find_orb_manifests(root))
}


#' Read and row-bind metadata parquets for one or more species
#'
#' @param sources Named character vector of metadata parquet paths, as
#'   listMetadataSpecies() returns; names are the species labels.
#' @return A tibble with a `species` column set to each source's label, or an
#'   empty tibble when no metadata parquet is found.
#' @keywords internal
#' @noRd
readSpeciesMetadata <- function(sources) {
  labels <- names(sources) %||% rep(NA_character_, length(sources))
  purrr::map_dfr(seq_along(sources), function(i) {
    fp <- sources[[i]]
    if (!is.na(fp) && file.exists(fp)) {
      arrow::read_parquet(fp) |> dplyr::mutate(species = labels[[i]])
    } else {
      tibble::tibble()
    }
  })
}



# ---- Tables ----

#' Tidy a free-text category: trimmed, lower case, blanks as NA
#'
#' Hosts and isolation sources are typed by hand in BV-BRC, so "Feces",
#' "feces " and "feces" are the same category.
#'
#' @param x Character vector.
#' @return `x` lower-cased with runs of spaces squished, and NA where blank.
#' @keywords internal
#' @noRd
.clean_category <- function(x) {
  x <- stringr::str_squish(stringr::str_to_lower(as.character(x)))
  x[!is.na(x) & !nzchar(x)] <- NA_character_
  x
}


#' Numbers behind the summary cards
#'
#' @param meta Metadata tibble from readSpeciesMetadata().
#' @return A list: `species`, `records` (genome-drug records), `genomes`,
#'   `drugs`, `drug_classes`, `resistant_records`, `susceptible_records`, and
#'   `top_drugs`, `top_drug_classes`, `top_countries` (up to 5 each, by
#'   records).
#' @keywords internal
#' @noRd
metadataSummary <- function(meta) {
  top <- function(x, n = 5) {
    x <- x[!is.na(x) & nzchar(trimws(x))]
    names(utils::head(sort(table(x), decreasing = TRUE), n))
  }
  phenotype <- meta$genome_drug.resistant_phenotype
  list(
    species = unique(meta$species)[1] %||% NA_character_,
    records = nrow(meta),
    genomes = length(unique(meta$genome.genome_id)),
    drugs = length(unique(meta$genome_drug.antibiotic)),
    drug_classes = length(unique(stats::na.omit(meta$drug_class))),
    resistant_records = sum(phenotype %in% "Resistant"),
    susceptible_records = sum(phenotype %in% "Susceptible"),
    top_drugs = top(meta$genome_drug.antibiotic),
    top_drug_classes = top(meta$drug_class),
    top_countries = top(filter_valid_country(meta)$genome.isolation_country)
  )
}


#' Genomes tested per drug and phenotype
#'
#' @param meta Metadata tibble from readSpeciesMetadata().
#' @return A tibble with `drug`, `phenotype` and `genomes`, the input to
#'   makeDataAvailabilityPlot().
#' @keywords internal
#' @noRd
metadataAvailability <- function(meta) {
  meta |>
    dplyr::distinct(
      .data$genome.genome_id, .data$genome_drug.antibiotic,
      .data$genome_drug.resistant_phenotype
    ) |>
    dplyr::count(
      drug = .data$genome_drug.antibiotic,
      phenotype = .data$genome_drug.resistant_phenotype,
      name = "genomes"
    )
}


#' Genomes and records per isolation country
#'
#' Countries are matched to ISO3 codes for the map; regions that aren't
#' countries (for example "Southeast Asia") get NA and aren't drawn.
#'
#' @param meta Metadata tibble from readSpeciesMetadata().
#' @return A tibble with `country`, `iso3`, `genomes` and `records`, the input
#'   to makeGeoChloroPlot().
#' @keywords internal
#' @noRd
metadataByCountry <- function(meta) {
  meta |>
    filter_valid_country() |>
    dplyr::group_by(country = .data$genome.isolation_country) |>
    dplyr::summarise(
      genomes = dplyr::n_distinct(.data$genome.genome_id),
      records = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      iso3 = suppressWarnings(countrycode::countrycode(
        .data$country,
        origin = "country.name", destination = "iso3c"
      )),
      .after = "country"
    )
}


#' Records per collection year and phenotype
#'
#' @param meta Metadata tibble from readSpeciesMetadata().
#' @return A tibble with `year`, `phenotype` and `records`, the input to
#'   makeTimeSeriesAMRPlot(). Rows without a collection year are left out.
#' @keywords internal
#' @noRd
metadataByYear <- function(meta) {
  meta |>
    dplyr::filter(!is.na(.data$genome.collection_year)) |>
    dplyr::count(
      year = .data$genome.collection_year,
      phenotype = .data$genome_drug.resistant_phenotype,
      name = "records"
    )
}


#' Records per drug and host
#'
#' @param meta Metadata tibble from readSpeciesMetadata().
#' @return A tibble with `drug`, `host` (lower case) and `records`, the input
#'   to makeHostIsolatePlot(). Rows without a host are left out.
#' @keywords internal
#' @noRd
metadataByHost <- function(meta) {
  meta |>
    dplyr::mutate(host = .clean_category(.data$genome.host_common_name)) |>
    dplyr::filter(!is.na(.data$host)) |>
    dplyr::count(drug = .data$genome_drug.antibiotic, .data$host, name = "records")
}


#' Records per drug and isolation source
#'
#' Keeps the `top_n` most recorded sources and groups the rest as "other";
#' rows without a source are "not recorded", so they aren't mistaken for
#' other sources.
#'
#' @param meta Metadata tibble from readSpeciesMetadata().
#' @param top_n Number of sources shown by name.
#' @return A tibble with `drug`, `source` (a factor: named sources by records,
#'   then "other", then "not recorded") and `records`, the input to
#'   makeIsolationSourcesPlot().
#' @keywords internal
#' @noRd
metadataBySource <- function(meta, top_n = 10) {
  source <- .clean_category(meta$genome.isolation_source)
  ranked <- names(sort(table(source), decreasing = TRUE))
  named <- utils::head(ranked, top_n)
  source <- dplyr::case_when(
    is.na(source) ~ "not recorded",
    source %in% named ~ source,
    TRUE ~ "other"
  )
  tibble::tibble(drug = meta$genome_drug.antibiotic, source = source) |>
    dplyr::count(.data$drug, .data$source, name = "records") |>
    dplyr::mutate(source = factor(
      .data$source,
      levels = intersect(c(named, "other", "not recorded"), .data$source)
    ))
}


#' Nodes and links for the resistance-flow Sankey
#'
#' Flows from phenotype to drug class, drug, country, host and isolation
#' source, for the chosen drug classes (or the top `max_classes` by records).
#' Records missing any tier are left out. Isolation sources beyond the 8 most
#' recorded are grouped as "other".
#'
#' @param meta Metadata tibble (one row per genome-drug record).
#' @param drug_classes Drug classes to keep (NULL = top `max_classes`).
#' @param max_classes Number of top drug classes to keep when none are given.
#' @return A list of `nodes` (`id`, unique across tiers, `label` and `tier`)
#'   and `links` (`from` and `to` node ids, `records`, and the `drug_class`
#'   that colors them), or NULL when there is nothing to draw.
#' @keywords internal
#' @noRd
metadataSankey <- function(meta, drug_classes = NULL, max_classes = 3) {
  required <- c(
    "genome_drug.resistant_phenotype", "genome_drug.antibiotic",
    "genome.isolation_country", "genome.host_common_name",
    "genome.isolation_source", "drug_class"
  )
  if (is.null(meta) || !is.data.frame(meta) || !nrow(meta) ||
    !all(required %in% names(meta))) {
    return(NULL)
  }
  df <- meta |>
    filter_valid_country() |>
    dplyr::transmute(
      phenotype = .data$genome_drug.resistant_phenotype,
      drug_class = .data$drug_class,
      drug = .data$genome_drug.antibiotic,
      country = .data$genome.isolation_country,
      host = .clean_category(.data$genome.host_common_name),
      source = .clean_category(.data$genome.isolation_source)
    ) |>
    dplyr::filter(
      .data$phenotype %in% c("Resistant", "Susceptible"),
      !is.na(.data$drug_class), !is.na(.data$drug),
      !is.na(.data$host), !is.na(.data$source)
    )
  if (!nrow(df)) {
    return(NULL)
  }

  if (is.null(drug_classes) || !length(drug_classes)) {
    drug_classes <- names(utils::head(
      sort(table(df$drug_class), decreasing = TRUE), max_classes
    ))
  }
  df <- df[df$drug_class %in% drug_classes, , drop = FALSE]
  if (!nrow(df)) {
    return(NULL)
  }
  top_sources <- names(utils::head(sort(table(df$source), decreasing = TRUE), 8))
  df$source[!df$source %in% top_sources] <- "other"

  # A node's id includes its tier, so a value appearing in two tiers (a host
  # and an isolation source both called "human") stays two nodes.
  tiers <- c("phenotype", "drug_class", "drug", "country", "host", "source")
  id <- function(tier, value) paste0(tier, ":", value)
  links <- dplyr::bind_rows(lapply(seq_len(length(tiers) - 1), function(i) {
    df |>
      dplyr::count(
        from = id(tiers[i], .data[[tiers[i]]]),
        to = id(tiers[i + 1], .data[[tiers[i + 1]]]),
        .data$drug_class,
        name = "records"
      )
  }))
  nodes <- dplyr::bind_rows(lapply(tiers, function(tier) {
    tibble::tibble(label = unique(df[[tier]]), tier = tier)
  })) |>
    dplyr::mutate(id = id(.data$tier, .data$label), .before = 1)

  list(nodes = nodes, links = links[c("from", "to", "records", "drug_class")])
}


#' Most-represented drug classes, the default sankey selection
#'
#' @param meta Metadata tibble from readSpeciesMetadata().
#' @param n Number of classes to return.
#' @return A character vector of up to `n` drug classes, or NULL when `meta`
#'   has no `drug_class` column or no rows.
#' @keywords internal
#' @noRd
topDrugClasses <- function(meta, n = 3) {
  if (!nrow(meta) || !"drug_class" %in% names(meta)) {
    return(NULL)
  }
  meta |>
    dplyr::filter(!is.na(.data$drug_class)) |>
    dplyr::count(.data$drug_class, name = "records") |>
    dplyr::arrange(dplyr::desc(.data$records)) |>
    dplyr::slice_head(n = n) |>
    dplyr::pull(.data$drug_class)
}
