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


#' Share of resistant, intermediate and susceptible results per group
#'
#' Horizontal bars filled to 100%, most resistant on top. The resistant part
#' starts at 0%, so its length reads straight off the axis.
#'
#' @param data Output of metadataResistance().
#' @param unit What `total` counts, for the hover text and axis
#'   ("genome-drug records" or "genomes").
#' @param top_n Number of groups shown, by `total`.
#' @return A plotly horizontal bar chart.
#' @keywords internal
#' @noRd
makeResistanceSharePlot <- function(data, unit = "genome-drug records",
                                    top_n = 30) {
  if (is.null(data) || !nrow(data)) {
    return(plotly_placeholder("No phenotype data"))
  }
  shown <- utils::head(data[order(-data$total), , drop = FALSE], top_n)
  shown$label <- gsub("_", " ", shown$label, fixed = TRUE)
  order <- shown$label[order(shown$pct_resistant, shown$total)]
  long <- tidyr::pivot_longer(
    shown, c("resistant", "intermediate", "susceptible"),
    names_to = "phenotype", values_to = "n"
  ) |>
    dplyr::filter(.data$n > 0) |>
    dplyr::mutate(
      phenotype = factor(
        stringr::str_to_sentence(.data$phenotype),
        levels = c("Susceptible", "Intermediate", "Resistant")
      ),
      share = .data$n / .data$total,
      label = factor(.data$label, levels = order)
    )
  g <- ggplot2::ggplot(
    long,
    ggplot2::aes(
      x = .data$share,
      y = .data$label,
      fill = .data$phenotype,
      text = paste0(
        .data$label,
        "<br>", .data$phenotype, ": ", .data$n, " of ", .data$total,
        " ", unit, " (", round(100 * .data$share), "%)"
      )
    )
  ) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::scale_x_continuous(
      labels = function(x) paste0(round(100 * x), "%"),
      breaks = seq(0, 1, 0.25), expand = c(0, 0)
    ) +
    ggplot2::scale_fill_manual(values = PHENOTYPE_COLORS, drop = TRUE) +
    .metadata_theme(vertical_x = FALSE) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5)) +
    ggplot2::labs(x = paste("Share of", unit), y = NULL, fill = "AMR phenotype")
  plotly::ggplotly(g, tooltip = "text") |>
    plotly::layout(legend = list(traceorder = "reversed"))
}


#' Rows x collection years heatmap, shared by the two metadata heatmaps
#'
#' Years sit at evenly spaced positions (labeled by year), so a gap in the
#' years doesn't stretch the plot, and plotly can't read them as numbers.
#' Cells with no data are blank. With more than `max_labeled_years` years the
#' numbers in the cells are left to the hover text, and about 20 years are
#' labeled on the axis.
#'
#' @param rows Row of each cell (a factor: its levels give the order, bottom
#'   to top).
#' @param years Year of each cell.
#' @param value Color value of each cell.
#' @param label Number written in each cell.
#' @param hover Hover text of each cell.
#' @param colorbar,colorscale,zmin,zmax Passed to plotly.
#' @param max_labeled_years Most years whose cells still show their number.
#' @return A plotly heatmap.
#' @keywords internal
#' @noRd
.year_heatmap <- function(rows, years, value, label, hover, colorbar,
                          colorscale, zmin = NULL, zmax = NULL,
                          max_labeled_years = 25) {
  row_levels <- levels(droplevels(rows))
  year_levels <- sort(unique(years))
  shape <- c(length(row_levels), length(year_levels))
  idx <- cbind(match(as.character(rows), row_levels), match(years, year_levels))
  z <- matrix(NA_real_, shape[1], shape[2])
  text <- matrix("", shape[1], shape[2])
  tip <- matrix("", shape[1], shape[2])
  z[idx] <- value
  text[idx] <- as.character(label)
  tip[idx] <- hover
  if (length(year_levels) > max_labeled_years) text[] <- ""
  ticks <- seq(1, length(year_levels), by = ceiling(length(year_levels) / 20))
  plotly::plot_ly(
    x = seq_along(year_levels), y = seq_along(row_levels), z = z,
    type = "heatmap", zmin = zmin, zmax = zmax, xgap = 1, ygap = 1,
    colorscale = colorscale, colorbar = colorbar,
    text = text, texttemplate = "%{text}", textfont = list(size = 10),
    hovertext = tip, hoverinfo = "text"
  ) |>
    plotly::layout(
      xaxis = list(
        title = "Collection year", tickmode = "array",
        tickvals = ticks, ticktext = year_levels[ticks],
        tickangle = -45, showgrid = FALSE, zeroline = FALSE
      ),
      yaxis = list(
        title = "", tickmode = "array",
        tickvals = seq_along(row_levels), ticktext = row_levels,
        showgrid = FALSE, zeroline = FALSE
      ),
      plot_bgcolor = "#ffffff"
    )
}


#' Heatmap of the share of resistant records per drug and year
#'
#' @param data Output of metadataDrugYear().
#' @return A plotly heatmap: color is the share resistant, from susceptible
#'   blue to resistant amber, and each cell shows its number of records.
#'   Years with no records for a drug are blank.
#' @keywords internal
#' @noRd
makeDrugYearHeatmap <- function(data) {
  if (is.null(data) || !nrow(data) || dplyr::n_distinct(data$year) < 2) {
    return(plotly_placeholder("Too few collection years recorded"))
  }
  drug <- factor(
    gsub("_", " ", data$drug, fixed = TRUE),
    levels = gsub("_", " ", levels(data$drug), fixed = TRUE)
  )
  .year_heatmap(
    rows = drug, years = data$year, value = data$pct_resistant,
    label = data$total,
    hover = paste0(
      drug, ", ", data$year,
      "<br>Resistant: ", data$resistant, " of ", data$total,
      " genome-drug records (", round(100 * data$pct_resistant), "%)"
    ),
    colorbar = list(title = "% resistant", tickformat = ".0%"),
    colorscale = list(
      c(0, PHENOTYPE_COLORS[["Susceptible"]]),
      c(0.5, "#f4f1ea"),
      c(1, PHENOTYPE_COLORS[["Resistant"]])
    ),
    zmin = 0, zmax = 1
  )
}


#' Heatmap of genomes per country and collection year
#'
#' @param data Output of metadataCountryYear().
#' @return A plotly heatmap of genome counts; each cell shows its count.
#' @keywords internal
#' @noRd
makeCountryYearHeatmap <- function(data) {
  if (is.null(data) || !nrow(data) ||
    dplyr::n_distinct(data$year) < 2 || dplyr::n_distinct(data$country) < 2) {
    return(plotly_placeholder("Too few countries or years recorded"))
  }
  .year_heatmap(
    rows = data$country, years = data$year, value = data$genomes,
    label = data$genomes,
    hover = paste0(data$country, ", ", data$year, "<br>Genomes: ", data$genomes),
    colorbar = list(title = "Genomes"),
    colorscale = list(c(0, "#fdf3e6"), c(0.5, "#e2ad6a"), c(1, "#9a5e1c"))
  )
}


#' Bar chart of genomes by number of resistant drug classes
#'
#' @param data Output of metadataMdr().
#' @return A plotly bar chart, with multidrug-resistant genomes in amber and
#'   the share of them in the title.
#' @keywords internal
#' @noRd
makeMdrPlot <- function(data) {
  if (is.null(data) || !nrow(data) || !sum(data$genomes)) {
    return(plotly_placeholder("No resistant-class counts recorded"))
  }
  mdr_label <- paste0("Multidrug resistant (", .MDR_MIN_CLASSES, " or more classes)")
  data$group <- factor(
    ifelse(data$mdr, mdr_label, "Not multidrug resistant"),
    levels = c("Not multidrug resistant", mdr_label)
  )
  n_mdr <- sum(data$genomes[data$mdr])
  total <- sum(data$genomes)
  g <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = factor(.data$classes),
      y = .data$genomes,
      fill = .data$group,
      text = paste0(
        "Resistant to ", .data$classes, " drug class",
        ifelse(.data$classes == 1, "", "es"),
        "<br>Genomes: ", .data$genomes
      )
    )
  ) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::scale_fill_manual(
      values = stats::setNames(
        c(PHENOTYPE_COLORS[["Susceptible"]], PHENOTYPE_COLORS[["Resistant"]]),
        c("Not multidrug resistant", mdr_label)
      ),
      drop = FALSE
    ) +
    .metadata_theme(vertical_x = FALSE) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5)) +
    ggplot2::labs(
      title = sprintf(
        "%s of %s genomes (%d%%) are multidrug resistant",
        n_mdr, total, round(100 * n_mdr / total)
      ),
      x = "Number of drug classes resistant to", y = "Genomes", fill = NULL
    )
  plotly::ggplotly(g, tooltip = "text") |>
    plotly::layout(legend = list(orientation = "h", y = -0.25))
}


#' Bar chart of the most common resistance profiles
#'
#' @param data Output of metadataProfiles().
#' @return A plotly horizontal bar chart, most common profile on top.
#' @keywords internal
#' @noRd
makeProfilesPlot <- function(data) {
  if (is.null(data) || !nrow(data)) {
    return(plotly_placeholder("No resistance profiles recorded"))
  }
  data$profile <- factor(data$profile, levels = rev(data$profile))
  g <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data$genomes,
      y = .data$profile,
      text = paste0(
        "Resistant to: ", .data$profile,
        "<br>Drug classes: ", .data$classes,
        "<br>Genomes: ", .data$genomes
      )
    )
  ) +
    ggplot2::geom_col(fill = PHENOTYPE_COLORS[["Resistant"]], width = 0.75) +
    .metadata_theme(vertical_x = FALSE) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5)) +
    ggplot2::labs(x = "Genomes", y = "Resistant drug classes")
  plotly::ggplotly(g, tooltip = "text")
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


# ---- Resistance tables (shares, heatmaps, multidrug resistance) ----

# Fewest records a drug needs to get its own row in the drug x year heatmap.
.METADATA_MIN_RECORDS <- 10L

# A genome is multidrug resistant (MDR) when it is resistant to at least this
# many drug classes (Magiorakos et al. 2012).
.MDR_MIN_CLASSES <- 3L


#' Resistant, intermediate and susceptible counts per group
#'
#' @param df Tibble with a `phenotype` column and the grouping columns.
#' @param ... Grouping columns.
#' @return One row per group, with `resistant`, `intermediate`,
#'   `susceptible`, `total` and `pct_resistant` (0 to 1).
#' @keywords internal
#' @noRd
.phenotype_counts <- function(df, ...) {
  df |>
    dplyr::group_by(...) |>
    dplyr::summarise(
      resistant = sum(.data$phenotype == "Resistant"),
      intermediate = sum(.data$phenotype == "Intermediate"),
      susceptible = sum(.data$phenotype == "Susceptible"),
      total = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::mutate(pct_resistant = .data$resistant / .data$total)
}


#' Records with a resistant, intermediate or susceptible phenotype
#' @keywords internal
#' @noRd
.phenotype_records <- function(meta) {
  meta |>
    dplyr::mutate(phenotype = .data$genome_drug.resistant_phenotype) |>
    dplyr::filter(.data$phenotype %in% c("Resistant", "Intermediate", "Susceptible"))
}


#' One row per genome: its year, country and resistant drug classes
#'
#' Genome-level fields are the same on each of a genome's records, so the
#' first recorded value is kept.
#'
#' @param meta Metadata tibble from readSpeciesMetadata().
#' @return A tibble with `genome`, `year`, `country`, `classes` (number of
#'   resistant drug classes) and `profile` (the resistant classes, joined by
#'   "_", as amRdata writes them).
#' @keywords internal
#' @noRd
.metadata_genomes <- function(meta) {
  first <- function(x) x[!is.na(x)][1]
  col <- function(name, empty) if (name %in% names(meta)) meta[[name]] else empty
  country <- as.character(col("genome.isolation_country", NA_character_))
  country[!is.na(country) & !nzchar(trimws(country))] <- NA_character_
  profile <- as.character(col("resistant_classes", NA_character_))
  profile[!is.na(profile) & !nzchar(trimws(profile))] <- NA_character_
  tibble::tibble(
    genome = meta$genome.genome_id,
    year = suppressWarnings(as.numeric(col("genome.collection_year", NA_real_))),
    country = country,
    classes = suppressWarnings(as.numeric(col("num_resistant_classes", NA_real_))),
    profile = profile
  ) |>
    dplyr::group_by(.data$genome) |>
    dplyr::summarise(dplyr::across(dplyr::everything(), first), .groups = "drop")
}


#' Share of resistant results per drug or per drug class
#'
#' Per drug, each genome-drug record counts once. Per drug class, a genome
#' counts once per class: resistant if it is resistant to any drug in the
#' class, otherwise intermediate if any is intermediate, otherwise
#' susceptible.
#'
#' @param meta Metadata tibble from readSpeciesMetadata().
#' @param by "drug" or "class".
#' @return A tibble with `label` (drug or drug class), `resistant`,
#'   `intermediate`, `susceptible`, `total` (records per drug, genomes per
#'   class) and `pct_resistant`, most resistant first. The input to
#'   makeResistanceSharePlot().
#' @keywords internal
#' @noRd
metadataResistance <- function(meta, by = c("drug", "class")) {
  by <- match.arg(by)
  df <- .phenotype_records(meta)
  if (by == "class") {
    if (!"drug_class" %in% names(df)) {
      return(.phenotype_counts(df[0, ], label = character(0)))
    }
    df <- df |>
      dplyr::filter(!is.na(.data$drug_class), nzchar(.data$drug_class)) |>
      dplyr::group_by(.data$genome.genome_id, label = .data$drug_class) |>
      dplyr::summarise(
        phenotype = dplyr::case_when(
          any(.data$phenotype == "Resistant") ~ "Resistant",
          any(.data$phenotype == "Intermediate") ~ "Intermediate",
          TRUE ~ "Susceptible"
        ),
        .groups = "drop"
      )
  } else {
    df <- dplyr::mutate(df, label = .data$genome_drug.antibiotic)
  }
  df |>
    dplyr::filter(!is.na(.data$label)) |>
    .phenotype_counts(.data$label) |>
    dplyr::arrange(dplyr::desc(.data$pct_resistant), dplyr::desc(.data$total))
}


#' Share of resistant records per drug and collection year
#'
#' @param meta Metadata tibble from readSpeciesMetadata().
#' @param min_records Drugs with fewer dated records than this are left out.
#' @return A tibble with `drug`, `year`, `resistant`, `total` (records) and
#'   `pct_resistant`, the input to makeDrugYearHeatmap(). `drug` is a factor,
#'   least resistant overall first. Rows without a year are left out.
#' @keywords internal
#' @noRd
metadataDrugYear <- function(meta, min_records = .METADATA_MIN_RECORDS) {
  dated <- .phenotype_records(meta) |>
    dplyr::mutate(
      drug = .data$genome_drug.antibiotic,
      year = suppressWarnings(as.numeric(.data$genome.collection_year))
    ) |>
    dplyr::filter(!is.na(.data$year), !is.na(.data$drug))
  overall <- .phenotype_counts(dated, .data$drug) |>
    dplyr::filter(.data$total >= min_records) |>
    dplyr::arrange(.data$pct_resistant, .data$total)
  dated |>
    dplyr::filter(.data$drug %in% overall$drug) |>
    .phenotype_counts(.data$drug, .data$year) |>
    dplyr::select("drug", "year", "resistant", "total", "pct_resistant") |>
    dplyr::mutate(drug = factor(.data$drug, levels = overall$drug))
}


#' Genomes per isolation country and collection year
#'
#' @param meta Metadata tibble from readSpeciesMetadata().
#' @param top_n Number of countries shown by name; the rest are "other".
#' @return A tibble with `country` (a factor: fewest genomes first, "other"
#'   before all of them), `year` and `genomes`, the input to
#'   makeCountryYearHeatmap(). Genomes without a country or year are left out.
#' @keywords internal
#' @noRd
metadataCountryYear <- function(meta, top_n = 15) {
  g <- .metadata_genomes(meta) |>
    dplyr::filter(!is.na(.data$country), !is.na(.data$year))
  ranked <- names(sort(table(g$country), decreasing = TRUE))
  named <- utils::head(ranked, top_n)
  g$country <- ifelse(g$country %in% named, g$country, "other")
  levels <- rev(intersect(c(named, "other"), g$country))
  g |>
    dplyr::count(.data$country, .data$year, name = "genomes") |>
    dplyr::mutate(country = factor(.data$country, levels = levels))
}


#' Genomes by number of resistant drug classes
#'
#' @param meta Metadata tibble from readSpeciesMetadata().
#' @return A tibble with `classes` (0 to the most seen), `genomes` and `mdr`
#'   (TRUE at `.MDR_MIN_CLASSES` or more), the input to makeMdrPlot(). Empty
#'   when the metadata has no resistant-class counts.
#' @keywords internal
#' @noRd
metadataMdr <- function(meta) {
  g <- .metadata_genomes(meta) |>
    dplyr::filter(!is.na(.data$classes))
  if (!nrow(g)) {
    return(tibble::tibble(classes = numeric(0), genomes = integer(0), mdr = logical(0)))
  }
  tibble::tibble(classes = seq(0, max(g$classes))) |>
    dplyr::left_join(dplyr::count(g, .data$classes, name = "genomes"), by = "classes") |>
    dplyr::mutate(
      genomes = dplyr::coalesce(.data$genomes, 0L),
      mdr = .data$classes >= .MDR_MIN_CLASSES
    )
}


#' Most common combinations of resistant drug classes
#'
#' @param meta Metadata tibble from readSpeciesMetadata().
#' @param top_n Number of combinations kept.
#' @return A tibble with `profile` (class abbreviations joined by " + "),
#'   `classes` and `genomes`, most genomes first, the input to
#'   makeProfilesPlot(). Genomes resistant to no class are left out.
#' @keywords internal
#' @noRd
metadataProfiles <- function(meta, top_n = 15) {
  .metadata_genomes(meta) |>
    dplyr::filter(!is.na(.data$profile), !is.na(.data$classes), .data$classes > 0) |>
    dplyr::count(.data$profile, .data$classes, name = "genomes", sort = TRUE) |>
    utils::head(top_n) |>
    dplyr::mutate(profile = gsub("_", " + ", .data$profile, fixed = TRUE))
}
