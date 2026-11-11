# Force-directed network visualizations.
#
# The drug-feature network draws one table, built by networkTables(): its
# nodes and links. The ego network is drawn on the Features tab.


#' Nodes and links for the drug-feature network
#'
#' Links each drug or drug class to its `top_n` most important features in
#' the baseline models (not stratified or cross-tested), and optionally each
#' feature to its protein cluster when cluster annotations are available.
#'
#' @param top_data Top-features tibble from loadTopFeat().
#' @param bug Species code.
#' @param top_n Number of top features per drug or drug class.
#' @param include_clusters Add cluster nodes from the annotations when TRUE.
#' @param results_root Path for the annotation lookup.
#' @return A list of `nodes` (`id`, unique across node types, `label` and
#'   `group`: "drug", a feature type, or "cluster") and `links` (`from` and
#'   `to` node ids, and `importance`), or NULL when there is nothing to draw.
#' @keywords internal
#' @noRd
networkTables <- function(top_data, bug, top_n = 10,
                          include_clusters = FALSE, results_root = NULL) {
  if (is.null(top_data) || !is.data.frame(top_data) || !nrow(top_data)) {
    return(NULL)
  }

  df <- top_data |>
    dplyr::mutate(species = normalize_species(.data$species)) |>
    dplyr::filter(.data$species %in% normalize_species(bug)) |>
    dplyr::filter(is.na(.data$strat_label) | !nzchar(.data$strat_label)) |>
    dplyr::filter(!.data$cross_test %in% TRUE) |>
    dplyr::filter(!is.na(.data$Variable), !is.na(.data$drug_or_class)) |>
    dplyr::mutate(
      Variable = dplyr::case_when(
        # Older amRml named domains "PF23840_IPR056912"; keep the Pfam ID.
        # Current Pfam names ("DDE_Tnp_Tn3") are kept whole.
        feature_type == "domains" & grepl("^PF[0-9]+_", Variable) ~
          sub("_.+$", "", Variable),
        feature_type == "proteins" ~ sub("fig.", "fig|", Variable, fixed = TRUE),
        feature_type == "args" ~ sub(
          "^X", "",
          gsub("\\.NCBIFAM", "", Variable)
        ),
        TRUE ~ Variable
      )
    )
  if (!nrow(df)) {
    return(NULL)
  }

  # Top N features per drug or class, by their highest importance across
  # seeds and encodings.
  top <- df |>
    dplyr::group_by(.data$drug_or_class, .data$Variable, .data$feature_type) |>
    dplyr::summarise(
      importance = max(.data$Importance, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::group_by(.data$drug_or_class) |>
    dplyr::slice_max(order_by = .data$importance, n = top_n) |>
    dplyr::ungroup()
  if (!nrow(top)) {
    return(NULL)
  }

  # A node's id includes its type, so a feature named like a drug, or the
  # same name in two feature types, stays a separate node.
  drug_id <- paste0("drug:", top$drug_or_class)
  feature_id <- paste0(top$feature_type, ":", top$Variable)
  nodes <- dplyr::bind_rows(
    tibble::tibble(id = drug_id, label = top$drug_or_class, group = "drug"),
    tibble::tibble(id = feature_id, label = top$Variable, group = top$feature_type)
  ) |>
    dplyr::distinct(.data$id, .keep_all = TRUE)
  links <- tibble::tibble(from = drug_id, to = feature_id, importance = top$importance)

  if (include_clusters) {
    ann <- load_feature_annotations(bug, results_root)
    if (!is.null(ann) && nrow(ann)) {
      clusters <- ann |>
        dplyr::filter(.data$feature %in% top$Variable, !is.na(.data$cluster)) |>
        dplyr::distinct(.data$feature, .data$cluster)
      if (nrow(clusters)) {
        features <- nodes[nodes$group != "drug", ]
        from <- features$id[match(clusters$feature, features$label)]
        nodes <- dplyr::bind_rows(nodes, tibble::tibble(
          id = paste0("cluster:", unique(clusters$cluster)),
          label = unique(clusters$cluster), group = "cluster"
        ))
        links <- dplyr::bind_rows(links, tibble::tibble(
          from = from, to = paste0("cluster:", clusters$cluster), importance = 0.5
        ))
      }
    }
  }

  list(nodes = nodes, links = links)
}


#' Interactive force-directed drug -> feature network
#'
#' Draws networkTables(): drugs (or drug classes) linked to their top
#' features, optionally extended to protein clusters.
#'
#' @inheritParams networkTables
#' @return A `networkD3` forceNetwork widget, or NULL when there is nothing to
#'   plot.
#' @keywords internal
#' @noRd
makeDrugFeatureNetwork <- function(top_data, bug, top_n = 10,
                                   include_clusters = FALSE,
                                   results_root = NULL) {
  if (!requireNamespace("networkD3", quietly = TRUE)) {
    stop(
      "Package 'networkD3' is required. ",
      "Install it with install.packages('networkD3')."
    )
  }
  net <- networkTables(top_data, bug, top_n, include_clusters, results_root)
  if (is.null(net)) {
    return(NULL)
  }
  nodes <- as.data.frame(net$nodes)
  edges <- data.frame(
    source = match(net$links$from, nodes$id) - 1,
    target = match(net$links$to, nodes$id) - 1,
    value = net$links$importance
  )

  # Colour domain: drug + every molecular scale + cluster, drawing scale colours
  # from the shared SCALE_COLORS palette so the network matches the other tabs.
  js_array <- function(x) paste0('["', paste(x, collapse = '","'), '"]')
  group_domain <- c("drug", names(SCALE_COLORS), "cluster")
  group_range <- c("#d4872a", unname(SCALE_COLORS), "#66a61e")
  color_scale <- sprintf(
    "d3.scaleOrdinal().domain(%s).range(%s)",
    js_array(group_domain), js_array(group_range)
  )

  .styleNetwork(networkD3::forceNetwork(
    Links = edges,
    Nodes = nodes,
    Source = "source",
    Target = "target",
    Value = "value",
    NodeID = "label",
    Group = "group",
    colourScale = networkD3::JS(color_scale),
    opacity = 0.9,
    zoom = TRUE,
    fontSize = 13,
    fontFamily = .NETWORK_FONT,
    linkDistance = 100,
    charge = -80,
    legend = TRUE
  ))
}


#' Ego network for a single feature
#'
#' Small force-directed graph for one selected feature, showing its
#' feature -> cluster -> COG links.
#'
#' @param enriched_tbl Annotation-enriched top-features tibble (see
#'   enrich_with_annotations()).
#' @param variable The Variable (feature id) to centre the graph on.
#' @return A `networkD3` forceNetwork widget, or NULL when the feature has no
#'   cluster/COG links to show.
#' @keywords internal
#' @noRd
makeFeatureEgoNetwork <- function(enriched_tbl, variable) {
  if (!requireNamespace("networkD3", quietly = TRUE)) {
    return(NULL)
  }
  if (is.null(enriched_tbl) || !nrow(enriched_tbl) ||
    is.null(variable) || !nzchar(variable)) {
    return(NULL)
  }

  row <- enriched_tbl |>
    dplyr::filter(.data$Variable == variable) |>
    dplyr::slice_head(n = 1)
  if (!nrow(row)) {
    return(NULL)
  }

  # cogs <- character(0)
  # if ("COG" %in% names(row) && !is.na(row$COG) && nzchar(row$COG)) {
  #   cogs <- trimws(strsplit(row$COG, ",", fixed = TRUE)[[1]])
  # }
  cluster <- if ("cluster" %in% names(row) && !is.na(row$cluster) &&
    nzchar(row$cluster)) {
    row$cluster
  } else {
    NA_character_
  }

  if (is.na(cluster)) {
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
    # for (cg in cogs) {
    #   nodes_name <- c(nodes_name, cg)
    #   nodes_type <- c(nodes_type, "cog")
    #   edges_src <- c(edges_src, cluster)
    #   edges_tgt <- c(edges_tgt, cg)
    # }
  } else {
    # for (cg in cogs) {
    #   nodes_name <- c(nodes_name, cg)
    #   nodes_type <- c(nodes_type, "cog")
    #   edges_src <- c(edges_src, variable)
    #   edges_tgt <- c(edges_tgt, cg)
    # }
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
    .domain(["variable","cluster"])
    .range(["#5b8db8","#66a61e"])'

  .styleNetwork(networkD3::forceNetwork(
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
    fontFamily = .NETWORK_FONT,
    linkDistance = 80,
    charge = -120,
    legend = TRUE
  ))
}


# Font family for the force-directed networks, matching the plotly/ggplot
# panels (which use "Arial, sans-serif").
.NETWORK_FONT <- "Arial, sans-serif"

# Post-render hook for forceNetwork: replace each node's default circle with a
# shape keyed to its node_type group so element types are distinguishable
# beyond colour (drug = diamond, variable = circle, cluster = square,
# cog = triangle). Wrapped in try/catch so the network falls back to circles
# if the D3 internals change.
.NETWORK_SHAPE_JS <- "
function(el, x) {
  try {
    var d3 = window.d3;
    if (!d3 || !d3.symbol) { return; }
    var shapeMap = {
      'drug': d3.symbolDiamond,
      'variable': d3.symbolCircle,
      'cluster': d3.symbolSquare
    };
    var sym = d3.symbol();
    d3.select(el).selectAll('.node').each(function(d) {
      var g = d3.select(this);
      var circle = g.select('circle');
      if (circle.empty()) { return; }
      var fill = circle.style('fill');
      var r = parseFloat(circle.attr('r')) || 8;
      var type = shapeMap[d.group] || d3.symbolCircle;
      circle.style('display', 'none');
      g.insert('path', ':first-child')
        .attr('class', 'node-shape')
        .attr('d', sym.type(type).size(Math.PI * r * r * 1.6)())
        .style('fill', fill)
        .style('opacity', 0.9)
        .style('stroke', '#ffffff')
        .style('stroke-width', 1.2);
    });
    // Mirror the shapes in the legend: swap each colour rect for its glyph.
    var legendSym = d3.symbol().size(150);
    d3.select(el).selectAll('.legend').each(function(d) {
      var g = d3.select(this);
      var rect = g.select('rect');
      if (rect.empty()) { return; }
      var fill = rect.style('fill');
      var type = shapeMap[d] || d3.symbolCircle;
      rect.style('display', 'none');
      g.insert('path', ':first-child')
        .attr('class', 'legend-shape')
        .attr('transform', 'translate(9,9)')
        .attr('d', legendSym.type(type)())
        .style('fill', fill)
        .style('stroke', fill);
    });
  } catch (e) { /* keep default circles on error */ }
}
"

# Apply the shared font + per-group node shapes to a forceNetwork widget.
.styleNetwork <- function(net) {
  if (is.null(net)) {
    return(NULL)
  }
  htmlwidgets::onRender(net, .NETWORK_SHAPE_JS)
}
