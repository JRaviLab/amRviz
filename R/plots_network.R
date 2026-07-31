# Force-directed network visualisations.


#' Interactive force-directed drug -> feature network
#'
#' Links drugs (or drug classes) to their top features (Variables), optionally
#' extending to cluster and COG tiers when an annotations parquet is available.
#'
#' @param top_data Pre-loaded top-features tibble from loadTopFeat().
#' @param bug 3-letter species code.
#' @param top_n Number of top features per drug to include as edges.
#' @param include_clusters,include_cogs Add annotation tiers when TRUE.
#' @param results_root Path for annotation lookup (falls back to extdata).
#' @return A `networkD3` forceNetwork widget, or NULL when there is nothing to
#'   plot.
#' @keywords internal
#' @noRd
makeDrugFeatureNetwork <- function(top_data, bug, top_n = 10,
                                   include_clusters = FALSE,
                                  #  include_cogs = FALSE,
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

  df <- top_data |>
    dplyr::mutate(species = normalize_species(.data$species)) |>
    dplyr::filter(.data$species %in% normalize_species(bug)) |>
    dplyr::filter(is.na(.data$strat_label) | !nzchar(.data$strat_label)) |>
    dplyr::filter(!is.na(.data$Variable), !is.na(.data$drug_or_class)) |>
    dplyr::mutate(
      Variable = dplyr::case_when(
        feature_type == "domains" ~ sub("_.+$", "", Variable),
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

  # Take top N features per drug/class by max importance. Keep feature_type so
  # variable nodes can be coloured by molecular scale downstream.
  edges_df <- df |>
    dplyr::group_by(.data$drug_or_class, .data$Variable, .data$feature_type) |>
    dplyr::summarise(
      Importance = max(.data$Importance, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::group_by(.data$drug_or_class) |>
    dplyr::slice_max(order_by = .data$Importance, n = top_n) |>
    dplyr::ungroup()

  if (!nrow(edges_df)) {
    return(NULL)
  }

  drugs <- unique(edges_df$drug_or_class)
  variables <- unique(edges_df$Variable)
  # Molecular scale (feature_type) for each variable node, used for colouring.
  variables_scale <- edges_df$feature_type[match(variables, edges_df$Variable)]

  drug_var_edges <- edges_df |>
    dplyr::transmute(
      source_name = .data$drug_or_class,
      target_name = .data$Variable,
      value = .data$Importance
    )

  extra_nodes_cluster <- character(0)
  # extra_nodes_cog <- character(0)
  var_cluster_edges <- NULL
  # cluster_cog_edges <- NULL

  if (include_clusters) {
    ann <- load_feature_annotations(bug, results_root)
    if (!is.null(ann) && nrow(ann)) {
      # Match Variable ("PF23840_IPR056912") against ann$feature ("PF23840")
      # variables_key <- stringr::str_split_i(variables, "_", 1)
      # var_key_map <- stats::setNames(variables, variables_key)
      ann_sub <- ann |>
        dplyr::filter(.data$feature %in% variables)

      if (include_clusters && nrow(ann_sub)) {
        cl_edges <- ann_sub |>
          dplyr::filter(!is.na(.data$cluster)) |>
          dplyr::distinct(.data$feature, .data$cluster) |>
          dplyr::transmute(
            source_name = .data$feature,
            target_name = .data$cluster,
            value = 0.5
          )
        if (nrow(cl_edges)) {
          var_cluster_edges <- cl_edges
          extra_nodes_cluster <- unique(cl_edges$target_name)
        }
      }

      # if (include_cogs && nrow(ann_sub)) {
  #       if (include_clusters && !is.null(var_cluster_edges)) {
  #         # cluster -> COG edges
  #         cc_edges <- ann_sub |>
  #           dplyr::filter(!is.na(.data$cluster), !is.na(.data$COG)) |>
  #           dplyr::distinct(.data$cluster, .data$COG) |>
  #           dplyr::transmute(
  #             source_name = .data$cluster,
  #             target_name = .data$COG,
  #             value = 0.3
  #           )
  #         if (nrow(cc_edges)) {
  #           cluster_cog_edges <- cc_edges
  #           extra_nodes_cog <- unique(cc_edges$target_name)
  #         }
  #       } else {
  #         # variable -> COG directly
  #         cc_edges <- ann_sub |>
  #           dplyr::filter(!is.na(.data$COG)) |>
  #           dplyr::distinct(.data$feature, .data$COG) |>
  #           dplyr::transmute(
  #             source_name = var_key_map[.data$feature],
  #             target_name = .data$COG,
  #             value = 0.3
  #           )
  #         if (nrow(cc_edges)) {
  #           cluster_cog_edges <- cc_edges
  #           extra_nodes_cog <- unique(cc_edges$target_name)
  #         }
  #       }
  #     }
    }
  }

  # Variable nodes are grouped by their molecular scale (feature_type) so each
  # scale - COGs included - gets its own colour from the shared SCALE_COLORS
  # palette. Drug and cluster keep their own groups.
  nodes <- data.frame(
    name = c(drugs, variables, extra_nodes_cluster),
    node_type = c(
      rep("drug", length(drugs)),
      variables_scale,
      rep("cluster", length(extra_nodes_cluster))
    ),
    stringsAsFactors = FALSE
  )

  all_edges <- dplyr::bind_rows(
    drug_var_edges, var_cluster_edges
  )
  edges <- all_edges |>
    dplyr::transmute(
      source = match(.data$source_name, nodes$name) - 1,
      target = match(.data$target_name, nodes$name) - 1,
      value = .data$value
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
    NodeID = "name",
    Group = "node_type",
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
