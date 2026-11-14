# Force-directed network visualizations.
#
# The drug-feature network draws one table, built by networkTables(): its
# nodes and links. The ego network is drawn on the Features tab.


#' A feature's name as the network shows it
#'
#' Applied to the importance tables' names and to the ORB's own names alike, so
#' the two meet.
#'
#' @param feature_type Feature types as the tabs name them.
#' @param variable Feature names.
#' @return `variable`, tidied per feature type.
#' @keywords internal
#' @noRd
.network_feature_name <- function(feature_type, variable) {
  dplyr::case_when(
    # Older amRml named domains "PF23840_IPR056912"; keep the Pfam ID.
    # Current Pfam names ("DDE_Tnp_Tn3") are kept whole.
    feature_type == "domains" & grepl("^PF[0-9]+_", variable) ~
      sub("_.+$", "", variable),
    feature_type == "proteins" ~ sub("fig.", "fig|", variable, fixed = TRUE),
    feature_type == "args" ~ sub("^X", "", gsub("\\.NCBIFAM", "", variable)),
    TRUE ~ variable
  )
}


#' Nodes and links for the drug-feature network
#'
#' Links each drug or drug class to its `top_n` most important features of each
#' feature type in the baseline models (not stratified or cross-tested), and
#' optionally each
#' feature to the dyads it belongs to, from the ORB the models were fitted on.
#'
#' @param top_data Top-features tibble from loadTopFeat().
#' @param bug Species code.
#' @param top_n Number of top features per feature type, per drug or drug class.
#' @param include_dyads Add dyad nodes when TRUE and the ORB has dyads.
#' @return A list of `nodes` (`id`, unique across node types (a drug's starts
#'   with its unit, "drug" or "drug_class"), `label`, `group`:
#'   "drug", a feature type, or "dyad", `key`: the ID shown on hover, and `info`:
#'   the hover text as HTML) and `links` (`from` and
#'   `to` node ids, `importance`, and `sign`: "POS", "NEG" or NA for a link
#'   that has none), or NULL when there is nothing to draw.
#' @keywords internal
#' @noRd
networkTables <- function(top_data, bug, top_n = 10, include_dyads = FALSE) {
  if (is.null(top_data) || !is.data.frame(top_data) || !nrow(top_data)) {
    return(NULL)
  }

  df <- top_data |>
    dplyr::mutate(species = normalize_species(.data$species)) |>
    dplyr::filter(.data$species %in% normalize_species(bug)) |>
    dplyr::filter(is.na(.data$strat_label) | !nzchar(.data$strat_label)) |>
    dplyr::filter(!.data$cross_test %in% TRUE) |>
    dplyr::filter(!is.na(.data$Variable), !is.na(.data$drug_or_class)) |>
    dplyr::mutate(Variable = .network_feature_name(.data$feature_type, .data$Variable))
  if (!nrow(df)) {
    return(NULL)
  }

  # Top N features of each type per drug or class, by their highest importance
  # across seeds and encodings; the link takes the sign of that best model.
  if (!"Sign" %in% names(df)) {
    df$Sign <- NA_character_
  }
  # A drug and a drug class can share an abbreviation (TET), so a node is one
  # target of one unit, "drug" or "drug_class".
  if (!"drug_label" %in% names(df)) {
    df$drug_label <- "drug"
  }
  top <- df |>
    dplyr::group_by(
      .data$drug_label, .data$drug_or_class, .data$Variable, .data$feature_type
    ) |>
    dplyr::summarise(
      sign = .data$Sign[which.max(.data$Importance)],
      importance = max(.data$Importance, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::group_by(.data$drug_label, .data$drug_or_class, .data$feature_type) |>
    dplyr::slice_max(order_by = .data$importance, n = top_n, with_ties = FALSE) |>
    dplyr::ungroup()
  if (!nrow(top)) {
    return(NULL)
  }

  # A node's id includes its type, so a feature named like a drug, or the
  # same name in two feature types, stays a separate node.
  drug_id <- paste0(top$drug_label, ":", top$drug_or_class)
  # Say which one only where an abbreviation is both a drug and a class.
  both <- unique(top$drug_or_class[top$drug_label == "drug"]) |>
    intersect(top$drug_or_class[top$drug_label == "drug_class"])
  drug_label <- ifelse(
    top$drug_or_class %in% both,
    paste0(top$drug_or_class, ifelse(top$drug_label == "drug", " (drug)", " (class)")),
    top$drug_or_class
  )
  feature_id <- paste0(top$feature_type, ":", top$Variable)
  nodes <- dplyr::bind_rows(
    tibble::tibble(
      id = drug_id, label = drug_label, group = "drug",
      key = top$drug_or_class
    ),
    tibble::tibble(
      id = feature_id, label = top$Variable, group = top$feature_type,
      key = top$Variable
    )
  ) |>
    dplyr::distinct(.data$id, .keep_all = TRUE)
  links <- tibble::tibble(
    from = drug_id, to = feature_id, importance = top$importance, sign = top$sign
  )

  # The ORBs the models were fitted on hold the dyads and the feature names.
  orb_dirs <- if ("metadata_path" %in% names(df)) {
    unique(dirname(df$metadata_path[!is.na(df$metadata_path)]))
  } else {
    character(0)
  }

  if (include_dyads && length(orb_dirs)) {
    dyads <- load_dyad_links(orb_dirs)
    if (!is.null(dyads)) {
      dyads$feature <- .network_feature_name(dyads$feature_type, dyads$feature)
      hits <- dplyr::inner_join(
        dplyr::distinct(top, .data$feature_type, .data$Variable),
        dplyr::distinct(dyads, .data$feature_type, .data$feature, .data$dyad_id),
        by = c("feature_type", "Variable" = "feature")
      )
      if (nrow(hits)) {
        dyad_id <- paste0("dyad:", hits$dyad_id)
        nodes <- dplyr::bind_rows(nodes, tibble::tibble(
          id = unique(dyad_id),
          label = paste("Dyad", unique(hits$dyad_id)), group = "dyad",
          key = as.character(unique(hits$dyad_id))
        ))
        links <- dplyr::bind_rows(links, tibble::tibble(
          from = paste0(hits$feature_type, ":", hits$Variable),
          to = dyad_id, importance = 0.5, sign = NA_character_
        ))
      }
    }
  }

  metadata_paths <- if ("metadata_path" %in% names(df)) {
    unique(df$metadata_path)
  } else {
    character(0)
  }
  nodes$info <- network_node_info(nodes, orb_dirs, metadata_paths)
  list(nodes = nodes, links = links)
}


#' Interactive force-directed drug -> feature network
#'
#' Draws networkTables(): drugs (or drug classes) linked to their top
#' features, optionally extended to the dyads they belong to.
#'
#' @inheritParams networkTables
#' @return A `networkD3` forceNetwork widget, or NULL when there is nothing to
#'   plot.
#' @keywords internal
#' @noRd
makeDrugFeatureNetwork <- function(top_data, bug, top_n = 10,
                                   include_dyads = FALSE) {
  if (!requireNamespace("networkD3", quietly = TRUE)) {
    stop(
      "Package 'networkD3' is required. ",
      "Install it with install.packages('networkD3')."
    )
  }
  net <- networkTables(top_data, bug, top_n, include_dyads)
  if (is.null(net)) {
    return(NULL)
  }
  nodes <- as.data.frame(net$nodes)
  edges <- data.frame(
    source = match(net$links$from, nodes$id) - 1,
    target = match(net$links$to, nodes$id) - 1,
    value = net$links$importance
  )
  # A link is coloured by the sign of the feature's importance for the drug;
  # links to dyads have no sign.
  link_colour <- unname(.NETWORK_SIGN_COLORS[net$links$sign])
  link_colour[is.na(link_colour)] <- .NETWORK_NO_SIGN_COLOR

  # Colour domain: drug + every molecular scale + dyad, drawing scale colours
  # from the shared SCALE_COLORS palette so the network matches the other tabs.
  js_array <- function(x) paste0('["', paste(x, collapse = '","'), '"]')
  group_domain <- c("drug", names(SCALE_COLORS), "dyad")
  group_range <- c("#d4872a", unname(SCALE_COLORS), "#66a61e")
  color_scale <- sprintf(
    "d3.scaleOrdinal().domain(%s).range(%s)",
    js_array(group_domain), js_array(group_range)
  )

  .styleNetwork(info = nodes$info, networkD3::forceNetwork(
    Links = edges,
    Nodes = nodes,
    Source = "source",
    Target = "target",
    Value = "value",
    NodeID = "label",
    Group = "group",
    colourScale = networkD3::JS(color_scale),
    linkColour = link_colour,
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


# Link colours: the sign of a feature's importance for a drug (amRml's "POS" and
# "NEG"), and links that have no sign, such as those to dyads.
.NETWORK_SIGN_COLORS <- c(POS = "#c0392b", NEG = "#2c7bb6")
.NETWORK_NO_SIGN_COLOR <- "#999999"

# Font family for the force-directed networks, matching the plotly/ggplot
# panels (which use "Arial, sans-serif").
.NETWORK_FONT <- "Arial, sans-serif"

# Post-render hook for forceNetwork: replace each node's default circle with a
# shape keyed to its node_type group so element types are distinguishable
# beyond colour (drug = diamond, variable = circle, cluster/dyad = square,
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
      'dyad': d3.symbolSquare,
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

# Post-render hook for forceNetwork: a box that follows the pointer over a node
# and shows that node's `info` (HTML, one string per node, in node order).
.NETWORK_TOOLTIP_JS <- "
function(el, x, info) {
  try {
    var d3 = window.d3;
    if (!d3 || !info) { return; }
    info = [].concat(info);
    var root = d3.select(el);
    if (root.style('position') === 'static') { root.style('position', 'relative'); }
    var tip = root.append('div')
      .attr('class', 'network-tip')
      .style('position', 'absolute')
      .style('display', 'none')
      .style('pointer-events', 'none')
      .style('z-index', 10)
      .style('max-width', '320px')
      .style('padding', '6px 8px')
      .style('background', '#ffffff')
      .style('color', '#222222')
      .style('border', '1px solid #bbbbbb')
      .style('border-radius', '4px')
      .style('box-shadow', '0 2px 6px rgba(0,0,0,0.2)')
      .style('font', '12px Arial, sans-serif')
      .style('line-height', '1.4');
    function place() {
      var e = d3.event, box = el.getBoundingClientRect();
      var w = tip.node().offsetWidth, h = tip.node().offsetHeight;
      var left = e.clientX - box.left + 14, top = e.clientY - box.top + 14;
      if (left + w > box.width) { left = e.clientX - box.left - w - 14; }
      if (top + h > box.height) { top = e.clientY - box.top - h - 14; }
      tip.style('left', Math.max(0, left) + 'px').style('top', Math.max(0, top) + 'px');
    }
    // Namespaced listeners keep the widget's own hover behaviour.
    root.selectAll('.node')
      .on('mouseover.tip', function(d) {
        tip.html(info[d.index] || '').style('display', 'block');
        place();
      })
      .on('mousemove.tip', place)
      .on('mouseout.tip', function() { tip.style('display', 'none'); });
  } catch (e) { /* no hover box on error */ }
}
"

# Apply the shared font + per-group node shapes to a forceNetwork widget, and
# the hover box when `info` is given.
.styleNetwork <- function(net, info = NULL) {
  if (is.null(net)) {
    return(NULL)
  }
  net <- htmlwidgets::onRender(net, .NETWORK_SHAPE_JS)
  if (!is.null(info)) {
    net <- htmlwidgets::onRender(net, .NETWORK_TOOLTIP_JS, data = as.list(info))
  }
  net
}
