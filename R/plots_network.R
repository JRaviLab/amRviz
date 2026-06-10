# Force-directed network visualisations.


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

  df <- top_data |>
    dplyr::mutate(species = normalize_species(.data$species)) |>
    dplyr::filter(.data$species %in% normalize_species(bug)) |>
    dplyr::filter(is.na(.data$strat_label) | !nzchar(.data$strat_label)) |>
    dplyr::filter(!is.na(.data$Variable), !is.na(.data$drug_or_class))

  if (!nrow(df)) {
    return(NULL)
  }

  # Take top N features per drug/class by max importance
  edges_df <- df |>
    dplyr::group_by(.data$drug_or_class, .data$Variable) |>
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

  drug_var_edges <- edges_df |>
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
      ann_sub <- ann |>
        dplyr::filter(.data$feature %in% variables_key)

      if (include_clusters && nrow(ann_sub)) {
        cl_edges <- ann_sub |>
          dplyr::filter(!is.na(.data$cluster)) |>
          dplyr::distinct(.data$feature, .data$cluster) |>
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
          cc_edges <- ann_sub |>
            dplyr::filter(!is.na(.data$cluster), !is.na(.data$COG)) |>
            dplyr::distinct(.data$cluster, .data$COG) |>
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
          cc_edges <- ann_sub |>
            dplyr::filter(!is.na(.data$COG)) |>
            dplyr::distinct(.data$feature, .data$COG) |>
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
  edges <- all_edges |>
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

  row <- enriched_tbl |>
    dplyr::filter(.data$Variable == variable) |>
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

