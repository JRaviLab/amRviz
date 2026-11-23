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
#' @param view Draw drug-to-feature or direct drug-to-dyad associations.
#'   Default keeps the existing feature view.
#' @param show_features Also show contributing features in dyad view.
#' @param rescoring Optional list from loadRescoring(); NULL uses raw features.
#' @param thresholds Rescoring filters, as from rescoringDefaults().
#' @param feature_types Feature types to include; NULL keeps all, character(0)
#'   keeps none. In dyad view, selects only the extra feature nodes;
#'   drug-dyad associations and support counts are unchanged.
#' @return A list of `nodes` (`id`, unique across node types (a drug's starts
#'   with its unit, "drug" or "drug_class"), `label`, `group`:
#'   "drug", a feature type, or "dyad", `key`: the ID shown on hover, and `info`:
#'   the hover text as HTML) and `links` (`from` and
#'   `to` node ids, `importance`, and `sign`: "POS", "NEG" or NA for a link
#'   that has none), or NULL when there is nothing to draw.
#' @keywords internal
#' @noRd
networkTables <- function(top_data, bug, top_n = 10, include_dyads = FALSE,
                          view = c("feature", "dyad"), show_features = FALSE,
                          rescoring = NULL, thresholds = rescoringDefaults(),
                          feature_types = NULL) {
  view <- match.arg(view)
  if (!is.null(rescoring)) {
    return(networkTablesRescored(rescoring, bug, top_n, include_dyads,
                                 thresholds, view = view, show_features = show_features,
                                 feature_types = feature_types))
  }
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
  if (view == "feature" && !is.null(feature_types)) df <- dplyr::filter(df, .data$feature_type %in% feature_types)
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

  orb_dirs <- if ("metadata_path" %in% names(df)) {
    unique(dirname(df$metadata_path[!is.na(df$metadata_path)]))
  } else {
    character(0)
  }
  metadata_paths <- if ("metadata_path" %in% names(df)) {
    unique(df$metadata_path)
  } else {
    character(0)
  }

  # Dyads: each shown feature to the dyads it belongs to, in the ORB.
  dyad_hits <- NULL
  if ((include_dyads || view == "dyad") && length(orb_dirs)) {
    dyads <- load_dyad_links(orb_dirs)
    if (!is.null(dyads)) {
      dyads$feature <- .network_feature_name(dyads$feature_type, dyads$feature)
      dyad_hits <- dplyr::inner_join(
        dplyr::distinct(top, .data$feature_type, .data$Variable),
        dplyr::distinct(dyads, .data$feature_type, .data$feature, .data$dyad_id),
        by = c("feature_type", "Variable" = "feature")
      ) |>
        dplyr::transmute(
          .data$feature_type, .data$Variable,
          dyad_key = as.character(.data$dyad_id),
          dyad_label = paste("Dyad", .data$dyad_id)
        )
    }
  }
  net <- .network_assemble(top, dyad_hits, orb_dirs, metadata_paths)
  if (view == "dyad") .network_dyad_view(net, top, dyad_hits, show_features,
                                        feature_types = feature_types) else net
}


#' Nodes and links from the features chosen for each drug
#'
#' @param top One row per drug or class and feature: `drug_label`,
#'   `drug_or_class`, `feature_type`, `Variable`, `importance` and `sign`.
#' @param dyad_hits Optional dyads to add: `feature_type`, `Variable` (a
#'   feature in `top`), `dyad_key` (the dyad's ID, as shown on hover) and
#'   `dyad_label`.
#' @param orb_dirs,metadata_paths The ORBs the models were fitted on and their
#'   metadata, for the hover text.
#' @param scores Optional rescoring scores to add to the hover text; see
#'   network_node_info().
#' @return The list described in networkTables().
#' @keywords internal
#' @noRd
.network_assemble <- function(top, dyad_hits, orb_dirs, metadata_paths,
                              scores = NULL) {
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

  if (!is.null(dyad_hits) && nrow(dyad_hits)) {
    dyad_hits <- dplyr::distinct(
      dyad_hits, .data$feature_type, .data$Variable, .data$dyad_key,
      .data$dyad_label
    )
    dyad_id <- paste0("dyad:", dyad_hits$dyad_key)
    keep <- !duplicated(dyad_id)
    nodes <- dplyr::bind_rows(nodes, tibble::tibble(
      id = dyad_id[keep], label = dyad_hits$dyad_label[keep], group = "dyad",
      key = dyad_hits$dyad_key[keep]
    ))
    links <- dplyr::bind_rows(links, tibble::tibble(
      from = paste0(dyad_hits$feature_type, ":", dyad_hits$Variable),
      to = dyad_id, importance = 0.5, sign = NA_character_
    ))
  }

  nodes$info <- network_node_info(nodes, orb_dirs, metadata_paths, scores)
  list(nodes = nodes, links = links)
}


#' Nodes and links for the network of rescored features
#'
#' Like networkTables(), but the features and dyads come from amRml's
#' rescoring of the run (see loadRescoring()): the features that pass the
#' thresholds, up to `top_n` of each type per drug or class by rank score, and
#' optionally the rescored dyads that reach them. The hover text also gives the
#' scores: each drug's MCC, and each feature's and dyad's scores per drug.
#'
#' @param rescoring A list of `features` and `dyads` from loadRescoring().
#' @param bug Species code.
#' @param top_n Number of top features per feature type, per drug or drug class.
#' @param include_dyads Add dyad nodes when TRUE.
#' @param thresholds A list as from rescoringDefaults().
#' @param view Network view: features or direct drug-to-dyad associations.
#' @param show_features Add the selected features through their dyads.
#' @param feature_types Feature types to include, as for networkTables().
#' @return The list described in networkTables(), or NULL when nothing passes.
#' @keywords internal
#' @noRd
networkTablesRescored <- function(rescoring, bug, top_n = 10,
                                  include_dyads = FALSE,
                                  thresholds = rescoringDefaults(),
                                  view = c("feature", "dyad"), show_features = FALSE,
                                  feature_types = NULL) {
  view <- match.arg(view)
  # Only features from model fits that pass the performance thresholds.
  models <- rescoredModels(rescoring$models, bug, thresholds)
  features <- rescoredFeatures(
    rescoring$features, bug, thresholds,
    models = if (!is.null(rescoring$models) && nrow(rescoring$models)) models
  )
  if (view == "feature" && !is.null(feature_types)) {
    features <- dplyr::filter(features, .data$feature_type %in% feature_types)
    if (nrow(models)) models <- dplyr::filter(models, .data$feature_type %in% feature_types)
  }
  if (!nrow(features)) {
    return(NULL)
  }
  # Variable names as the network shows them, then the top N of each type.
  features$Variable <- .network_feature_name(features$feature_type, features$Variable)
  top <- features |>
    dplyr::group_by(.data$drug_label, .data$drug_or_class, .data$feature_type) |>
    dplyr::slice_max(
      order_by = .data$rank_score * 10 + .data$contribution,
      n = top_n, with_ties = FALSE
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      # Link width follows the feature's share of the model's importance.
      importance = pmax(.data$contribution * 10, 0.25)
    )

  orb_dirs <- unique(dirname(top$metadata_path[!is.na(top$metadata_path)]))
  metadata_paths <- unique(top$metadata_path)

  scores <- list(
    drugs = rescoredDrugScores(models),
    features = top[c(
      "drug_label", "drug_or_class", "feature_type", "Variable", "rank_score",
      "contribution", "seed_ratio"
    )]
  )
  dyad_hits <- NULL
  if (include_dyads || view == "dyad") {
    reached <- rescoredDyads(rescoring$dyads, top, thresholds)
    if (nrow(reached)) {
      # Name a dyad by its ID in the ORB where it is there, else as rescored.
      orb <- dplyr::bind_rows(lapply(orb_dirs, .read_dyads))
      known <- if (nrow(orb)) {
        stats::setNames(orb$dyad_id, paste(orb$protein, orb$gene, sep = "|"))
      } else {
        integer(0)
      }
      # A fused gene is "a...b" in the rescoring and "a~~~b" in the ORB.
      name <- gsub("...", "~~~", reached$dyad, fixed = TRUE)
      id <- unname(known[name])
      dyad_hits <- tibble::tibble(
        drug_label = reached$drug_label, drug_or_class = reached$drug_or_class,
        feature_type = reached$feature_type, Variable = reached$Variable,
        dyad_key = ifelse(is.na(id), name, as.character(id)),
        dyad_label = paste("Dyad", ifelse(is.na(id), name, id))
      )
      scores$dyads <- dplyr::distinct(tibble::tibble(
        drug_label = reached$drug_label, drug_or_class = reached$drug_or_class,
        dyad_key = dyad_hits$dyad_key, rank_score = reached$dyad_rank_score,
        contribution = reached$dyad_contribution, frequency = reached$frequency
      ))
    }
  }
  net <- .network_assemble(top, dyad_hits, orb_dirs, metadata_paths, scores)
  if (view == "dyad") {
    .network_dyad_view(net, top, dyad_hits, show_features, scores$dyads, feature_types)
  } else net
}


# Show a simplified graph of drug -> dyad + optionally -> features
.network_dyad_view <- function(net, top, hits, show_features, scores = NULL,
                               feature_types = NULL) {
  if (is.null(hits) || !nrow(hits)) return(NULL)
  keys <- c("drug_label", "drug_or_class", "feature_type", "Variable")
  if (all(c("drug_label", "drug_or_class") %in% names(hits))) {
    hits <- dplyr::inner_join(hits, top[c(keys, "importance", "sign")], by = keys)
  } else {
    hits <- dplyr::inner_join(
      hits, top[c(keys, "importance", "sign")], by = c("feature_type", "Variable")
    )
  }
  # Counting R/S association.
  ### One vote per typed feature, regardless of duplicates
  hits <- dplyr::distinct(hits, .data$drug_label, .data$drug_or_class,
    .data$dyad_key, .data$feature_type, .data$Variable, .keep_all = TRUE)
  associations <- hits |>
    dplyr::group_by(.data$drug_label, .data$drug_or_class, .data$dyad_key) |>
    dplyr::summarise(
      features = paste(sort(unique(paste0(.data$feature_type, ":", .data$Variable))), collapse = ", "),
      n_r = sum(.data$sign %in% "NEG"),
      n_s = sum(.data$sign %in% "POS"),
      n_mixed = sum(.data$sign %in% "MIXED"),
      n_unknown = sum(!.data$sign %in% c("NEG", "POS", "MIXED")),
      importance = max(.data$importance), .groups = "drop"
    ) |>
    dplyr::mutate(sign = dplyr::case_when(
      .data$n_r > .data$n_s ~ "NEG",
      .data$n_s > .data$n_r ~ "POS",
      TRUE ~ NA_character_
    ))
  if (!is.null(scores)) {
    associations <- dplyr::left_join(
      associations, scores[c("drug_label", "drug_or_class", "dyad_key", "contribution")],
      by = c("drug_label", "drug_or_class", "dyad_key")
    ) |>
      dplyr::mutate(importance = pmax(.data$contribution * 10, 0.25))
  }
  net$links <- associations |>
    dplyr::transmute(
      from = paste0(.data$drug_label, ":", .data$drug_or_class),
      to = paste0("dyad:", .data$dyad_key), .data$importance,
      .data$sign, .data$features, .data$n_r, .data$n_s,
      .data$n_mixed, .data$n_unknown
    ) |>
    dplyr::distinct()
  if (show_features) {
    shown <- if (is.null(feature_types)) hits else dplyr::filter(hits, .data$feature_type %in% feature_types)
    net$links <- dplyr::bind_rows(net$links, shown |>
      dplyr::transmute(
        from = paste0("dyad:", .data$dyad_key),
        to = paste0(.data$feature_type, ":", .data$Variable),
        importance = 0.5, sign = NA_character_, features = NA_character_
      ) |>
      dplyr::distinct())
  }
  net$nodes <- net$nodes[net$nodes$id %in% c(net$links$from, net$links$to), ]
  for (i in which(net$nodes$group == "dyad")) {
    a <- associations[associations$dyad_key == net$nodes$key[i], ]
    target <- paste0(a$drug_or_class, ifelse(a$drug_label == "drug_class", " (class)", ""))
    lines <- htmltools::htmlEscape(paste0(target, ": ", a$features))
    support <- vapply(seq_len(nrow(a)), function(j) {
      counts <- c(R = a$n_r[j], S = a$n_s[j], mixed = a$n_mixed[j], unknown = a$n_unknown[j])
      counts <- counts[counts > 0]
      paste0(target[j], ": ", paste(counts, names(counts), collapse = ", "))
    }, character(1))
    support <- htmltools::htmlEscape(support)
    net$nodes$info[i] <- paste0(net$nodes$info[i],
      "<br><b>Support by drug:</b><br>", paste(support, collapse = "<br>"),
      "<br><b>Features by drug:</b><br>", paste(lines, collapse = "<br>"))
  }
  net
}


#' Interactive force-directed drug -> feature network
#'
#' Draws networkTables(): drugs (or drug classes) linked to their top
#' features, optionally extended to the dyads they belong to.
#'
#' @inheritParams networkTables
#' @param rescoring A list from loadRescoring(); when given, the network draws
#'   the rescored features and dyads instead of `top_data`.
#' @param thresholds Thresholds for the rescoring, as from rescoringDefaults().
#' @param view,show_features As for networkTables().
#' @param feature_types Feature types to include, as for networkTables().
#' @param weighted_edges Vary edge thickness by importance/contribution;
#'   FALSE draws every edge at the same width without changing the scores.
#' @return A `networkD3` forceNetwork widget, or NULL when there is nothing to
#'   plot.
#' @keywords internal
#' @noRd
makeDrugFeatureNetwork <- function(top_data, bug, top_n = 10,
                                   include_dyads = FALSE, rescoring = NULL,
                                   thresholds = rescoringDefaults(),
                                   view = c("feature", "dyad"), show_features = FALSE,
                                   feature_types = NULL, weighted_edges = TRUE) {
  if (!requireNamespace("networkD3", quietly = TRUE)) {
    stop(
      "Package 'networkD3' is required. ",
      "Install it with install.packages('networkD3')."
    )
  }
  view <- match.arg(view)
  net <- networkTables(top_data, bug, top_n, include_dyads,
                       view = view, show_features = show_features,
                       rescoring = rescoring, thresholds = thresholds,
                       feature_types = feature_types)
  if (is.null(net)) {
    return(NULL)
  }
  nodes <- as.data.frame(net$nodes)
  # The widget's legend lists these group names, so use the names shown.
  nodes$group <- .network_group_label(nodes$group)
  edges <- data.frame(
    source = match(net$links$from, nodes$id) - 1,
    target = match(net$links$to, nodes$id) - 1,
    value = net$links$importance
  )
  # Feature edges use their sign, direct dyad edge uses R/S majority in case of
  # split evidence. Membership edges and dyad ties have no sign
  sign_colors <- .network_sign_colors()
  link_colour <- unname(sign_colors[net$links$sign])
  link_colour[is.na(link_colour)] <- .NETWORK_NO_SIGN_COLOR
  # One legend entry per link type
  legend_rows <- list(
    POS = c(if (view == "dyad") "Majority towards susceptibility" else "Towards susceptibility", sign_colors[["POS"]]),
    NEG = c(if (view == "dyad") "Majority towards resistance" else "Towards resistance", sign_colors[["NEG"]]),
    none = c(if (view == "dyad") "No majority / NA" else "Feature in a dyad",
             .NETWORK_NO_SIGN_COLOR)
  )
  present <- c(
    POS = "POS" %in% net$links$sign, NEG = "NEG" %in% net$links$sign,
    none = any(is.na(net$links$sign))
  )
  edge_legend <- unname(lapply(legend_rows[present], function(r) {
    list(label = r[[1]], colour = r[[2]])
  }))

  # Colour domain: drug + every molecular scale + dyad, drawing scale colours
  # from the shared SCALE_COLORS palette so the network matches the other tabs.
  js_array <- function(x) paste0('["', paste(x, collapse = '","'), '"]')
  group_domain <- .network_group_label(c("drug", names(SCALE_COLORS), "dyad"))
  group_range <- c(.NETWORK_DRUG_COLOR, unname(SCALE_COLORS), .NETWORK_DYAD_COLOR)
  color_scale <- sprintf(
    "d3.scaleOrdinal().domain(%s).range(%s)",
    js_array(group_domain), js_array(group_range)
  )

  .styleNetwork(info = nodes$info, edge_legend = edge_legend, networkD3::forceNetwork(
    Links = edges,
    Nodes = nodes,
    Source = "source",
    Target = "target",
    Value = "value",
    NodeID = "label",
    Group = "group",
    colourScale = networkD3::JS(color_scale),
    linkColour = link_colour,
    linkWidth = networkD3::JS(if (weighted_edges) {
      "function(d) { return Math.sqrt(d.value); }"
    } else {
      "function(d) { return 1.5; }"
    }),
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


# Link colours: the sign of a feature's effect for a drug (amRml's "POS" and
# "NEG"), and links that have no sign, such as those to dyads. A negative sign
# points toward resistance (blaTEM for ampicillin, tet genes for tetracycline),
# so the links take the Metadata tab's colours for Resistant and Susceptible.
# A function, since the palette is defined in another file.
.network_sign_colors <- function() {
  c(
    POS = PHENOTYPE_COLORS[["Susceptible"]],
    NEG = PHENOTYPE_COLORS[["Resistant"]]
  )
}
.NETWORK_NO_SIGN_COLOR <- "#999999"

# Dyad nodes, and the dyad sliders to match.
.NETWORK_DYAD_COLOR <- "#66a61e"

# How the widget's legend names the node groups, where it isn't the group's
# own name.
.NETWORK_GROUP_LABELS <- c(
  drug = "drugs/drug classes", args = "AMRFinder", cogs = "COGs", DefenseCas = "Defense/Cas"
)

#' Legend names for node groups
#'
#' @param group Node groups: "drug", a feature type or "dyad".
#' @return The names the legend shows.
#' @keywords internal
#' @noRd
.network_group_label <- function(group) {
  named <- group %in% names(.NETWORK_GROUP_LABELS)
  group[named] <- unname(.NETWORK_GROUP_LABELS[group[named]])
  group
}

# Drug nodes: neutral, so they don't read as either phenotype colour.
.NETWORK_DRUG_COLOR <- "#3d3d3d"

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
    var svg = d3.select(el).select('svg');
    svg.property('__zoom', d3.zoomIdentity);
    svg.select('.zoom-layer').attr('transform', d3.zoomIdentity);
    var shapeMap = {
      'drug': d3.symbolDiamond,
      'drugs/drug classes': d3.symbolDiamond,
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
    // A re-render keeps the widget's div, so drop the previous box.
    root.selectAll('.network-tip').remove();
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

# Post-render hook for forceNetwork: a legend for the link colours, under the
# widget's own node legend. `items` is a list of `label` and `colour`.
.NETWORK_EDGE_LEGEND_JS <- "
function(el, x, items) {
  try {
    var d3 = window.d3;
    if (!d3 || !items) { return; }
    items = [].concat(items);
    var svg = d3.select(el).select('svg');
    // Start below the node legend: 22px per entry, as networkD3 spaces it.
    var top = svg.selectAll('.legend').size() * 22 + 18;
    var g = svg.append('g').attr('class', 'edge-legend')
      .attr('transform', 'translate(18,' + top + ')');
    items.forEach(function(d, i) {
      var row = g.append('g').attr('transform', 'translate(0,' + (i * 22) + ')');
      row.append('line')
        .attr('x1', 0).attr('x2', 18).attr('y1', 9).attr('y2', 9)
        .style('stroke', d.colour).style('stroke-width', 3);
      row.append('text')
        .attr('x', 22).attr('y', 14)
        .style('font-size', '12px').style('fill', '#333333')
        .style('font-family', x.options.fontFamily)
        .text(d.label);
    });
  } catch (e) { /* no link legend on error */ }
}
"

# Apply the shared font + per-group node shapes to a forceNetwork widget, and
# the hover box when `info` is given.
.styleNetwork <- function(net, info = NULL, edge_legend = NULL) {
  if (is.null(net)) {
    return(NULL)
  }
  net <- htmlwidgets::onRender(net, .NETWORK_SHAPE_JS)
  if (!is.null(info)) {
    net <- htmlwidgets::onRender(net, .NETWORK_TOOLTIP_JS, data = as.list(info))
  }
  if (length(edge_legend)) {
    net <- htmlwidgets::onRender(net, .NETWORK_EDGE_LEGEND_JS, data = edge_legend)
  }
  net
}
