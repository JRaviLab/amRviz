# Feature importance visualisations.


#' Rank wide importance rows and convert to a labelled matrix
#'
#' Orders rows by coverage (number of non-NA groups) then peak importance, drops
#' the scratch ranking columns, and returns a feature x group numeric matrix.
#'
#' @param vi_wider Wide importance tibble with a `COG_name` column.
#' @param group_cols Names of the group (value) columns.
#' @return A numeric matrix with features as row names and groups as columns.
#' @keywords internal
#' @noRd
.vi_matrix <- function(vi_wider, group_cols) {
  vi_wider |>
    dplyr::rowwise() |>
    dplyr::mutate(
      .n = sum(!is.na(c_across(all_of(group_cols)))),
      .mx = max(c_across(all_of(group_cols)), na.rm = TRUE)
    ) |>
    dplyr::ungroup() |>
    dplyr::arrange(dplyr::desc(.data$.n), dplyr::desc(.data$.mx)) |>
    dplyr::select(-c(".n", ".mx")) |>
    tibble::column_to_rownames("COG_name") |>
    as.matrix()
}


#' Feature-importance heatmap across species or drugs
#'
#' Heatmap of top features across bugs (across_bug) or drugs (across_drug),
#' min-max normalised within each group. Feature labels come from an annotation
#' join when available, otherwise the raw Variable id. Uses the amRml columns
#' drug_or_class, feature_subtype (encoding), feature_type (scale), and
#' strat_label (NA for baseline).
#'
#' @param data Top-features tibble from loadTopFeat() / topFeatures().
#' @param bug Species code(s) to include.
#' @param amr_drug Drug/class identifier(s) to include.
#' @param model_scale Molecular scale (feature_type) to include.
#' @param data_type_ Data encoding(s) (feature_subtype) to include.
#' @param top_n_features Number of features per group, or "all".
#' @param feature_importance_tabset "across_bug" or "across_drug".
#' @param annotated_dir Optional directory of annotated parquet files.
#' @param amrdata_root,results_root Optional roots for the name-map lookup.
#' @return A plotly heatmap, or NULL when there is nothing to plot.
#' @keywords internal
#' @noRd
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
  top_features_df <- data |>
    dplyr::mutate(species = normalize_species(.data$species)) |>
    dplyr::filter(.data$species %in% bug_norm) |>
    dplyr::filter(.data$drug_or_class %in% amr_drug) |>
    dplyr::filter(.data$feature_type %in% c(model_scale, "struct")) |>
    dplyr::filter(.data$feature_subtype %in% data_type_) |>
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
        top_features_df |>
          dplyr::inner_join(annotated_table, by = join_by_expr) |>
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
    top_features_df <- top_features_df |>
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
  top_features_df <- top_features_df |>
    dplyr::group_by(!!rlang::sym(group_column), .data$COG_name) |>
    dplyr::summarize(Importance = max(.data$Importance, na.rm = TRUE), .groups = "drop")

  # Min-max normalise within each group
  top_features_df <- top_features_df |>
    dplyr::group_by(!!rlang::sym(group_column)) |>
    dplyr::mutate(
      Importance = (.data$Importance - min(.data$Importance, na.rm = TRUE)) /
        (max(.data$Importance, na.rm = TRUE) - min(.data$Importance, na.rm = TRUE))
    ) |>
    dplyr::ungroup()

  # Slice top N features per group
  if (!identical(top_n_features, "all")) {
    top_features_df <- top_features_df |>
      dplyr::group_by(!!rlang::sym(group_column)) |>
      dplyr::slice_max(order_by = .data$Importance, n = top_n_features) |>
      dplyr::ungroup()
  }

  if (!nrow(top_features_df)) {
    return(NULL)
  }

  # Build wide matrix
  if (feature_importance_tabset == "across_bug") {
    vi_wider <- top_features_df |>
      dplyr::select("COG_name", "Importance", "species") |>
      dplyr::distinct() |>
      tidyr::pivot_wider(names_from = "species", values_from = "Importance")

    group_cols <- setdiff(colnames(vi_wider), "COG_name")
    if (!length(group_cols)) {
      return(NULL)
    }

    vi_mat <- .vi_matrix(vi_wider, group_cols)

    eskape_order <- c("Efa", "Sau", "Kpn", "Aba", "Pae", "Esp")
    col_ord <- intersect(eskape_order, colnames(vi_mat))
    if (length(col_ord)) vi_mat <- vi_mat[, col_ord, drop = FALSE]
  }

  if (feature_importance_tabset == "across_drug") {
    top_features_df <- top_features_df |>
      dplyr::mutate(
        drug_or_class = stringr::str_trim(as.character(.data$drug_or_class))
      )

    vi_wider <- top_features_df |>
      dplyr::select("COG_name", "Importance", "drug_or_class") |>
      dplyr::group_by(.data$COG_name, .data$drug_or_class) |>
      dplyr::summarise(
        Importance = max(.data$Importance, na.rm = TRUE), .groups = "drop"
      ) |>
      tidyr::pivot_wider(
        names_from = "drug_or_class", values_from = "Importance"
      )

    group_cols <- setdiff(colnames(vi_wider), "COG_name")
    if (!length(group_cols)) {
      return(NULL)
    }

    vi_mat <- .vi_matrix(vi_wider, group_cols)
  }

  if (!exists("vi_mat") || !length(vi_mat)) {
    return(NULL)
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
  ) |>
    plotly::layout(
      xaxis = list(title = "", tickangle = -45, side = "top"),
      yaxis = list(title = "", autorange = "reversed"),
      margin = list(l = 200, b = 20, t = 80)
    )
}


#' Drop literal "NA" tokens from a COG_name string
#'
#' The annotation source records unnamed COGs as the literal string "NA" and
#' space-joins them into `COG_name` (e.g. "Alanine racemase NA NA NA"). This
#' removes those standalone tokens, returning `NA` when nothing meaningful is
#' left so callers can fall back to the bare COG id.
#'
#' @param x Character vector of `COG_name` values.
#' @return A character vector with `NA` tokens removed; `NA_character_` where no
#'   real name remains.
#' @keywords internal
#' @noRd
.clean_cog_name <- function(x) {
  vapply(
    x,
    function(s) {
      if (is.na(s)) {
        return(NA_character_)
      }
      toks <- setdiff(strsplit(s, "[[:space:]]+")[[1]], c("NA", ""))
      if (!length(toks)) NA_character_ else paste(toks, collapse = " ")
    },
    character(1),
    USE.NAMES = FALSE
  )
}


#' Horizontal bar chart of the most common COGs
#'
#' Counts COG occurrences across the features in an annotation-enriched
#' top-features tibble and shows the `top_n` most frequent.
#'
#' @param enriched_tbl Annotation-enriched top-features tibble (needs a `COG`
#'   column; see enrich_with_annotations()).
#' @param top_n Number of COGs to display.
#' @return A horizontal plotly bar chart (empty placeholder when there are no
#'   annotations).
#' @keywords internal
#' @noRd
makeCogBarChart <- function(enriched_tbl, top_n = 15) {
  if (is.null(enriched_tbl) || !nrow(enriched_tbl) ||
    !"COG" %in% names(enriched_tbl)) {
    return(plotly_placeholder("No annotations available"))
  }

  # Split the comma-separated COG cells and count occurrences per Variable.
  cog_df <- enriched_tbl |>
    dplyr::filter(!is.na(.data$COG), nzchar(.data$COG)) |>
    dplyr::select("Variable", "COG", dplyr::any_of("COG_name")) |>
    dplyr::distinct()

  if (!nrow(cog_df)) {
    return(plotly_placeholder("No COGs in selection"))
  }

  rows <- do.call(rbind, lapply(seq_len(nrow(cog_df)), function(i) {
    cogs <- trimws(strsplit(cog_df$COG[i], ",", fixed = TRUE)[[1]])
    # The annotation source stores unnamed COGs as the literal token "NA" and
    # space-joins them into COG_name (e.g. "Alanine racemase NA NA NA");
    # .clean_cog_name() strips those so they don't leak into the bar labels.
    clean <- if ("COG_name" %in% names(cog_df)) {
      .clean_cog_name(cog_df$COG_name[i])
    } else {
      NA_character_
    }
    names <- if (!is.na(clean)) {
      n <- trimws(strsplit(clean, ";", fixed = TRUE)[[1]])
      rep_len(n, length(cogs))
    } else {
      rep(NA_character_, length(cogs))
    }
    data.frame(COG = cogs, COG_name = names, stringsAsFactors = FALSE)
  }))

  counts <- rows |>
    dplyr::filter(nzchar(.data$COG)) |>
    dplyr::count(.data$COG, .data$COG_name, name = "n") |>
    dplyr::arrange(dplyr::desc(.data$n)) |>
    dplyr::slice_head(n = top_n) |>
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
  ) |>
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


#' Render comma-separated ids as HTML links
#'
#' @param ids A single string of comma-separated ids (or NA).
#' @param make_url Function mapping one id to its URL.
#' @return The ids rendered as comma-separated `<a>` links, or `ids` unchanged
#'   when empty/NA.
#' @keywords internal
#' @noRd
.link_ids <- function(ids, make_url) {
  if (is.na(ids) || !nzchar(ids)) {
    return(ids)
  }
  parts <- trimws(strsplit(ids, ",", fixed = TRUE)[[1]])
  linked <- vapply(parts, function(id) {
    paste0(
      "<a href='", make_url(id), "' target='_blank' ",
      "style='color:#1a73e8; text-decoration: underline;'>",
      id, "</a>"
    )
  }, character(1))
  paste(linked, collapse = ", ")
}


#' Build an interactive feature-importance table
#'
#' Formats numeric columns, reorders to a preferred column order, and renders
#' clickable links (accession -> NCBI, cluster -> BV-BRC, COG -> NCBI COG) as a
#' DT datatable.
#'
#' @param feature_import_table Feature-importance tibble to display.
#' @return A `DT::datatable` HTML widget.
#' @keywords internal
#' @noRd
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
  # Link cluster (fig IDs) to BV-BRC, one <a> per unique id (comma-sep cells).
  if ("cluster" %in% names(tbl)) {
    tbl$cluster <- vapply(tbl$cluster, function(ids) {
      .link_ids(ids, function(id) {
        paste0(
          "https://www.bv-brc.org/view/Feature/",
          utils::URLencode(id, reserved = TRUE)
        )
      })
    }, character(1))
  }
  # Link COG ids (comma-separated) to the NCBI COG page.
  if ("COG" %in% names(tbl)) {
    tbl$COG <- vapply(tbl$COG, function(ids) {
      .link_ids(ids, function(id) {
        paste0("https://www.ncbi.nlm.nih.gov/research/cog/cog/", id)
      })
    }, character(1))
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
