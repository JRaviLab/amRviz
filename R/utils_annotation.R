# Feature/cluster/COG + drug-class annotation lookups and enrichment.


#' Load the drug -> drug-class lookup table
#'
#' Reads the packaged `drug_class_map.tsv` and returns the distinct
#' antibiotic-name / drug-class pairs.
#'
#' @return A tibble with columns `drug.antibiotic_name` and `drug_class`.
#' @keywords internal
#' @noRd
loadDrugClassMap <- function() {
  drug_class_map_fp <- system.file(
    "extdata", "drug_class_map.tsv",
    package = "amRviz"
  )
  message(stringr::str_glue(
    "loadDrugClassMap(): Looking for TSV at: {drug_class_map_fp}"
  ))
  readr::read_tsv(here(drug_class_map_fp), show_col_types = FALSE) |>
    dplyr::select(drug.antibiotic_name, drug_class) |>
    dplyr::distinct()
}


#' Load a feature-id -> human-readable name mapping
#'
#' Searches the amRdata-style layout
#' (`{root}/{species_dir}/{gene|domain|protein}_names.parquet`) across the
#' amrdata_root, results_root, and packaged extdata, then normalises the result
#' to `{Variable, label}`.
#'
#' @param species_code Species code (currently unused in the lookup; kept for
#'   call-site consistency).
#' @param model_scale Molecular scale ("genes", "proteins", or "domains").
#' @param amrdata_root Optional amRdata data root, searched first.
#' @param results_root Optional user results root, searched next.
#' @return A tibble with columns `Variable` and `label`, or NULL when no
#'   matching parquet is found.
#' @keywords internal
#' @noRd
load_feature_name_map <- function(model_scale,
                                  amrdata_root = NULL,
                                  results_root = NULL,
                                  species_dir = NULL) {
  scale <- dplyr::case_when(
    model_scale == "proteins" ~ "protein",
    model_scale == "domains" ~ "domain",
    model_scale == "genes" ~ "gene",
    TRUE ~ model_scale
  )
  # fname <- paste0(scale, "_names.parquet")

  # roots <- c(
  #   amrdata_root, results_root,
  #   system.file("extdata", package = "amRviz")
  # )
  # roots <- roots[!is.null(roots) & nzchar(roots) & dir.exists(roots)]

  # fp <- NULL
  # for (r in roots) {
  #   fp <- .find_file_in_subdirs(r, fname)
  #   if (!is.null(fp)) break
  # }
  # if (is.null(fp)) {
  #   return(NULL)
  # }

  # df <- tryCatch(arrow::read_parquet(fp), error = function(e) NULL)
  # if (is.null(df) || !nrow(df)) {
  #   return(NULL)
  # }

  # # Normalise to {Variable, label}
  # if (scale == "gene" && all(c("Gene", "Annotation") %in% names(df))) {
  #   return(tibble::tibble(Variable = df$Gene, label = df$Annotation))
  # }
  # if (scale == "protein" &&
  #   all(c("proteinID", "proteinName") %in% names(df))) {
  #   return(tibble::tibble(Variable = df$proteinID, label = df$proteinName))
  # }
  # if (scale == "domain" && all(c("DB.ID", "SignDesc") %in% names(df))) {
  #   # domain ids need deduping since one Pfam can occur many times
  #   agg <- df |>
  #     dplyr::distinct(.data$DB.ID, .data$SignDesc) |>
  #     dplyr::group_by(.data$DB.ID) |>
  #     dplyr::summarise(
  #       label = dplyr::first(.data$SignDesc),
  #       .groups = "drop"
  #     )
  #   return(tibble::tibble(Variable = agg$DB.ID, label = agg$label))
  # }
  # NULL
  df <- .load_one_species_annotations(species_dir) |>
    dplyr::filter(.data$feature_type == model_scale) |>
    dplyr::distinct()

  if (!nrow(df)) {
    return(NULL)
  }
  return(df)
}


#' Load cluster/COG annotations for a species
#'
#' Searches the species subdirectories under `results_root`, then the packaged
#' extdata, for `cluster_feature.parquet`.
#'
#' @param species_code Species code (currently unused; kept for call-site
#'   consistency).
#' @param results_root Optional user results root, searched first.
#' @return The annotations tibble, or NULL when no parquet is found.
#' @keywords internal
#' @noRd
load_feature_annotations <- function(species_code, results_root = NULL) {
  feature_fname <- "cluster_feature.parquet"
  protein_fname <- "protein_names.parquet"

  feature_fp <- .find_file_in_subdirs(results_root, feature_fname)
  if (is.null(feature_fp)) {
    feature_fp <- .find_file_in_subdirs(
      system.file("extdata", package = "amRviz"),
      feature_fname
    )
  }

  protein_fp <- .find_file_in_subdirs(results_root, protein_fname)
  if (is.null(protein_fp)) {
    protein_fp <- .find_file_in_subdirs(
      system.file("extdata", package = "amRviz"),
      protein_fname
    )
  }

  if (is.null(feature_fp) || is.null(protein_fp)) {
    return(NULL)
  }

  feature_tbl <- arrow::read_parquet(feature_fp)
  protein_tbl <- arrow::read_parquet(protein_fp)

  dplyr::left_join(
    feature_tbl,
    protein_tbl,
    by = c("cluster" = "proteinID") 
  ) |>
    dplyr::rename(cluster_name = "proteinName")
}


#' Enrich a top-features tibble with cluster annotations
#'
#' Joins annotations on the feature key extracted from `Variable` (the part
#' before the first "_"), collapsing multiple clusters per feature into a single
#' comma-separated cell. Returns `tbl` unchanged when no annotations are found.
#'
#' @param tbl Top-features tibble (must contain a `Variable` column).
#' @param species_dir Species directory passed to .load_one_species_cluster_features().
#' @param results_root Optional user results root.
#' @return `tbl` with cluster columns joined on, or unchanged when no
#'   annotations are available.
#' @keywords internal
#' @noRd
enrich_with_annotations <- function(tbl, species_dir, results_root = NULL) {
  if (is.null(tbl) || !nrow(tbl) || !"Variable" %in% names(tbl)) {
    return(tbl)
  }
  ann <- .load_one_species_cluster_features(species_dir = species_dir)
  if (is.null(ann) || !nrow(ann)) {
    return(tbl)
  }

  # Top features Variable can be "PF23840_IPR056912" (domains) or the raw
  # feature id (genes/proteins). Split on first "_" to extract the join key.
  # tbl <- tbl |>
  #   dplyr::mutate(
  #     .join_key = stringr::str_split_i(.data$Variable, "_", 1)
  #   )

  tbl <- tbl |>
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

  ann_collapsed <- ann |>
    # dplyr::filter(.data$feature %in% join_keys) |>
    dplyr::group_by(.data$feature) |>
    dplyr::summarise(
      cluster = dplyr::first(.data$cluster),
      cluster_name = dplyr::first(.data$cluster_name),
      # COG = paste(unique(stats::na.omit(.data$COG)), collapse = ", "),
      # COG_name = paste(unique(stats::na.omit(.data$COG_name)),
      #   collapse = "; "
      # ),
      .groups = "drop"
    )

  tbl |>
    dplyr::left_join(ann_collapsed, by = c("Variable" = "feature")) 
}
