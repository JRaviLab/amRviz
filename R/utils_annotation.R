# Feature/cluster/COG + drug-class annotation lookups and enrichment.


loadDrugClassMap <- function() {
  cwd <- getwd()
  # drug_class_map_fp <- file.path(cwd, "data", "drug_class_map.tsv")
  drug_class_map_fp <- system.file("extdata", "drug_class_map.tsv", package = "amRviz")
  message(stringr::str_glue("loadDrugClassMap(): Looking for TSV at: {drug_class_map_fp}"))
  drug_class_map_df <- readr::read_tsv(
    here(drug_class_map_fp),
    show_col_types = FALSE
  ) |>
    dplyr::select(drug.antibiotic_name, drug_class) |>
    dplyr::distinct()
  return(drug_class_map_df)
}


# Load a feature-id -> human-readable name mapping from the amRdata-style
# directory layout. Returns a tibble with columns `Variable` and `label`,
# or NULL if no matching parquet is found. Looks in:
#   {amrdata_root}/{species_dir}/{gene|domain|protein}_names.parquet
#   {results_root}/{species_dir}/{gene|domain|protein}_names.parquet
#   {extdata}/{species_dir}/{gene|domain|protein}_names.parquet
load_feature_name_map <- function(species_code, model_scale,
                                  amrdata_root = NULL,
                                  results_root = NULL) {
  scale <- dplyr::case_when(
    model_scale == "proteins" ~ "protein",
    model_scale == "domains" ~ "domain",
    model_scale == "genes" ~ "gene",
    TRUE ~ model_scale
  )
  fname <- paste0(scale, "_names.parquet")

  # Species-code and species-label maps use the same lookup paths as
  # cluster_feature_COG.parquet.
  roots <- c(
    amrdata_root, results_root,
    system.file("extdata", package = "amRviz")
  )
  roots <- roots[!is.null(roots) & nzchar(roots) & dir.exists(roots)]

  fp <- NULL
  for (r in roots) {
    for (d in list.dirs(r, full.names = TRUE, recursive = FALSE)) {
      cand <- file.path(d, fname)
      if (file.exists(cand)) {
        fp <- cand
        break
      }
    }
    if (!is.null(fp)) break
  }
  if (is.null(fp)) {
    return(NULL)
  }

  df <- tryCatch(arrow::read_parquet(fp), error = function(e) NULL)
  if (is.null(df) || !nrow(df)) {
    return(NULL)
  }

  # Normalise to {Variable, label}
  if (scale == "gene" && all(c("Gene", "Annotation") %in% names(df))) {
    return(tibble::tibble(
      Variable = df$Gene, label = df$Annotation
    ))
  }
  if (scale == "protein" && all(
    c("proteinID", "proteinName") %in% names(df)
  )) {
    return(tibble::tibble(
      Variable = df$proteinID, label = df$proteinName
    ))
  }
  if (scale == "domain" && all(c("DB.ID", "SignDesc") %in% names(df))) {
    # domain ids need deduping since one Pfam can occur many times
    agg <- df |>
      dplyr::distinct(.data$DB.ID, .data$SignDesc) |>
      dplyr::group_by(.data$DB.ID) |>
      dplyr::summarise(
        label = dplyr::first(.data$SignDesc),
        .groups = "drop"
      )
    return(tibble::tibble(Variable = agg$DB.ID, label = agg$label))
  }
  NULL
}


# Load cluster/COG annotations for a species if the parquet exists.
# Searches in results_root/<species_dir>/cluster_feature_COG.parquet first,
# then falls back to extdata/<species_dir>/cluster_feature_COG.parquet.
load_feature_annotations <- function(species_code, results_root = NULL) {
  fname <- "cluster_feature_COG.parquet"
  if (!is.null(results_root) && nzchar(results_root)) {
    for (d in list.dirs(results_root, full.names = TRUE, recursive = FALSE)) {
      fp <- file.path(d, fname)
      if (file.exists(fp)) {
        return(arrow::read_parquet(fp))
      }
    }
  }
  extdata <- system.file("extdata", package = "amRviz")
  if (nzchar(extdata)) {
    for (d in list.dirs(extdata, full.names = TRUE, recursive = FALSE)) {
      fp <- file.path(d, fname)
      if (file.exists(fp)) {
        return(arrow::read_parquet(fp))
      }
    }
  }
  NULL
}


# Enrich a top-features tibble with cluster/COG annotations joined on
# Variable -> feature. Collapses multiple COGs per feature into one comma-
# separated cell. Returns the input unchanged if no annotations found.
enrich_with_annotations <- function(tbl, species_code, results_root = NULL) {
  if (is.null(tbl) || !nrow(tbl) || !"Variable" %in% names(tbl)) {
    return(tbl)
  }
  ann <- load_feature_annotations(species_code, results_root)
  if (is.null(ann) || !nrow(ann)) {
    return(tbl)
  }

  # Top features Variable can be "PF23840_IPR056912" (domains) or the raw
  # feature id (genes/proteins). Split on first "_" to extract the join key.
  tbl <- tbl |>
    dplyr::mutate(
      .join_key = stringr::str_split_i(.data$Variable, "_", 1)
    )

  join_keys <- unique(tbl$.join_key)

  ann_collapsed <- ann |>
    dplyr::filter(.data$feature %in% join_keys) |>
    dplyr::group_by(.data$feature) |>
    dplyr::summarise(
      cluster = dplyr::first(.data$cluster),
      cluster_name = dplyr::first(.data$cluster_name),
      COG = paste(unique(stats::na.omit(.data$COG)), collapse = ", "),
      COG_name = paste(unique(stats::na.omit(.data$COG_name)),
        collapse = "; "
      ),
      .groups = "drop"
    )

  tbl |>
    dplyr::left_join(ann_collapsed, by = c(".join_key" = "feature")) |>
    dplyr::select(-".join_key")
}
