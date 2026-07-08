# Data loaders: ML results, top features, metadata, file discovery.


#' Normalise a results-root path
#'
#' @param results_root Single path string, or NULL.
#' @return The normalised absolute path, or NULL when `results_root` is missing,
#'   empty, or not a single string.
#' @keywords internal
#' @noRd
.normalize_results_root <- function(results_root) {
  if (is.null(results_root) || length(results_root) != 1 ||
    is.na(results_root) || !nzchar(results_root)) {
    return(NULL)
  }
  normalizePath(results_root, winslash = "/", mustWork = FALSE)
}


#' Read a parquet file, returning an empty tibble on failure
#'
#' @param fp Path to a parquet file.
#' @param verbose Whether to message on missing or unreadable files.
#' @return A tibble of the file contents, or an empty tibble when `fp` is
#'   missing or cannot be read.
#' @keywords internal
#' @noRd
.read_parquet_safe <- function(fp, verbose = TRUE) {
  if (is.null(fp) || !file.exists(fp)) {
    if (isTRUE(verbose)) message("File not found: ", fp)
    return(tibble::tibble())
  }
  tryCatch(
    arrow::read_parquet(fp),
    error = function(e) {
      if (isTRUE(verbose)) {
        message("Failed to read parquet: ", fp, " (", conditionMessage(e), ")")
      }
      tibble::tibble()
    }
  )
}


#' Discover species subdirectories that contain amRml output
#'
#' Files live in per-species subdirectories as
#' `{root}/{SpeciesDir}/{code}_perf.parquet`. Only directories holding at
#' least one baseline `all_perf.parquet` (not a country/year/MDR/cross variant)
#' are kept.
#'
#' @param results_root Root directory to scan.
#' @param verbose Unused; kept for signature consistency with the loaders.
#' @return A named character vector (names = directory basenames used as display
#'   labels, values = full species-subdirectory paths), or empty when nothing
#'   matches.
#' @keywords internal
#' @noRd
listAmRmlSpeciesFolders <- function(results_root, verbose = TRUE) {
  rr <- .normalize_results_root(results_root)
  if (is.null(rr) || !dir.exists(rr)) {
    return(character(0))
  }

  subdirs <- list.dirs(rr, full.names = TRUE, recursive = FALSE)
  if (!length(subdirs)) {
    return(character(0))
  }

  # Keep subdirs holding at least one baseline `all_perf.parquet` file, i.e. a
  # `*_perf.parquet` that is not a country/year/MDR/cross_drug variant.
  has_perf <- vapply(subdirs, function(d) {
    fps <- list.files(d, pattern = "_perf\\.parquet$", full.names = FALSE)
    any(!grepl("(country|year|MDR|cross_drug)_perf\\.parquet$", fps))
  }, logical(1))

  ok <- subdirs[has_perf]
  if (!length(ok)) {
    return(character(0))
  }
  setNames(ok, basename(ok))
}


#' Load all performance parquets from one species directory
#'
#' Reads baseline (`all_perf`), stratified (`country_perf`/`year_perf`), and
#' cross-tested (`cross_*_perf`) parquet files, excluding the MDR and
#' leave-one-out summaries. Cross-tested rows come in their own files in the
#' amRml schema, so each row is tagged with `cross_test` from its source
#' filename. Also tags `species_label = basename(species_dir)`.
#'
#' @param species_dir Full path to a species subdirectory.
#' @param verbose Passed through to .read_parquet_safe().
#' @return A tibble of combined performance rows, or an empty tibble.
#' @keywords internal
#' @noRd
.load_one_species_perf <- function(species_dir, verbose = TRUE) {
  fps <- list.files(
    species_dir,
    pattern = "_perf\\.parquet$", full.names = TRUE
  )
  fps <- fps[!grepl("(MDR|LOO)_perf\\.parquet$", basename(fps))]
  if (!length(fps)) {
    return(tibble::tibble())
  }
  df <- dplyr::bind_rows(lapply(fps, function(fp) {
    d <- .read_parquet_safe(fp, verbose = verbose)
    if (nrow(d)) d$cross_test <- grepl("^cross", basename(fp))
    d
  }))
  if (!nrow(df)) {
    return(df)
  }
  # Guarantee the columns the cross-model views read, even when a species ships
  # only baseline / self-evaluation files (no cross_* parquets).
  if (!"cross_test" %in% names(df)) df$cross_test <- FALSE
  df$cross_test[is.na(df$cross_test)] <- FALSE
  if (!"strat_value_test" %in% names(df)) df$strat_value_test <- NA
  df$species_label <- basename(species_dir)
  df
}


#' Load all top-feature parquets from one species directory
#'
#' Reads baseline (`all_top_features`), stratified (`country`/`year`), and
#' cross-tested (`cross_*`) `*_top_features.parquet` files, excluding the MDR
#' and leave-one-out summaries. Each row is tagged with `cross_test` from its
#' filename and with `species_label = basename(species_dir)`.
#'
#' @param species_dir Full path to a species subdirectory.
#' @param verbose Passed through to .read_parquet_safe().
#' @return A tibble of combined top-feature rows, or an empty tibble.
#' @keywords internal
#' @noRd
.load_one_species_top <- function(species_dir, verbose = TRUE) {
  fps <- list.files(
    species_dir,
    pattern = "_top_features\\.parquet$", full.names = TRUE
  )
  fps <- fps[!grepl("(MDR|LOO)_top_features\\.parquet$", basename(fps))]
  if (!length(fps)) {
    return(tibble::tibble())
  }
  df <- dplyr::bind_rows(lapply(fps, function(fp) {
    d <- .read_parquet_safe(fp, verbose = verbose)
    if (nrow(d)) d$cross_test <- grepl("^cross", basename(fp))
    d
  }))
  if (!nrow(df)) {
    return(df)
  }
  if (!"cross_test" %in% names(df)) df$cross_test <- FALSE
  df$cross_test[is.na(df$cross_test)] <- FALSE
  df$species_label <- basename(species_dir)
  df
}

#' Load all annotation parquets from one species directory
#'
#' Reads all `*_names.parquet`, `protein_*.parquet`
#' @param species_dir
#' @param verbose
#'
#' @return A tibble of combined annotation rows, or an empty tibble.
#' @keywords internal
#' @noRd
.load_one_species_annotations <- function(species_dir, verbose = TRUE) {
  ann_files <- list.files(
    species_dir,
    pattern = "_names\\.parquet$|^protein_.*\\.parquet$",
    full.names = TRUE
  )
  if (!length(ann_files)) {
    return(tibble::tibble())
  }
  # dplyr::bind_rows(lapply(ann_files, function(fp) {
  #   d <- .read_parquet_safe(fp, verbose = verbose)
  #   d$species_label <- basename(species_dir)
  #   d
  # }))

  dplyr::bind_rows(
    lapply(ann_files, function(fp) {
      d <- .read_parquet_safe(fp, verbose = verbose)

      ftype <- tools::file_path_sans_ext(basename(fp))

      d <- switch(ftype,
        gene_names = dplyr::transmute(
          d,
          Variable = Gene,
          description = Annotation,
          feature_type = "genes"
        ),
        protein_names = dplyr::transmute(
          d,
          Variable = proteinID,
          description = proteinName,
          feature_type = "proteins"
        ),
        domain_names = dplyr::transmute(
          d,
          Variable = DB.ID,
          description = SignDesc,
          feature_type = "domains"
        ),
        protein_COG = dplyr::transmute(
          d,
          Variable = name,
          description = description,
          feature_type = "cogs"
        ),
        protein_ResFinder = dplyr::transmute(
          d,
          Variable = name,
          description = description,
          feature_type = "args"
        ),
        NULL
      )

      if (!is.null(d)) {
        d$species_label <- basename(species_dir)
      }

      d
    })
  )
}

#' Load cluster feature parquet from one species directory
#'
#' @param species_dir
#' @param verbose
#'
#' @returns
#'
#' @keywords internal
#' @noRd
.load_one_species_cluster_features <- function(species_dir, verbose = TRUE) {
  cluster_file <- file.path(species_dir, "cluster_feature.parquet")
  cluster_annotation <- file.path(species_dir, "protein_names.parquet")
  if (!file.exists(cluster_file)) {
    if (isTRUE(verbose)) message("Cluster feature file not found: ", cluster_file)
    return(tibble::tibble())
  }
  d <- .read_parquet_safe(cluster_file, verbose = verbose)
  d$species_label <- basename(species_dir)
  if (file.exists(cluster_annotation)) {
    ann <- .read_parquet_safe(cluster_annotation, verbose = verbose)
    d <- dplyr::left_join(
      d,
      dplyr::transmute(
        ann,
        cluster = proteinID,
        cluster_name = proteinName
      ),
      by = c("cluster" = "cluster")
    )
  }
  d
}

#' Resolve and load multi-species amRml results
#'
#' Shared mode-selection used by loadMLResults() and loadTopFeat(): user mode
#' loads from the selected species subdirectories; if a results_root is set but
#' nothing is selected, returns empty; otherwise falls back to the packaged demo
#' parquets in extdata.
#'
#' @param loader Per-species loader, e.g. .load_one_species_perf().
#' @param results_root User results root, or NULL for demo mode.
#' @param species_dirs Selected species-subdirectory paths, or NULL.
#' @param verbose Passed through to `loader`.
#' @param demo_msg Message emitted when falling back to the demo parquets.
#' @return A tibble of combined results across the resolved directories.
#' @keywords internal
#' @noRd
.load_species_results <- function(loader, results_root, species_dirs,
                                  verbose, demo_msg) {
  rr <- .normalize_results_root(results_root)

  # User mode: results_root + selected species -> load from those subdirs
  if (!is.null(rr) && !is.null(species_dirs) && length(species_dirs) > 0) {
    return(dplyr::bind_rows(lapply(species_dirs, loader, verbose = verbose)))
  }

  # User mode: results_root set but nothing selected yet
  if (!is.null(rr) && is.null(species_dirs)) {
    return(tibble::tibble())
  }

  # Demo fallback: scan packaged extdata subdirectories
  extdata <- system.file("extdata", package = "amRviz")
  if (!nzchar(extdata)) {
    return(tibble::tibble())
  }
  subdirs <- list.dirs(extdata, full.names = TRUE, recursive = FALSE)
  if (isTRUE(verbose)) message(demo_msg)
  dplyr::bind_rows(lapply(subdirs, loader, verbose = verbose))
}


#' Load multi-species model performance results
#'
#' @param results_root User results root, or NULL for the packaged demo data.
#' @param species_dirs Full paths to the species subdirectories selected by the
#'   user, or NULL.
#' @param verbose Whether to emit progress/diagnostic messages.
#' @return A tibble of combined performance rows.
#' @keywords internal
#' @noRd
loadMLResults <- function(results_root = NULL, species_dirs = NULL,
                          verbose = TRUE) {
  .load_species_results(
    .load_one_species_perf, results_root, species_dirs, verbose,
    "loadMLResults(): using packaged demo parquets"
  )
}


#' Load multi-species top-feature results
#'
#' @param results_root User results root, or NULL for the packaged demo data.
#' @param species_dirs Full paths to the species subdirectories selected by the
#'   user, or NULL.
#' @param verbose Whether to emit progress/diagnostic messages.
#' @return A tibble of combined top-feature rows.
#' @keywords internal
#' @noRd
loadTopFeat <- function(results_root = NULL, species_dirs = NULL,
                        verbose = TRUE) {
  .load_species_results(
    .load_one_species_top, results_root, species_dirs, verbose,
    "loadTopFeat(): using packaged demo parquets"
  )
}


#' Find a file by name within the immediate subdirectories of a root
#'
#' @param root Directory whose immediate subdirectories are searched, or NULL.
#' @param fname File name to look for in each subdirectory.
#' @return The full path to the first match, or NULL when none is found.
#' @keywords internal
#' @noRd
.find_file_in_subdirs <- function(root, fname) {
  if (is.null(root) || !nzchar(root)) {
    return(NULL)
  }
  for (d in list.dirs(root, full.names = TRUE, recursive = FALSE)) {
    fp <- file.path(d, fname)
    if (file.exists(fp)) {
      return(fp)
    }
  }
  NULL
}


#' Locate a species metadata parquet
#'
#' Each species subdirectory holds a single `metadata.parquet`, so a species is
#' addressed by its subdirectory name. Looks for
#' `{results_root}/{species_dir}/metadata.parquet` (user mode) and then
#' `{extdata}/{species_dir}/metadata.parquet` (demo mode). When `species_dir` is
#' NULL, falls back to the first `metadata.parquet` found under any subdir.
#'
#' @param species_dir Species subdirectory basename (e.g. "Shigella_flexneri"),
#'   or NULL to return the first metadata parquet found.
#' @param results_root User results root, or NULL to search demo data only.
#' @return The full path to the metadata parquet, or NULL when not found.
#' @keywords internal
#' @noRd
get_metadata_path <- function(species_dir = NULL, results_root = NULL) {
  fname <- "metadata.parquet"

  if (!is.null(species_dir) && nzchar(species_dir)) {
    roots <- c(results_root, system.file("extdata", package = "amRviz"))
    for (root in roots) {
      if (is.null(root) || !nzchar(root)) next
      fp <- file.path(root, species_dir, fname)
      if (file.exists(fp)) {
        return(fp)
      }
    }
    return(NULL)
  }

  fp <- .find_file_in_subdirs(results_root, fname)
  if (!is.null(fp)) {
    return(fp)
  }
  .find_file_in_subdirs(system.file("extdata", package = "amRviz"), fname)
}
