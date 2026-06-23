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
#' `{root}/{SpeciesDir}/{code}_ML_perf.parquet`. Only directories holding at
#' least one baseline `*_ML_perf.parquet` (not a country/year/MDR/cross variant)
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

  # Keep subdirs with at least one baseline *_ML_perf.parquet file
  has_perf <- vapply(subdirs, function(d) {
    fps <- list.files(d, pattern = "_ML_perf\\.parquet$", full.names = FALSE)
    any(!grepl("_(country|year|MDR|cross)_ML_perf\\.parquet$", fps))
  }, logical(1))

  ok <- subdirs[has_perf]
  if (!length(ok)) {
    return(character(0))
  }
  setNames(ok, basename(ok))
}


#' Load all performance parquets from one species directory
#'
#' Reads baseline + country + year + cross `*_ML_perf.parquet` files (excluding
#' MDR) and tags rows with `species_label = basename(species_dir)`.
#'
#' @param species_dir Full path to a species subdirectory.
#' @param verbose Passed through to .read_parquet_safe().
#' @return A tibble of combined performance rows, or an empty tibble.
#' @keywords internal
#' @noRd
.load_one_species_perf <- function(species_dir, verbose = TRUE) {
  fps <- list.files(
    species_dir,
    pattern = "_ML_perf\\.parquet$", full.names = TRUE
  )
  fps <- fps[!grepl("_MDR_ML_perf\\.parquet$", fps)]
  if (!length(fps)) {
    return(tibble::tibble())
  }
  df <- dplyr::bind_rows(lapply(fps, .read_parquet_safe, verbose = verbose))
  if (nrow(df)) df$species_label <- basename(species_dir)
  df
}


#' Load all top-feature parquets from one species directory
#'
#' Reads baseline + country + year `*_ML_top_features.parquet` files (excluding
#' MDR) and tags rows with `species_label = basename(species_dir)`.
#'
#' @param species_dir Full path to a species subdirectory.
#' @param verbose Passed through to .read_parquet_safe().
#' @return A tibble of combined top-feature rows, or an empty tibble.
#' @keywords internal
#' @noRd
.load_one_species_top <- function(species_dir, verbose = TRUE) {
  fps <- list.files(
    species_dir,
    pattern = "_ML_top_features\\.parquet$", full.names = TRUE
  )
  fps <- fps[!grepl("_MDR_ML_top_features\\.parquet$", fps)]
  if (!length(fps)) {
    return(tibble::tibble())
  }
  df <- dplyr::bind_rows(lapply(fps, .read_parquet_safe, verbose = verbose))
  if (nrow(df)) df$species_label <- basename(species_dir)
  df
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
#' Searches the species subdirectories under `results_root` (user mode) and then
#' under the packaged extdata (demo mode) for `{species_code}_metadata.parquet`.
#'
#' @param species_code Species code prefix of the metadata file.
#' @param results_root User results root, or NULL to search demo data only.
#' @return The full path to the metadata parquet, or NULL when not found.
#' @keywords internal
#' @noRd
get_metadata_path <- function(species_code, results_root = NULL) {
  fname <- paste0(species_code, "_metadata.parquet")
  fp <- .find_file_in_subdirs(results_root, fname)
  if (!is.null(fp)) {
    return(fp)
  }
  .find_file_in_subdirs(system.file("extdata", package = "amRviz"), fname)
}
