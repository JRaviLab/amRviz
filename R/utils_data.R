# Data loaders: ML results, top features, metadata, file discovery.


# helpers to discover species folders and load to combine results generated from amRml
.normalize_results_root <- function(results_root) {
  if (is.null(results_root) || length(results_root) != 1 || is.na(results_root) || !nzchar(results_root)) {
    return(NULL)
  }
  normalizePath(results_root, winslash = "/", mustWork = FALSE)
}


.read_parquet_safe <- function(fp, verbose = TRUE) {
  if (is.null(fp) || !file.exists(fp)) {
    if (isTRUE(verbose)) message("File not found: ", fp)
    return(tibble::tibble())
  }
  tryCatch(
    arrow::read_parquet(fp),
    error = function(e) {
      if (isTRUE(verbose)) message("Failed to read parquet: ", fp, " (", conditionMessage(e), ")")
      tibble::tibble()
    }
  )
}


# Discover species subdirectories under results_root that contain amRml output.
# Files live inside per-species subdirectories: {root}/{SpeciesDir}/{code}_ML_perf.parquet
# Returns a named character vector: names = directory basename (display label),
# values = full path to the species subdirectory.
listAmRmlSpeciesFolders <- function(results_root, verbose = TRUE) {
  rr <- .normalize_results_root(results_root)
  if (is.null(rr) || !dir.exists(rr)) {
    return(character(0))
  }

  subdirs <- list.dirs(rr, full.names = TRUE, recursive = FALSE)
  if (!length(subdirs)) {
    return(character(0))
  }

  # Keep subdirs that contain at least one baseline *_ML_perf.parquet file
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


# Load all performance parquets (baseline + country + year + cross) from a species directory.
# species_dir: full path to the species subdirectory (e.g. "/results/Shigella_flexneri")
# Attaches species_label = basename(species_dir) so the full name is available for display.
.load_one_species_perf <- function(species_dir, verbose = TRUE) {
  fps <- list.files(species_dir, pattern = "_ML_perf\\.parquet$", full.names = TRUE)
  fps <- fps[!grepl("_MDR_ML_perf\\.parquet$", fps)]
  if (!length(fps)) {
    return(tibble::tibble())
  }
  df <- dplyr::bind_rows(lapply(fps, .read_parquet_safe, verbose = verbose))
  if (nrow(df)) df$species_label <- basename(species_dir)
  df
}


# Load all top-feature parquets (baseline + country + year) from a species directory.
.load_one_species_top <- function(species_dir, verbose = TRUE) {
  fps <- list.files(species_dir, pattern = "_ML_top_features\\.parquet$", full.names = TRUE)
  fps <- fps[!grepl("_MDR_ML_top_features\\.parquet$", fps)]
  if (!length(fps)) {
    return(tibble::tibble())
  }
  df <- dplyr::bind_rows(lapply(fps, .read_parquet_safe, verbose = verbose))
  if (nrow(df)) df$species_label <- basename(species_dir)
  df
}


# Public loaders used by app.R (multi-species aware).
# species_dirs: character vector of full paths to species subdirectories selected by the user.
loadMLResults <- function(results_root = NULL, species_dirs = NULL, verbose = TRUE) {
  rr <- .normalize_results_root(results_root)

  # User mode: results_root + species selected -> load from selected subdirectories
  if (!is.null(rr) && !is.null(species_dirs) && length(species_dirs) > 0) {
    dfs <- lapply(species_dirs, .load_one_species_perf, verbose = verbose)
    return(dplyr::bind_rows(dfs))
  }

  # User mode: results_root provided but nothing selected yet
  if (!is.null(rr) && is.null(species_dirs)) {
    return(tibble::tibble())
  }

  # Demo fallback: scan extdata subdirectories recursively for *_ML_perf.parquet
  extdata <- system.file("extdata", package = "amRviz")
  if (!nzchar(extdata)) {
    return(tibble::tibble())
  }
  subdirs <- list.dirs(extdata, full.names = TRUE, recursive = FALSE)
  if (isTRUE(verbose)) message("loadMLResults(): using packaged demo parquets")
  dplyr::bind_rows(lapply(subdirs, .load_one_species_perf, verbose = verbose))
}


loadTopFeat <- function(results_root = NULL, species_dirs = NULL, verbose = TRUE) {
  rr <- .normalize_results_root(results_root)

  if (!is.null(rr) && !is.null(species_dirs) && length(species_dirs) > 0) {
    dfs <- lapply(species_dirs, .load_one_species_top, verbose = verbose)
    return(dplyr::bind_rows(dfs))
  }

  if (!is.null(rr) && is.null(species_dirs)) {
    return(tibble::tibble())
  }

  # Demo fallback: scan extdata subdirectories for *_ML_top_features.parquet
  extdata <- system.file("extdata", package = "amRviz")
  if (!nzchar(extdata)) {
    return(tibble::tibble())
  }
  subdirs <- list.dirs(extdata, full.names = TRUE, recursive = FALSE)
  if (isTRUE(verbose)) message("loadTopFeat(): using packaged demo parquets")
  dplyr::bind_rows(lapply(subdirs, .load_one_species_top, verbose = verbose))
}


# Return path to a species metadata parquet.
# Searches inside the species subdirectory under results_root (user mode) or
# extdata (demo mode). The metadata file follows the pattern {code}_metadata.parquet
# and lives alongside the other parquets for that species.
get_metadata_path <- function(species_code, results_root = NULL) {
  fname <- paste0(species_code, "_metadata.parquet")

  # User mode: search species subdirectories under results_root
  if (!is.null(results_root) && nzchar(results_root)) {
    subdirs <- list.dirs(results_root, full.names = TRUE, recursive = FALSE)
    for (d in subdirs) {
      fp <- file.path(d, fname)
      if (file.exists(fp)) {
        return(fp)
      }
    }
  }

  # Demo mode: search species subdirectories under extdata
  extdata <- system.file("extdata", package = "amRviz")
  if (nzchar(extdata)) {
    subdirs <- list.dirs(extdata, full.names = TRUE, recursive = FALSE)
    for (d in subdirs) {
      fp <- file.path(d, fname)
      if (file.exists(fp)) {
        return(fp)
      }
    }
  }
  NULL
}
