# Data loaders: ML results, top features, and finding them.
#
# A results source is one folder: an amRml run (read by utils_amrml.R) or an
# older amRml species folder (read by utils_legacy.R).


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


#' Display label for an amRml run
#'
#' @param run_dir An amRml run folder.
#' @return The ORB's species label and the run's start time (the folder name
#'   when the run has no run.json).
#' @keywords internal
#' @noRd
.amrml_run_label <- function(run_dir) {
  record <- .amrml_run_record(run_dir)
  when <- record$started_at
  if (is.na(when)) when <- basename(run_dir)
  paste0(record$species_label, " (amRml run ", when, ")")
}


#' Find the ML results under a root
#'
#' Finds amRml runs, inside the ORBs whose manifests sit in `results_root` or
#' up to two levels below it, or in `results_root` or its subfolders directly;
#' and older-format species folders among its subfolders.
#'
#' @param results_root Root directory to search.
#' @param verbose Unused; kept for signature consistency with the loaders.
#' @return A named character vector (names = display labels, values = result
#'   folders), or empty when nothing is found.
#' @keywords internal
#' @noRd
listResultSources <- function(results_root, verbose = TRUE) {
  rr <- .normalize_results_root(results_root)
  if (is.null(rr) || !dir.exists(rr)) {
    return(character(0))
  }

  subdirs <- list.dirs(rr, full.names = TRUE, recursive = FALSE)
  orbs <- .list_available_orbs(.find_orb_manifests(rr))
  runs <- unique(c(
    unlist(lapply(orbs$orb_dir, .orb_amrml_runs)),
    Filter(.is_amrml_run, c(rr, subdirs))
  ))
  legacy <- Filter(.is_legacy_results, setdiff(subdirs, runs))

  if (!length(runs) && !length(legacy)) {
    return(character(0))
  }
  labels <- c(
    vapply(runs, .amrml_run_label, character(1), USE.NAMES = FALSE),
    if (length(legacy)) {
      paste0(.legacy_label(legacy), " (older amRml results)")
    }
  )
  stats::setNames(c(runs, legacy), make.unique(labels, " #"))
}


#' Load one kind of result from one source
#'
#' @param dir An amRml run or older-format species folder.
#' @param kind `"perf"` or `"top"`.
#' @param verbose Passed through to the reader.
#' @keywords internal
#' @noRd
.load_source <- function(dir, kind, verbose = TRUE) {
  if (.is_amrml_run(dir)) {
    if (kind == "perf") {
      .load_amrml_perf(dir, verbose)
    } else {
      .load_amrml_top(dir, verbose)
    }
  } else {
    .load_legacy(dir, kind, verbose)
  }
}


#' Resolve and load results from several sources
#'
#' Shared mode-selection used by loadMLResults() and loadTopFeat(): user mode
#' loads from the selected result folders; if a results_root is set but nothing
#' is selected, returns empty; otherwise falls back to the packaged demo
#' parquets in extdata. Models fitted on shuffled labels (amRml's null
#' models) are left out; no tab shows them yet.
#'
#' @param kind `"perf"` or `"top"`.
#' @param results_root User results root, or NULL for demo mode.
#' @param result_dirs Selected result folders, or NULL.
#' @param verbose Passed through to the readers.
#' @param demo_msg Message emitted when falling back to the demo parquets.
#' @return A tibble of combined results across the resolved folders.
#' @keywords internal
#' @noRd
.load_results <- function(kind, results_root, result_dirs, verbose, demo_msg) {
  rr <- .normalize_results_root(results_root)

  # User mode: results_root + selected folders -> load from those
  if (!is.null(rr) && !is.null(result_dirs) && length(result_dirs) > 0) {
    return(.real_labels(dplyr::bind_rows(lapply(result_dirs, .load_source,
      kind = kind, verbose = verbose
    ))))
  }

  # User mode: results_root set but nothing selected yet
  if (!is.null(rr) && is.null(result_dirs)) {
    return(tibble::tibble())
  }

  # Demo fallback: the packaged older-format species folders
  extdata <- system.file("extdata", package = "amRviz")
  if (!nzchar(extdata)) {
    return(tibble::tibble())
  }
  subdirs <- list.dirs(extdata, full.names = TRUE, recursive = FALSE)
  if (isTRUE(verbose)) message(demo_msg)
  .real_labels(dplyr::bind_rows(lapply(subdirs, .load_legacy,
    kind = kind, verbose = verbose
  )))
}


#' Drop models fitted on shuffled labels
#'
#' @param df Result rows; older results may lack `shuffled` or leave it NA.
#' @return `df` without the rows whose `shuffled` is TRUE.
#' @keywords internal
#' @noRd
.real_labels <- function(df) {
  if (!"shuffled" %in% names(df)) {
    return(df)
  }
  df[!df$shuffled %in% TRUE, , drop = FALSE]
}


#' Load model performance results
#'
#' @param results_root User results root, or NULL for the packaged demo data.
#' @param result_dirs Result folders selected by the user (from
#'   listResultSources()), or NULL.
#' @param verbose Whether to emit progress/diagnostic messages.
#' @return A tibble of combined performance rows.
#' @keywords internal
#' @noRd
loadMLResults <- function(results_root = NULL, result_dirs = NULL,
                          verbose = TRUE) {
  .load_results(
    "perf", results_root, result_dirs, verbose,
    "loadMLResults(): using packaged demo parquets"
  )
}


#' Load top-feature results
#'
#' @param results_root User results root, or NULL for the packaged demo data.
#' @param result_dirs Result folders selected by the user (from
#'   listResultSources()), or NULL.
#' @param verbose Whether to emit progress/diagnostic messages.
#' @return A tibble of combined top-feature rows.
#' @keywords internal
#' @noRd
loadTopFeat <- function(results_root = NULL, result_dirs = NULL,
                        verbose = TRUE) {
  .load_results(
    "top", results_root, result_dirs, verbose,
    "loadTopFeat(): using packaged demo parquets"
  )
}


#' Find a file by name within the immediate subdirectories of a root
#'
#' Used only by the annotation lookups, which still read older-format files.
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
