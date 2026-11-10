# Older amRml results: the per-species folders old amRml wrote
# (`all_perf.parquet`, `country_top_features.parquet`, ...), which the packaged
# demo data still uses. This file is the only place amRviz reads them, and goes
# once the demo data is rebuilt from an ORB.
#
# Each file is told apart by its columns, not its name: performance files have
# `mcc` and `model`, top-feature files `Variable` and `Importance`, and
# cross-tested files their test group (`test_drug` or `strat_value_test`). Old
# files don't record their species' full name, so the folder name is the one
# thing still used for it.


#' The older-format result files in a folder, and what each holds
#'
#' MDR files (no `species` column) and leave-one-out files
#' (`leaveout_strat_value`) are left out; no tab shows them yet.
#'
#' @param dir A species folder.
#' @return A tibble with `path`, `kind` (`"perf"` or `"top"`) and `cross_test`.
#' @keywords internal
#' @noRd
.legacy_result_files <- function(dir) {
  paths <- list.files(dir, pattern = "[.]parquet$", full.names = TRUE)
  columns <- lapply(paths, function(fp) {
    tryCatch(
      arrow::open_dataset(fp, format = "parquet")$schema$names,
      error = function(e) character(0)
    )
  })
  kind <- vapply(columns, function(cols) {
    if (!"species" %in% cols || "leaveout_strat_value" %in% cols) {
      NA_character_
    } else if (all(c("mcc", "model") %in% cols)) {
      "perf"
    } else if (all(c("Variable", "Importance") %in% cols)) {
      "top"
    } else {
      NA_character_
    }
  }, character(1))
  cross_test <- vapply(columns, function(cols) {
    any(c("test_drug", "strat_value_test") %in% cols)
  }, logical(1))

  files <- tibble::tibble(path = paths, kind = kind, cross_test = cross_test)
  files[!is.na(files$kind), , drop = FALSE]
}


#' Does a folder hold older-format results?
#' @keywords internal
#' @noRd
.is_legacy_results <- function(dir) {
  "perf" %in% .legacy_result_files(dir)$kind
}


#' Species label of an older-format folder: its name, with spaces
#' @keywords internal
#' @noRd
.legacy_label <- function(dir) {
  gsub("_", " ", basename(dir))
}


#' Load one kind of older-format result from a species folder
#'
#' Tags each row with `cross_test` from its file's columns, and with
#' `species_label` from the folder name and `metadata_path` from the folder's
#' `metadata.parquet`, the layout old amRml used.
#'
#' @param dir A species folder.
#' @param kind `"perf"` or `"top"`.
#' @param verbose Passed through to .read_parquet_safe().
#' @return A tibble of combined rows, or an empty tibble.
#' @keywords internal
#' @noRd
.load_legacy <- function(dir, kind, verbose = TRUE) {
  files <- .legacy_result_files(dir)
  files <- files[files$kind == kind, , drop = FALSE]
  df <- dplyr::bind_rows(lapply(seq_len(nrow(files)), function(i) {
    d <- .read_parquet_safe(files$path[i], verbose = verbose)
    if (nrow(d)) d$cross_test <- files$cross_test[i]
    d
  }))
  if (!nrow(df)) {
    return(df)
  }
  # Guarantee the columns the cross-model views read, even when a species ships
  # only baseline / self-evaluation files (no cross-tested files).
  if (!"strat_value_test" %in% names(df) && kind == "perf") {
    df$strat_value_test <- NA
  }
  metadata <- file.path(dir, "metadata.parquet")
  df$species_label <- .legacy_label(dir)
  df$metadata_path <- if (file.exists(metadata)) metadata else NA_character_
  df
}


#' Older-format metadata folders under a root
#'
#' @param root Folder whose subfolders are searched.
#' @return A named character vector (names = folder names with spaces,
#'   values = metadata paths).
#' @keywords internal
#' @noRd
.legacy_metadata_sources <- function(root) {
  if (is.null(root) || !nzchar(root) || !dir.exists(root)) {
    return(character(0))
  }
  dirs <- list.dirs(root, full.names = TRUE, recursive = FALSE)
  paths <- file.path(dirs, "metadata.parquet")
  keep <- file.exists(paths)
  stats::setNames(paths[keep], .legacy_label(dirs[keep]))
}
