# ORB discovery: find amRdata output bundles (ORBs) registered on this machine.
#
# amRdata records every dataset's manifest (`{dataset}/orb/manifest_run_*.json`)
# in the shared BiocFileCache under the resource name `amR_dataset_manifest_*`.
# That registry, not a directory scan, is how amRdata itself finds ORBs (see
# `.amr_resolve_export_orb()` in amRdata), so the Metadata tab does the same.
# An ORB is addressed in the UI by an `orb:<orb_dir>` id so that it stays
# unique even when two ORBs share a species name.


#' Empty ORB listing
#' @keywords internal
#' @noRd
.empty_orbs <- function() {
  tibble::tibble(
    label = character(),
    species_dir = character(),
    dataset_id = character(),
    orb_dir = character(),
    manifest_path = character(),
    metadata_path = character(),
    modified = as.POSIXct(character())
  )
}


#' Manifest paths registered in the amRdata BiocFileCache
#'
#' @return A character vector of manifest paths (possibly stale or missing on
#'   disk), or `character(0)` when the cache is empty or cannot be opened.
#' @keywords internal
#' @noRd
.registered_orb_manifests <- function() {
  tryCatch(
    {
      bfc <- BiocFileCache::BiocFileCache(ask = FALSE)
      hits <- BiocFileCache::bfcquery(
        bfc,
        query = "^amR_dataset_manifest_",
        field = "rname",
        exact = FALSE
      )
      as.character(hits$rpath)
    },
    error = function(e) character(0)
  )
}


#' Describe one ORB from its manifest
#'
#' Applies the same checks amRdata does (parseable `amR_dataset` manifest,
#' schema version 1, living in an `orb/` folder), plus that `metadata.parquet`
#' exists, since that is all the Metadata tab needs.
#'
#' @param manifest_path Path to a `manifest_run_*.json`.
#' @return A one-row tibble in the shape of `.empty_orbs()`, or NULL when the
#'   manifest is missing, invalid, or its ORB has no `metadata.parquet`.
#' @keywords internal
#' @noRd
.describe_orb <- function(manifest_path) {
  if (is.null(manifest_path) || length(manifest_path) != 1 ||
    is.na(manifest_path) || !nzchar(manifest_path) ||
    !file.exists(manifest_path)) {
    return(NULL)
  }

  manifest <- tryCatch(
    jsonlite::read_json(manifest_path, simplifyVector = FALSE),
    error = function(e) NULL
  )
  if (!is.list(manifest) ||
    !identical(manifest$manifest_type, "amR_dataset") ||
    !identical(suppressWarnings(as.integer(manifest$schema_version)), 1L)) {
    return(NULL)
  }

  manifest_path <- normalizePath(manifest_path, winslash = "/", mustWork = TRUE)
  orb_dir <- dirname(manifest_path)
  metadata_path <- file.path(orb_dir, "metadata.parquet")
  if (!identical(basename(orb_dir), "orb") || !file.exists(metadata_path)) {
    return(NULL)
  }

  species_dir <- basename(dirname(orb_dir))
  user_bacs <- as.character(unlist(
    manifest$dataset$selection$user_bacs,
    use.names = FALSE
  ))
  label <- if (length(user_bacs)) {
    paste(user_bacs, collapse = ", ")
  } else {
    gsub("_", " ", species_dir)
  }

  tibble::tibble(
    label = label,
    species_dir = species_dir,
    dataset_id = as.character(manifest$dataset_id %||% ""),
    orb_dir = orb_dir,
    manifest_path = manifest_path,
    metadata_path = metadata_path,
    modified = file.info(manifest_path)$mtime
  )
}


#' List the ORBs available on this machine
#'
#' Stale registry entries (deleted files, temp-dir runs, pre-`orb/` layouts)
#' are dropped silently. When several manifests point at the same ORB, the
#' newest wins.
#'
#' @param manifest_paths Manifest paths to inspect. Defaults to those
#'   registered in the amRdata BiocFileCache.
#' @return A tibble with columns `label`, `species_dir`, `dataset_id`,
#'   `orb_dir`, `manifest_path`, `metadata_path` and `modified`, newest first.
#' @keywords internal
#' @noRd
.list_available_orbs <- function(manifest_paths = .registered_orb_manifests()) {
  orbs <- dplyr::bind_rows(lapply(manifest_paths, .describe_orb))
  if (!nrow(orbs)) {
    return(.empty_orbs())
  }
  orbs |>
    dplyr::arrange(dplyr::desc(.data$modified)) |>
    dplyr::distinct(.data$orb_dir, .keep_all = TRUE)
}


#' Selector id for an ORB
#' @keywords internal
#' @noRd
.orb_id <- function(orb_dir) paste0("orb:", orb_dir)


#' Is `x` an ORB selector id?
#' @keywords internal
#' @noRd
.is_orb_id <- function(x) {
  is.character(x) && length(x) == 1 && !is.na(x) && startsWith(x, "orb:")
}


#' Selector choices for a set of ORBs
#'
#' Names are the display labels; labels shared by several ORBs get the
#' manifest's modification time appended (and a counter if even that ties) so
#' they can be told apart.
#'
#' @param orbs Output of `.list_available_orbs()`.
#' @return A named character vector of ORB ids, sorted by label.
#' @keywords internal
#' @noRd
.orb_choices <- function(orbs) {
  if (is.null(orbs) || !nrow(orbs)) {
    return(character(0))
  }
  labels <- orbs$label
  dup <- labels %in% labels[duplicated(labels)]
  labels[dup] <- paste0(
    labels[dup], " [", format(orbs$modified[dup], "%Y-%m-%d %H:%M"), "]"
  )
  choices <- stats::setNames(.orb_id(orbs$orb_dir), make.unique(labels, " #"))
  choices[order(names(choices))]
}


#' Metadata parquet behind an ORB selector id
#'
#' @param id An ORB selector id (`orb:<orb_dir>`).
#' @return The path to the ORB's `metadata.parquet`, or NULL when `id` is not an
#'   ORB id or the file is missing.
#' @keywords internal
#' @noRd
.orb_metadata_path <- function(id) {
  if (!.is_orb_id(id)) {
    return(NULL)
  }
  fp <- file.path(sub("^orb:", "", id), "metadata.parquet")
  if (file.exists(fp)) fp else NULL
}


#' Species folder name behind a metadata selector id
#'
#' Plain species folders are their own id; an ORB id maps to the dataset
#' folder that holds its `orb/` directory.
#'
#' @keywords internal
#' @noRd
.metadata_species_name <- function(id) {
  if (.is_orb_id(id)) {
    basename(dirname(sub("^orb:", "", id)))
  } else {
    id
  }
}


#' Does a species folder hold a metadata parquet, flat or under `orb/`?
#' @keywords internal
#' @noRd
.species_dir_has_metadata <- function(d) {
  any(file.exists(c(
    file.path(d, "metadata.parquet"),
    file.path(d, "orb", "metadata.parquet")
  )))
}
