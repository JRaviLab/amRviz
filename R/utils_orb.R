# ORB discovery: find amRdata output bundles (ORBs) and describe them from
# their manifests.
#
# amRdata records every dataset's manifest (`manifest_run_*.json`) in the shared
# BiocFileCache under the resource name `amR_dataset_manifest_*`, and writes it
# next to the ORB's files. Everything amRviz knows about an ORB comes from that
# manifest: its species (`dataset_id`, `dataset.selection.user_bacs`) and where
# its metadata is (`artifacts.amRml_input.metadata_parquet`). Folder names are
# never read for meaning, the same rule amRml's `readORB()` follows.


# Names and formats set by amRdata, kept in one place.
# The manifest type and schema version this code can read.
.AMRDATA_MANIFEST_TYPE <- "amR_dataset"
.AMRDATA_SCHEMA_VERSION <- 1L
# The manifest entry that lists the files amRml reads.
.AMRDATA_ML_ARTIFACT <- "amRml_input"
# File names amRdata gives its manifests.
.MANIFEST_PATTERN <- "^manifest_.*[.]json$"


#' Empty ORB listing
#' @keywords internal
#' @noRd
.empty_orbs <- function() {
  tibble::tibble(
    label = character(),
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


#' Manifests in a folder and the two levels below it
#'
#' Covers a folder that is an ORB, one that holds ORBs, and amRdata's data
#' folder (`{data}/{dataset}/orb/`).
#'
#' @param root Folder to search.
#' @return A character vector of manifest paths.
#' @keywords internal
#' @noRd
.find_orb_manifests <- function(root) {
  if (is.null(root) || !nzchar(root) || !dir.exists(root)) {
    return(character(0))
  }
  level1 <- list.dirs(root, full.names = TRUE, recursive = FALSE)
  level2 <- unlist(lapply(level1, list.dirs,
    full.names = TRUE, recursive = FALSE
  ))
  unlist(lapply(c(root, level1, level2), list.files,
    pattern = .MANIFEST_PATTERN, full.names = TRUE
  ))
}


#' A path a manifest records, resolved against the manifest's folder
#'
#' amRdata writes paths relative to the manifest; older ORBs and ORBs written
#' on Windows may record absolute or backslashed paths, so only the file name
#' is kept, as amRml does.
#'
#' @keywords internal
#' @noRd
.orb_file_path <- function(path, manifest_path) {
  file.path(dirname(manifest_path), basename(chartr("\\", "/", path)))
}


#' Describe one ORB from its manifest
#'
#' Applies the same checks amRml's `readORB()` does on the manifest (parseable
#' `amR_dataset` manifest, schema version 1, an `amRml_input` entry naming the
#' metadata file), plus that the metadata file exists, since every tab needs it.
#'
#' @param manifest_path Path to a `manifest_run_*.json`.
#' @return A one-row tibble in the shape of `.empty_orbs()`, or NULL when the
#'   manifest is missing or invalid, or names no metadata file that exists.
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
    !identical(manifest[["manifest_type"]], .AMRDATA_MANIFEST_TYPE) ||
    !identical(
      suppressWarnings(as.integer(manifest[["schema_version"]])),
      .AMRDATA_SCHEMA_VERSION
    )) {
    return(NULL)
  }

  manifest_path <- normalizePath(manifest_path, winslash = "/", mustWork = TRUE)
  artifact <- manifest[["artifacts"]][[.AMRDATA_ML_ARTIFACT]]
  metadata <- artifact[["metadata_parquet"]]
  if (!is.character(metadata) || length(metadata) != 1 || !nzchar(metadata)) {
    return(NULL)
  }
  metadata_path <- .orb_file_path(metadata, manifest_path)
  if (!file.exists(metadata_path)) {
    return(NULL)
  }

  dataset_id <- as.character(manifest[["dataset_id"]] %||% "")
  user_bacs <- as.character(unlist(
    manifest[["dataset"]][["selection"]][["user_bacs"]],
    use.names = FALSE
  ))
  label <- if (length(user_bacs)) {
    paste(user_bacs, collapse = ", ")
  } else {
    dataset_id
  }
  if (!nzchar(label)) {
    return(NULL)
  }

  tibble::tibble(
    label = label,
    dataset_id = dataset_id,
    orb_dir = dirname(manifest_path),
    manifest_path = manifest_path,
    metadata_path = metadata_path,
    modified = file.info(manifest_path)$mtime
  )
}


#' List ORBs from their manifests
#'
#' Stale registry entries (deleted files, temp-dir runs, manifests without a
#' metadata file) are dropped silently. When several manifests describe the
#' same ORB, the newest wins.
#'
#' @param manifest_paths Manifest paths to inspect. Defaults to those
#'   registered in the amRdata BiocFileCache.
#' @return A tibble with columns `label`, `dataset_id`, `orb_dir`,
#'   `manifest_path`, `metadata_path` and `modified`, newest first.
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


#' Selector choices for a set of ORBs
#'
#' Names are the display labels; labels shared by several ORBs get the
#' manifest's modification time appended (and a counter if even that ties) so
#' they can be told apart.
#'
#' @param orbs Output of `.list_available_orbs()`.
#' @return A named character vector of metadata paths, sorted by label.
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
  choices <- stats::setNames(orbs$metadata_path, make.unique(labels, " #"))
  choices[order(names(choices))]
}
