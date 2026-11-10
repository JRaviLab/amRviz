# Reading amRml runs: the results of one amRml run on one amRdata ORB.
#
# amRml writes each run inside its ORB, at `{orb}/amRml/{run}/`, with a
# `manifest.json` whose `amRviz_input` artifact lists where the run's tables
# are. amRviz finds runs only there: the ORB is the folder two levels up, and
# its manifest gives the run's species and metadata.
#
# Everything about a model comes from the run's own tables, never from file or
# folder names: each result row's `task_id` joins to the tasks table (labels
# and seed), whose `matrix_id` joins to the matrices table (mode, drug or
# class, feature type, encoding, groups).
#
# The tabs still read old amRml's column names, so rows are returned under
# those, with the same mapping as amRml's temporary `exportMergedResults()`.


# A run's manifest: its file name, its type, and the artifact listing the
# tables amRviz reads.
.AMRML_MANIFEST <- "manifest.json"
.AMRML_MANIFEST_TYPE <- "amR_ml"
.AMRML_INPUT <- "amRviz_input"

# Where each table amRviz reads is listed in the amRviz_input artifact, as
# (section, entry).
.AMRML_TABLES <- list(
  tasks = c("results", "tasks"),
  matrices = c("metadata", "matrices"),
  performance = c("results", "performance"),
  importance = c("results", "importance")
)

# Features kept per model: the most any tab shows (the Features slider's
# maximum). Large runs hold tens of millions of importance rows.
.AMRML_TOP_N <- 100

# amRml's modeling modes, as `amRml::modelingModes()` lists them: what each
# splits genomes by, and how it tests.
.AMRML_MODES <- tibble::tribble(
  ~mode_id,             ~grouping, ~LOO,  ~cross_test,
  "vanilla",            "none",    FALSE, FALSE,
  "stratified_year",    "year",    FALSE, FALSE,
  "stratified_country", "country", FALSE, FALSE,
  "cross_year",         "year",    FALSE, TRUE,
  "cross_country",      "country", FALSE, TRUE,
  "cross_drug",         "drug",    FALSE, TRUE,
  "loto",               "year",    TRUE,  FALSE,
  "logo",               "country", TRUE,  FALSE,
  "lodo",               "drug",    TRUE,  FALSE
)

# amRml's feature types under the names the tabs use; others keep their name.
.AMRML_FEATURE_TYPES <- c(
  gene = "genes", protein = "proteins", Pfam = "domains", COG = "cogs",
  AMRFinder = "args", struct = "struct"
)

# Units the tabs show. MDR targets and leave-one-out modes have no tab yet.
.AMRML_UNITS <- c("drug", "drug_class")


#' A run's manifest, when the folder holds an amRml run
#'
#' @param run_dir Folder to read.
#' @return The parsed `manifest.json`, or NULL when there is none, it isn't an
#'   amRml manifest, or it has no `amRviz_input` artifact.
#' @keywords internal
#' @noRd
.amrml_manifest <- function(run_dir) {
  path <- file.path(run_dir, .AMRML_MANIFEST)
  if (length(run_dir) != 1 || is.na(run_dir) || !file.exists(path)) {
    return(NULL)
  }
  manifest <- tryCatch(jsonlite::read_json(path), error = function(e) NULL)
  if (!identical(manifest[["manifest_type"]], .AMRML_MANIFEST_TYPE) ||
    !is.list(manifest[["artifacts"]][[.AMRML_INPUT]])) {
    return(NULL)
  }
  manifest
}


#' Path to one of a run's tables, from its manifest
#'
#' @param run_dir An amRml run folder.
#' @param table A name in `.AMRML_TABLES`.
#' @param manifest The run's manifest, if already read.
#' @return The table's path, or NA when the manifest doesn't list it.
#' @keywords internal
#' @noRd
.amrml_table_path <- function(run_dir, table,
                              manifest = .amrml_manifest(run_dir)) {
  input <- manifest[["artifacts"]][[.AMRML_INPUT]]
  where <- .AMRML_TABLES[[table]]
  file <- input[[where[1]]][[where[2]]]
  if (!is.character(file) || length(file) != 1) {
    return(NA_character_)
  }
  # Paths are relative to the manifest, under the artifact's directory.
  file.path(run_dir, input[["directory"]] %||% ".", file)
}


#' Is a folder an amRml run?
#'
#' @param dir Folder to check.
#' @return TRUE when `dir` has an amRml manifest listing tasks, matrices and
#'   performance tables that exist.
#' @keywords internal
#' @noRd
.is_amrml_run <- function(dir) {
  manifest <- .amrml_manifest(dir)
  if (is.null(manifest)) {
    return(FALSE)
  }
  paths <- vapply(c("tasks", "matrices", "performance"), .amrml_table_path,
    character(1),
    run_dir = dir, manifest = manifest
  )
  all(!is.na(paths) & file.exists(paths))
}


#' amRml runs inside an ORB
#'
#' @param orb_dir An ORB folder.
#' @return A character vector of the run folders in `{orb}/amRml/`.
#' @keywords internal
#' @noRd
.orb_amrml_runs <- function(orb_dir) {
  runs <- list.dirs(file.path(orb_dir, "amRml"),
    full.names = TRUE, recursive = FALSE
  )
  runs[vapply(runs, .is_amrml_run, logical(1))]
}


#' Read one of a run's tables
#'
#' @param run_dir An amRml run folder.
#' @param table A name in `.AMRML_TABLES`.
#' @param verbose Passed through to .read_parquet_safe().
#' @keywords internal
#' @noRd
.read_amrml_table <- function(run_dir, table, verbose = TRUE) {
  .read_parquet_safe(.amrml_table_path(run_dir, table), verbose)
}


#' Where a run came from: its manifest and its amRdata ORB
#'
#' The ORB is the folder the run sits in, `{orb}/amRml/{run}/`.
#'
#' @param run_dir An amRml run folder.
#' @return A list with `species`, `species_label`, `metadata_path`,
#'   `started_at` and `orb` (a `.describe_orb()` row, or NULL when the ORB
#'   can't be read).
#' @keywords internal
#' @noRd
.amrml_run_record <- function(run_dir) {
  manifest <- .amrml_manifest(run_dir)
  orb_dir <- dirname(dirname(run_dir))
  orb <- .list_available_orbs(
    list.files(orb_dir, pattern = .MANIFEST_PATTERN, full.names = TRUE)
  )
  orb <- if (nrow(orb)) orb[1, ] else NULL
  run <- manifest[["runs"]][[1]]
  started_at <- run[["started_at"]] %||% manifest[["manifest_created_at"]]

  list(
    species = orb$dataset_id %||% manifest[["dataset_id"]] %||% NA_character_,
    species_label = orb$label %||% manifest[["dataset"]][["label"]] %||%
      manifest[["dataset_id"]] %||% NA_character_,
    metadata_path = orb$metadata_path %||% NA_character_,
    # Drop fractions of a second: "2026-11-09 12:48:10.52974" -> "...:10".
    started_at = if (is.character(started_at)) {
      sub("[.][0-9]+$", "", started_at)
    } else {
      NA_character_
    },
    orb = orb
  )
}


#' Every model in a run, under the column names the tabs read
#'
#' One row per task: the tasks table joined to the matrices table, and each
#' mode looked up in `.AMRML_MODES`. Keeps the units and modes the tabs show.
#'
#' @param run_dir An amRml run folder.
#' @param verbose Passed through to .read_parquet_safe().
#' @return A tibble keyed by `task_id`.
#' @keywords internal
#' @noRd
.amrml_models <- function(run_dir, verbose = TRUE) {
  tasks <- .read_amrml_table(run_dir, "tasks", verbose)
  matrices <- .read_amrml_table(run_dir, "matrices", verbose)
  if (!nrow(tasks) || !nrow(matrices)) {
    return(tibble::tibble())
  }
  m <- matrices[match(tasks$matrix_id, matrices$matrix_id), ]
  mode <- .AMRML_MODES[match(m$mode_id, .AMRML_MODES$mode_id), ]
  feature_type <- unname(.AMRML_FEATURE_TYPES[m$feature_type])

  models <- tibble::tibble(
    task_id = tasks$task_id,
    shuffled = tasks$labels == "shuffled",
    seed = tasks$seed,
    drug_label = m$unit,
    drug_or_class = m$target,
    strat_label = ifelse(mode$grouping %in% c("year", "country"),
      mode$grouping, NA_character_
    ),
    strat_value = m$train_group,
    strat_value_test = ifelse(mode$cross_test %in% TRUE,
      m$test_group, NA_character_
    ),
    feature_type = ifelse(is.na(feature_type), m$feature_type, feature_type),
    feature_subtype = m$encoding,
    cross_test = mode$cross_test %in% TRUE,
    LOO = mode$LOO
  )
  unknown <- unique(m$mode_id[is.na(mode$mode_id)])
  if (length(unknown) && isTRUE(verbose)) {
    message("Skipping amRml modes amRviz doesn't know: ", toString(unknown))
  }
  models |>
    dplyr::filter(
      .data$drug_label %in% .AMRML_UNITS,
      .data$LOO %in% FALSE
    ) |>
    dplyr::select(-"LOO")
}


#' Join a run's result rows to its models and source
#'
#' @param results Performance or importance rows.
#' @param models Output of `.amrml_models()`.
#' @param record Output of `.amrml_run_record()`.
#' @keywords internal
#' @noRd
.amrml_rows <- function(results, models, record) {
  if (!nrow(results) || !nrow(models)) {
    return(tibble::tibble())
  }
  # The seed and species come from the tasks and the ORB.
  results <- results[setdiff(names(results), c("seed", "species"))]
  rows <- dplyr::inner_join(models, results, by = "task_id")
  rows$species <- record$species
  rows$species_label <- record$species_label
  rows$metadata_path <- record$metadata_path
  rows
}


#' Load a run's performance rows
#' @inheritParams .amrml_models
#' @keywords internal
#' @noRd
.load_amrml_perf <- function(run_dir, verbose = TRUE) {
  .amrml_rows(
    .read_amrml_table(run_dir, "performance", verbose),
    .amrml_models(run_dir, verbose),
    .amrml_run_record(run_dir)
  )
}


#' Load a run's feature importance rows
#'
#' Reads only the models the tabs show (see `.amrml_models()`), and keeps each
#' one's `.AMRML_TOP_N` most important features.
#'
#' @inheritParams .amrml_models
#' @keywords internal
#' @noRd
.load_amrml_top <- function(run_dir, verbose = TRUE) {
  models <- .amrml_models(run_dir, verbose)
  .amrml_rows(
    .read_top_features(run_dir, models$task_id, verbose),
    models,
    .amrml_run_record(run_dir)
  )
}


#' The most important features of some of a run's models
#'
#' @param run_dir An amRml run folder.
#' @param task_ids The models to read.
#' @param verbose Whether to message when the table can't be read.
#' @param n Features kept per model, by absolute importance.
#' @return Importance rows, at most `n` per model.
#' @keywords internal
#' @noRd
.read_top_features <- function(run_dir, task_ids, verbose = TRUE,
                               n = .AMRML_TOP_N) {
  path <- .amrml_table_path(run_dir, "importance")
  if (is.na(path) || !file.exists(path) || !length(task_ids)) {
    if (isTRUE(verbose) && !is.na(path)) message("File not found: ", path)
    return(tibble::tibble())
  }
  rows <- tryCatch(
    arrow::open_dataset(path, format = "parquet") |>
      dplyr::filter(.data$task_id %in% task_ids) |>
      dplyr::select(dplyr::any_of(c("task_id", "Variable", "Importance", "Sign"))) |>
      dplyr::collect(),
    error = function(e) {
      if (isTRUE(verbose)) {
        message("Failed to read parquet: ", path, " (", conditionMessage(e), ")")
      }
      tibble::tibble()
    }
  )
  if (!nrow(rows)) {
    return(rows)
  }
  # Rank within each model by absolute importance. Sorting on integer codes
  # for the models, not their IDs, keeps this to seconds on tens of millions
  # of rows (sorting the text IDs takes minutes).
  model <- match(rows$task_id, task_ids)
  ord <- order(model, -abs(rows$Importance), method = "radix")
  rank <- sequence(rle(model[ord])$lengths)
  rows[ord[rank <= n], , drop = FALSE]
}
