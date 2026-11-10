# Reading amRml runs: the results of one amRml run on one amRdata ORB.
#
# Everything about a model comes from the run's own tables, never from file or
# folder names: each result row's `task_id` joins to tasks.parquet (labels and
# seed), whose `matrix_id` joins to matrices.parquet (mode, drug or class,
# feature type, encoding, groups). The run's species and metadata come from its
# amRdata ORB's manifest, found through run.json.
#
# The tabs still read old amRml's column names, so rows are returned under
# those, with the same mapping as amRml's temporary `exportMergedResults()`.


# File names amRml gives a run's tables, kept in one place.
.AMRML_FILES <- c(
  tasks = "tasks.parquet",
  matrices = "matrices.parquet",
  performance = "results/performance.parquet",
  importance = "results/importance.parquet",
  run = "run.json"
)

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


#' Is a folder an amRml run?
#'
#' @param dir Folder to check.
#' @return TRUE when `dir` holds a run's tasks, matrices and performance tables.
#' @keywords internal
#' @noRd
.is_amrml_run <- function(dir) {
  length(dir) == 1 && !is.na(dir) && dir.exists(dir) &&
    all(file.exists(file.path(
      dir, .AMRML_FILES[c("tasks", "matrices", "performance")]
    )))
}


#' amRml runs inside an ORB
#'
#' amRml writes a run to `{orb}/amRml/{date-time}/` unless told otherwise.
#'
#' @param orb_dir An ORB folder.
#' @return A character vector of run folders.
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
#' @param table A name in `.AMRML_FILES`.
#' @param verbose Passed through to .read_parquet_safe().
#' @keywords internal
#' @noRd
.read_amrml_table <- function(run_dir, table, verbose = TRUE) {
  .read_parquet_safe(file.path(run_dir, .AMRML_FILES[[table]]), verbose)
}


#' Where a run came from: its record and its amRdata ORB
#'
#' Reads `run.json` when the run has one (amRml from PR #53 on). Its `orb`
#' entry names the ORB's folder and manifest; a relative folder, or a run
#' without `run.json`, falls back to the folder amRml writes runs under by
#' default, `{orb}/amRml/{run}/`.
#'
#' @param run_dir An amRml run folder.
#' @return A list with `species`, `species_label`, `metadata_path`,
#'   `started_at` and `orb` (a `.describe_orb()` row, or NULL when the ORB
#'   can't be found).
#' @keywords internal
#' @noRd
.amrml_run_record <- function(run_dir) {
  run_json <- file.path(run_dir, .AMRML_FILES[["run"]])
  run <- if (file.exists(run_json)) {
    tryCatch(jsonlite::read_json(run_json), error = function(e) NULL)
  }
  recorded <- run[["orb"]]

  orb_dir <- recorded[["directory"]]
  absolute <- is.character(orb_dir) && length(orb_dir) == 1 &&
    grepl("^(/|~|[A-Za-z]:)", orb_dir)
  if (!absolute || !dir.exists(orb_dir)) {
    orb_dir <- dirname(dirname(run_dir))
  }
  manifests <- if (is.character(recorded[["manifest"]])) {
    file.path(orb_dir, recorded[["manifest"]])
  } else {
    list.files(orb_dir, pattern = .MANIFEST_PATTERN, full.names = TRUE)
  }
  orb <- .list_available_orbs(manifests)
  orb <- if (nrow(orb)) orb[1, ] else NULL

  list(
    species = run[["species"]] %||% orb$dataset_id %||% NA_character_,
    species_label = orb$label %||% run[["orb"]][["dataset_label"]] %||%
      run[["species"]] %||% NA_character_,
    metadata_path = orb$metadata_path %||% NA_character_,
    started_at = run[["started_at"]] %||% NA_character_,
    orb = orb
  )
}


#' Every model in a run, under the column names the tabs read
#'
#' One row per task: tasks.parquet joined to matrices.parquet, and each mode
#' looked up in `.AMRML_MODES`. Keeps the units and modes the tabs show.
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
#' @param results performance.parquet or importance.parquet rows.
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
#' @inheritParams .amrml_models
#' @keywords internal
#' @noRd
.load_amrml_top <- function(run_dir, verbose = TRUE) {
  .amrml_rows(
    .read_amrml_table(run_dir, "importance", verbose),
    .amrml_models(run_dir, verbose),
    .amrml_run_record(run_dir)
  )
}
