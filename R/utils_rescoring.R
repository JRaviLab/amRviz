# Rescored features and dyads: amRml's feature rescoring of a run.
#
# After a run finishes, amRml can rescore its features across seeds (models
# that don't beat shuffled labels are dropped first) and map them onto the
# ORB's dyads. The run's manifest lists where the tables are, as its
# `feature_rescoring` artifact. amRviz reads two of them:
#
# * `features_all_seeds`: one row per drug or class, feature and encoding that
#   is in a model's core group, with its seed ratio, median rank score and
#   share of the model's importance.
# * `dyads_top`: one row per drug or class and dyad reached by those features,
#   with the dyad's median rank score and the features that reach it.
# * `models_qualified`: the model fits those features were scored from, the
#   ones that passed amRml's performance gate, with their MCC and the MCC of
#   the shuffled-label fit.

# The manifest entry that lists the rescoring tables, and the files read.
.AMRML_RESCORING <- "feature_rescoring"
.RESCORING_FEATURES <- "features_all_seeds.parquet"
.RESCORING_DYADS <- "dyads_top.parquet"
.RESCORING_MODELS <- "models_qualified.parquet"

# The feature types the rescoring tables use, under the names the tabs use,
# and the prefixes a dyad's `variables_csv` gives them.
.RESCORING_TYPES <- c(
  AMRFinder = "args", COG = "cogs", Pfam = "domains",
  gene = "genes", protein = "proteins"
)
.RESCORING_PREFIXES <- c(
  amr = "args", cog = "cogs", pfam = "domains",
  gene = "genes", protein = "proteins"
)


#' Folder of a run's rescoring tables, from its manifest
#'
#' @param run_dir An amRml run folder.
#' @return The folder, or NA when the manifest lists no finished rescoring or
#'   its tables aren't there.
#' @keywords internal
#' @noRd
.amrml_rescoring_dir <- function(run_dir) {
  entry <- .amrml_manifest(run_dir)[["artifacts"]][[.AMRML_RESCORING]]
  dir <- entry[["directory"]]
  if (!is.list(entry) || !identical(entry[["status"]], "ready") ||
    !is.character(dir) || length(dir) != 1) {
    return(NA_character_)
  }
  # Paths are relative to the manifest.
  dir <- file.path(run_dir, dir)
  files <- file.path(dir, c(.RESCORING_FEATURES, .RESCORING_DYADS))
  if (all(file.exists(files))) dir else NA_character_
}


#' Load one run's rescoring
#'
#' @param run_dir An amRml run folder.
#' @param verbose Whether to message when a table can't be read.
#' @return A list of `features`, `dyads` and `models` tibbles, with each row's
#'   `species` and `metadata_path` from the run's ORB, or NULL when the run has
#'   no rescoring. Only drugs and drug classes are kept; MDR targets have no
#'   tab yet. `models` is empty when the run has no `models_qualified` table.
#' @keywords internal
#' @noRd
.load_rescoring_source <- function(run_dir, verbose = TRUE) {
  dir <- .amrml_rescoring_dir(run_dir)
  if (is.na(dir)) {
    return(NULL)
  }
  record <- .amrml_run_record(run_dir)
  say <- verbose
  read <- function(file, verbose = say) {
    rows <- .read_parquet_safe(file.path(dir, file), verbose)
    if (!nrow(rows)) {
      return(rows)
    }
    rows <- rows[rows$drug_label %in% .AMRML_UNITS, , drop = FALSE]
    rows$species <- record$species
    rows$metadata_path <- record$metadata_path
    rows
  }
  named <- function(rows) {
    if (nrow(rows)) {
      rows$feature_type <- unname(.RESCORING_TYPES[rows$feature_type])
    }
    rows
  }
  list(
    features = named(read(.RESCORING_FEATURES)),
    dyads = read(.RESCORING_DYADS),
    # Optional: the thresholds on model performance need it.
    models = named(read(.RESCORING_MODELS, verbose = FALSE))
  )
}


#' Load the rescoring of the selected runs
#'
#' Follows the same modes as loadTopFeat(): the selected result folders when a
#' results root is set, nothing until one is selected, and the packaged demo
#' run without a root.
#'
#' @param results_root User results root, or NULL for the packaged demo data.
#' @param result_dirs Selected result folders, or NULL.
#' @param verbose Passed through to the readers.
#' @return A list of `features`, `dyads` and `models` tibbles bound across the
#'   runs, or NULL when none of them has rescoring.
#' @keywords internal
#' @noRd
loadRescoring <- function(results_root = NULL, result_dirs = NULL,
                          verbose = TRUE) {
  rr <- .normalize_results_root(results_root)
  dirs <- if (!is.null(rr)) {
    result_dirs
  } else {
    extdata <- system.file("extdata", package = "amRviz")
    if (nzchar(extdata)) unname(listResultSources(extdata))
  }
  found <- Filter(Negate(is.null), lapply(dirs, .load_rescoring_source,
    verbose = verbose
  ))
  if (!length(found)) {
    return(NULL)
  }
  list(
    features = dplyr::bind_rows(lapply(found, `[[`, "features")),
    dyads = dplyr::bind_rows(lapply(found, `[[`, "dyads")),
    models = dplyr::bind_rows(lapply(found, `[[`, "models"))
  )
}


#' Default thresholds for filtering the rescoring
#'
#' The package's own selection. Model fits must reach an MCC of 0.5 and beat the
#' shuffled labels, which every fit in the rescoring already does, so these two
#' can only be raised. A feature must be in the core group, within the
#' cumulative share limit, and reach the 75th percentile of its drug's rank
#' scores (which is what amRml's `features_high_contribution` holds).
#'
#' @return A list of `mcc`, `delta_mcc` (the least MCC above the shuffled
#'   labels), `quantile`, `core`, `cumulative`, `both_subtypes` and
#'   `dyad_rank_score`.
#' @keywords internal
#' @noRd
rescoringDefaults <- function() {
  list(
    mcc = 0.5, delta_mcc = 0, quantile = 0.75, core = 0.9, cumulative = 0.75,
    both_subtypes = FALSE, dyad_rank_score = 0
  )
}


#' Thresholds that keep everything in the rescoring
#'
#' The Network tab shows one slider per category and leaves the other
#' thresholds at these values, so only the sliders filter: any fit that passed
#' amRml's own gate (MCC 0.5 and beating the shuffled labels), every feature
#' row, and every dyad.
#'
#' @return A list as from rescoringDefaults().
#' @keywords internal
#' @noRd
rescoringAll <- function() {
  list(
    mcc = 0.5, delta_mcc = 0, quantile = 0, core = 1, cumulative = 1,
    both_subtypes = FALSE, dyad_rank_score = 0
  )
}


#' Model fits of the rescoring that pass the performance thresholds
#'
#' @param models The `models` of loadRescoring().
#' @param bug Species code.
#' @param thresholds A list as from rescoringDefaults().
#' @return The fits whose MCC reaches `mcc` and is at least `delta_mcc` above
#'   the shuffled-label fit's; zero rows when there are none.
#' @keywords internal
#' @noRd
rescoredModels <- function(models, bug, thresholds = rescoringDefaults()) {
  if (is.null(models) || !nrow(models)) {
    return(tibble::tibble())
  }
  models |>
    dplyr::filter(
      normalize_species(.data$species) %in% normalize_species(bug),
      .data$nonshuffled_MCC >= thresholds$mcc - 1e-12,
      .data$MCC_diff >= thresholds$delta_mcc - 1e-12
    )
}


#' MCC of each drug's model fits that pass the thresholds
#'
#' @param models Output of rescoredModels().
#' @return One row per drug or class: `drug_label`, `drug_or_class`, `mcc` and
#'   `delta_mcc` (medians over its fits), and `n_fits`.
#' @keywords internal
#' @noRd
rescoredDrugScores <- function(models) {
  if (!nrow(models)) {
    return(tibble::tibble(
      drug_label = character(), drug_or_class = character(), mcc = numeric(),
      delta_mcc = numeric(), n_fits = integer()
    ))
  }
  models |>
    dplyr::group_by(.data$drug_label, .data$drug_or_class) |>
    dplyr::summarise(
      mcc = stats::median(.data$nonshuffled_MCC),
      delta_mcc = stats::median(.data$MCC_diff),
      n_fits = dplyr::n(),
      .groups = "drop"
    )
}


#' Features of the rescoring that pass the thresholds
#'
#' Follows the selection in amRml's rescoring. Per drug or class, a feature row
#' is kept when its rank score reaches `quantile` of the drug's rows, it is in
#' the core group (the share before it is below `core`) and its cumulative
#' share is at most `cumulative`. Rows are then collapsed across encodings,
#' keeping the best, and optionally kept only when they pass in both.
#'
#' @param features The `features` of loadRescoring().
#' @param bug Species code.
#' @param thresholds A list as from rescoringDefaults().
#' @param models Output of rescoredModels(). When given, a feature row is only
#'   kept if its drug and feature type and encoding still has a model fit that
#'   passes the performance thresholds.
#' @return One row per drug or class and feature: `drug_label`,
#'   `drug_or_class`, `feature_type`, `Variable`, `rank_score`, `contribution`,
#'   `seed_ratio` and `sign`, with `metadata_path`; zero rows when nothing
#'   passes.
#' @keywords internal
#' @noRd
rescoredFeatures <- function(features, bug, thresholds = rescoringDefaults(),
                             models = NULL) {
  empty <- tibble::tibble(
    drug_label = character(), drug_or_class = character(),
    feature_type = character(), Variable = character(),
    rank_score = numeric(), contribution = numeric(), seed_ratio = numeric(),
    sign = character(), metadata_path = character()
  )
  if (is.null(features) || !is.data.frame(features) || !nrow(features)) {
    return(empty)
  }
  rows <- features |>
    dplyr::filter(normalize_species(.data$species) %in% normalize_species(bug))
  if (!is.null(models)) {
    rows <- dplyr::semi_join(rows, models, by = c(
      "drug_label", "drug_or_class", "feature_type", "feature_subtype"
    ))
  }
  rows <- rows |>
    dplyr::group_by(.data$drug_label, .data$drug_or_class) |>
    # The cutoff comes from every row of the drug, before the other filters.
    dplyr::mutate(
      cutoff = stats::quantile(.data$median_rank_score, thresholds$quantile)
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(
      .data$median_rank_score >= .data$cutoff - 1e-12,
      .data$median_prior_cum_contrib < thresholds$core,
      .data$median_cum_contrib <= thresholds$cumulative
    )
  if (isTRUE(thresholds$both_subtypes)) {
    rows <- rows |>
      dplyr::group_by(
        .data$drug_label, .data$drug_or_class, .data$feature_type,
        .data$variable
      ) |>
      dplyr::filter(dplyr::n_distinct(.data$feature_subtype) > 1) |>
      dplyr::ungroup()
  }
  if (!nrow(rows)) {
    return(empty)
  }
  rows |>
    dplyr::group_by(
      .data$drug_label, .data$drug_or_class, .data$feature_type, .data$variable
    ) |>
    dplyr::slice_max(
      order_by = .data$median_rank_score * 10 + .data$median_contribution,
      n = 1, with_ties = FALSE
    ) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      .data$drug_label, .data$drug_or_class, .data$feature_type,
      Variable = .data$variable, rank_score = .data$median_rank_score,
      contribution = .data$median_contribution,
      seed_ratio = .data$seed_ratio, sign = .data$sign, .data$metadata_path
    )
}


#' Dyads of the rescoring that reach some features
#'
#' @param dyads The `dyads` of loadRescoring().
#' @param features Output of rescoredFeatures(), cut to the features shown and
#'   with `Variable` named as the network names it.
#' @param thresholds A list as from rescoringDefaults().
#' @return One row per drug or class, dyad and feature the dyad is reached by:
#'   `drug_label`, `drug_or_class`, `dyad` (as "protein|gene"), `feature_type`,
#'   `Variable`, `dyad_rank_score`, `dyad_contribution` and `frequency` (how
#'   many features reach the dyad); zero rows when there are none.
#' @keywords internal
#' @noRd
rescoredDyads <- function(dyads, features, thresholds = rescoringDefaults()) {
  empty <- tibble::tibble(
    drug_label = character(), drug_or_class = character(), dyad = character(),
    feature_type = character(), Variable = character(),
    dyad_rank_score = numeric(), dyad_contribution = numeric(),
    frequency = numeric()
  )
  if (is.null(dyads) || !nrow(dyads) || !nrow(features)) {
    return(empty)
  }
  # The features that reach each dyad are "prefix:name", comma-separated.
  reach <- dyads |>
    dplyr::semi_join(features, by = c("drug_label", "drug_or_class")) |>
    dplyr::filter(
      .data$dyad_median_rank_score >= thresholds$dyad_rank_score - 1e-12
    ) |>
    dplyr::transmute(
      .data$drug_label, .data$drug_or_class, .data$dyad,
      dyad_rank_score = .data$dyad_median_rank_score,
      dyad_contribution = .data$dyad_median_contribution,
      frequency = .data$frequency,
      variable = strsplit(.data$variables_csv, ",", fixed = TRUE)
    ) |>
    tidyr::unnest("variable")
  prefix <- sub(":.*$", "", reach$variable)
  reach$feature_type <- unname(.RESCORING_PREFIXES[prefix])
  reach$Variable <- sub("^[a-z]+:", "", reach$variable)
  # Only the features shown, as the network names them.
  shown <- dplyr::distinct(
    features, .data$drug_label, .data$drug_or_class, .data$feature_type,
    .data$Variable
  )
  reach |>
    dplyr::filter(!is.na(.data$feature_type)) |>
    dplyr::mutate(
      Variable = .network_feature_name(.data$feature_type, .data$Variable)
    ) |>
    dplyr::inner_join(
      shown,
      by = c("drug_label", "drug_or_class", "feature_type", "Variable")
    ) |>
    dplyr::select(dplyr::all_of(names(empty)))
}
