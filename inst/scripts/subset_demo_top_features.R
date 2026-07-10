# Regenerate the bundled top-features parquets by trimming them to only the rows
# the dashboard can actually surface. Run from the package source root.
#
# amRml emits two flavours of `*_top_features.parquet` per species, both far
# larger than a Bioconductor demo needs:
#
#   * baseline  `all_top_features.parquet`   -- has a `shuffled` column and a
#     `seed` column; `strat_label` is NA (no country/year stratification).
#   * stratified `country`/`year_top_features.parquet` -- one row per
#     resampling replicate (encoded in `filename`/`seed_from_name`), with
#     `strat_label`/`strat_value` naming the trained-on stratum.
#
# Two facts let us shrink both losslessly from the dashboard's point of view:
#
#   1. `shuffled` rows are the permutation-null baseline. No amRviz code path
#      reads the `shuffled` column, so shuffled == TRUE rows can never appear in
#      a plot and are dropped. (The `filename` and `seed_from_name` columns on
#      the stratified files are likewise read nowhere, so they are dropped too.)
#
#   2. Every view collapses replicates with max/mean(Importance) per feature and
#      then keeps at most the top N. The Feature Importance and Cross-model
#      sliders both cap N at 100, so no view can ask for more than the top 100
#      features of a model. A "model" at the finest user-selectable granularity
#      is species x drug_or_class x feature_type (scale) x feature_subtype
#      (encoding), plus strat_label x strat_value for the stratified files.
#
# We rank features by their max-across-replicate importance within each model
# and keep every replicate row of the top 100. Keeping all replicate rows (not
# an aggregated summary) means the app's own max()/mean(Importance) still return
# the exact values they would on the full file -- the trim is invisible to the
# dashboard. Ranking per (scale, subtype[, stratum]) independently makes the
# kept set a superset of what any single view slices (Feature Importance folds
# `struct` in alongside the chosen scale; the Network and Cross-model views fold
# scales/subtypes together), so no view loses a feature it would have shown.
#
# Abhirupa (upstream data producer) has confirmed the `shuffled` rows are not
# intended as a dashboard input.

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
})

top_n <- 100
# Columns present on the stratified files that no amRviz code path reads.
drop_cols <- c("filename", "seed_from_name")

# All baseline + stratified top-features files; MDR/LOO/cross are loaded
# differently (or not at all) and are left untouched.
tf_files <- list.files(
  "inst/extdata",
  pattern = "^(all|country|year)_top_features[.]parquet$",
  recursive = TRUE, full.names = TRUE
)

for (tf_file in tf_files) {
  orig_size <- file.size(tf_file)
  full <- arrow::read_parquet(tf_file)

  # Drop the permutation-null rows and the unread replicate-bookkeeping columns.
  kept <- full
  if ("shuffled" %in% names(kept)) {
    kept <- dplyr::filter(kept, is.na(.data$shuffled) | !.data$shuffled)
  }
  kept <- dplyr::select(kept, -dplyr::any_of(drop_cols))

  # Model granularity: add the stratification columns when present.
  model_cols <- c("species", "drug_or_class", "feature_type", "feature_subtype")
  strat_cols <- intersect(c("strat_label", "strat_value"), names(kept))
  model_cols <- c(model_cols, strat_cols)

  # Rank features by max-across-replicate importance within each model, then
  # keep every row belonging to the top `top_n` features of that model.
  top_features <- kept |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(model_cols, "Variable")))) |>
    dplyr::summarise(
      .maxImp = max(.data$Importance, na.rm = TRUE), .groups = "drop"
    ) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(model_cols))) |>
    dplyr::slice_max(order_by = .data$.maxImp, n = top_n, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::select(dplyr::all_of(c(model_cols, "Variable")))

  kept <- kept |>
    dplyr::semi_join(top_features, by = c(model_cols, "Variable"))

  arrow::write_parquet(kept, tf_file)

  message(sprintf(
    "%s: %d -> %d rows (%.2f -> %.2f MB)",
    sub("inst/extdata/", "", tf_file), nrow(full), nrow(kept),
    orig_size / 1e6, file.size(tf_file) / 1e6
  ))
}
