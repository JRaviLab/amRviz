# Regenerate the bundled feature-name maps ({gene,protein,domain}_names.parquet)
# by keeping only the ids that a top-features Variable can actually reference.
# Run from the package source root, AFTER subset_demo_top_features.R.
#
# amRml ships whole-proteome/-genome name maps (e.g. Shigella sonnei's
# `protein_names.parquet` is ~9.6M rows / ~10 MB). The dashboard only ever uses
# them to relabel features that appear in a top-features plot: Feature
# Importance looks up `name_map$Variable` for the features it is about to draw
# and ignores the rest. So any id that never appears as a top-features
# `Variable` is dead weight and is dropped.
#
# Id formats differ by scale and must be reconciled with the top-features
# `Variable` column:
#   * protein -- name map `proteinID` uses "fig|624...", top-features `Variable`
#     uses "fig.624..."; we normalise "|" -> "." before matching.
#   * domain  -- top-features `Variable` is "PF21279_IPR..."; the name map key
#     is the part before the first "_" ("PF21279").
#   * gene    -- ids match directly.
# Matching on the normalised form (rather than the raw proteinID) keeps the map
# correct even though the current runtime lookup does not itself normalise "|".

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(stringr)
})

# name-map file stem -> (id column, top-features feature_type it labels)
specs <- list(
  gene    = list(id = "Gene",      feature_type = "genes"),
  protein = list(id = "proteinID", feature_type = "proteins"),
  domain  = list(id = "DB.ID",     feature_type = "domains")
)

# Normalise a top-features Variable vector into the name-map's id space.
to_key <- function(variable, scale) {
  if (scale == "domain") {
    stringr::str_split_i(variable, "_", 1)
  } else {
    variable
  }
}
# Normalise a name-map id vector into the top-features Variable space.
from_id <- function(id, scale) {
  if (scale == "protein") gsub("[|]", ".", id) else id
}

species_dirs <- list.dirs("inst/extdata", recursive = FALSE)

for (sp_dir in species_dirs) {
  # Variables actually present in this species' trimmed top-features, by scale.
  tf_files <- list.files(
    sp_dir,
    pattern = "^(all|country|year)_top_features[.]parquet$", full.names = TRUE
  )
  if (!length(tf_files)) next
  tf <- dplyr::bind_rows(lapply(tf_files, function(f) {
    arrow::read_parquet(f, col_select = c("Variable", "feature_type"))
  }))

  for (scale in names(specs)) {
    id_col <- specs[[scale]]$id
    ftype <- specs[[scale]]$feature_type
    fp <- file.path(sp_dir, paste0(scale, "_names.parquet"))
    if (!file.exists(fp)) next

    keys <- unique(to_key(tf$Variable[tf$feature_type == ftype], scale))
    orig_size <- file.size(fp)
    full <- arrow::read_parquet(fp)
    if (!id_col %in% names(full)) next

    kept <- full[from_id(full[[id_col]], scale) %in% keys, , drop = FALSE]
    arrow::write_parquet(kept, fp)

    message(sprintf(
      "%s: %d -> %d rows (%.2f -> %.2f MB)",
      sub("inst/extdata/", "", fp), nrow(full), nrow(kept),
      orig_size / 1e6, file.size(fp) / 1e6
    ))
  }
}
