# Regenerate the bundled COG annotation by trimming it to only the features
# referenced by the demo top-features parquets. Run from the package source
# root

dir <- "inst/extdata/Shigella_flexneri"
ann_file <- file.path(dir, "cluster_feature.parquet")

# Each top-features file's Variable column has the form "<prefix>_<suffix>"
# (e.g. "PF21279_IPR056912"). enrich_with_annotations() joins on the prefix,
# so any annotation row whose `feature` doesn't appear as a prefix here can
# never match and is safe to drop.
tf_files <- list.files(
  dir,
  pattern = "_ML_top_features\\.parquet$",
  full.names = TRUE
)

join_keys <- tf_files |>
  lapply(function(fp) arrow::read_parquet(fp, col_select = "Variable")$Variable) |>
  unlist() |>
  stringr::str_split_i("_", 1) |>
  unique()

arrow::read_parquet(ann_file) |>
  dplyr::filter(feature %in% join_keys) |>
  arrow::write_parquet(ann_file)
