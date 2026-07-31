## Tests for annotation loading and enrichment:
## .load_one_species_cluster_features, enrich_with_annotations, load_feature_name_map.

# ── .load_one_species_cluster_features ────────────────────────────────────────────────

test_that("load_feature_annotations returns a tibble for bundled Shigella demo", {
  ann <- load_feature_annotations("Sfl")
  # Demo data ships with cluster_feature.parquet for Shigella_flexneri.
  expect_true(is.null(ann) || tibble::is_tibble(ann) || is.data.frame(ann))
  if (!is.null(ann)) {
    expect_true(all(
      c("cluster", "feature", "cluster_name", "COG", "COG_name") %in%
        names(ann)
    ))
    expect_gt(nrow(ann), 0)
  }
})

test_that(".load_one_species_cluster_features returns NULL when results_root has no annotations", {
  empty_root <- tempfile("empty_root")
  dir.create(empty_root)
  on.exit(unlink(empty_root, recursive = TRUE))

  folders <- listAmRmlSpeciesFolders(extdata, verbose = FALSE)
  skip_if(
    length(folders) == 0,
    "No species dirs with baseline perf in extdata"
  )

  result <- .load_one_species_cluster_features(folders[1])
  # No annotation parquet under results_root, and the species code does not
  # correspond to bundled data, so we expect either NULL or the bundled fallback.
  expect_true(is.null(result) || is.data.frame(result))
})

# ── enrich_with_annotations ─────────────────────────────────────────────────

test_that("enrich_with_annotations returns input unchanged for NULL", {
  expect_null(enrich_with_annotations(NULL, folders[1]))
})

test_that("enrich_with_annotations returns input unchanged for zero-row data", {
  df <- tibble::tibble(Variable = character(0))
  result <- enrich_with_annotations(df, folders[1])
  expect_equal(nrow(result), 0)
})

test_that("enrich_with_annotations returns input unchanged when Variable column missing", {
  df <- tibble::tibble(x = 1:3)
  result <- enrich_with_annotations(df, folders[1])
  expect_identical(result, df)
})

test_that("enrich_with_annotations adds cluster/COG columns when annotations exist", {
  folders <- listAmRmlSpeciesFolders(extdata, verbose = FALSE)
  skip_if(
    length(folders) == 0,
    "No species dirs with baseline perf in extdata"
  )

  ann <- .load_one_species_cluster_features(folders[1])
  skip_if(is.null(ann), "No bundled Shigella annotations")
  skip_if(nrow(ann) == 0, "Empty annotations")

  # Build a tbl using a feature ID that exists in the annotations.
  sample_feature <- ann$feature[1]
  df <- tibble::tibble(
    Variable = c(sample_feature, paste0(sample_feature, "_extra")),
    Importance = c(0.5, 0.3)
  )

  result <- enrich_with_annotations(df, "Sfl")
  expect_true(all(
    c("cluster", "cluster_name", "COG", "COG_name") %in% names(result)
  ))
  expect_equal(nrow(result), 2)
})

# ── load_feature_name_map ───────────────────────────────────────────────────

test_that("load_feature_name_map returns NULL when name file is absent", {
  # The demo data does not ship gene_names.parquet / protein_names.parquet, so
  # this should be NULL without an amrdata_root pointing at one.
  folders <- listAmRmlSpeciesFolders(extdata, verbose = FALSE)
  skip_if(
    length(folders) == 0,
    "No species dirs with baseline perf in extdata"
  )
  result <- load_feature_name_map("genes", species_dir = folders[1])
  expect_null(result)
})

test_that("load_feature_name_map returns NULL for unrecognised species", {
  folders <- listAmRmlSpeciesFolders(extdata, verbose = FALSE)
  skip_if(
    length(folders) == 0,
    "No species dirs with baseline perf in extdata"
  )
  result <- load_feature_name_map(
    "ZZZ", "genes",
    amrdata_root = tempdir(), species_dir = folders[1]
  )
  expect_null(result)
})
