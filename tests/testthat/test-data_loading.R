## Tests for data loading functions: loadMLResults, loadTopFeat,
## listAmRmlSpeciesFolders, get_metadata_path, loadDrugClassMap

# ── listAmRmlSpeciesFolders ──────────────────────────────────────────────────

test_that("listAmRmlSpeciesFolders returns empty for NULL/invalid root", {
  expect_equal(listAmRmlSpeciesFolders(NULL), character(0))
  expect_equal(listAmRmlSpeciesFolders(""), character(0))
  expect_equal(
    listAmRmlSpeciesFolders("/nonexistent/path"),
    character(0)
  )
})

test_that("listAmRmlSpeciesFolders discovers species dirs in extdata", {
  extdata <- system.file("extdata", package = "amRviz")
  skip_if(!nzchar(extdata), "No extdata directory found")

  result <- listAmRmlSpeciesFolders(extdata, verbose = FALSE)
  expect_type(result, "character")
  expect_true(length(result) > 0)
  expect_true(all(nzchar(names(result))))
})

test_that("listAmRmlSpeciesFolders ignores dirs without perf parquets", {
  tmp <- tempfile("test_species_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))

  empty_dir <- file.path(tmp, "Empty_species")
  dir.create(empty_dir)

  result <- listAmRmlSpeciesFolders(tmp, verbose = FALSE)
  expect_equal(length(result), 0)
})

# ── loadMLResults ────────────────────────────────────────────────────────────

test_that("loadMLResults returns empty tibble when results_root given but no species_dirs", {
  tmp <- tempdir()
  result <- suppressMessages(
    loadMLResults(
      results_root = tmp,
      species_dirs = NULL,
      verbose = FALSE
    )
  )
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 0)
})

test_that("loadMLResults loads demo data in fallback mode", {
  result <- suppressMessages(loadMLResults(verbose = FALSE))
  expect_s3_class(result, "tbl_df")
  expect_gt(nrow(result), 0)
  expect_true("species" %in% names(result))
})

test_that("loadMLResults loads from specific species_dirs", {
  extdata <- system.file("extdata", package = "amRviz")
  skip_if(!nzchar(extdata), "No extdata directory found")

  folders <- listAmRmlSpeciesFolders(extdata, verbose = FALSE)
  skip_if(
    length(folders) == 0,
    "No species dirs with baseline perf in extdata"
  )

  result <- suppressMessages(
    loadMLResults(
      results_root = extdata,
      species_dirs = folders[1],
      verbose = FALSE
    )
  )
  expect_s3_class(result, "tbl_df")
  expect_gt(nrow(result), 0)
  expect_true("species_label" %in% names(result))
})

# ── loadTopFeat ──────────────────────────────────────────────────────────────

test_that("loadTopFeat returns empty tibble when results_root given but no species_dirs", {
  tmp <- tempdir()
  result <- suppressMessages(
    loadTopFeat(
      results_root = tmp,
      species_dirs = NULL,
      verbose = FALSE
    )
  )
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 0)
})

test_that("loadTopFeat loads demo data in fallback mode", {
  result <- suppressMessages(loadTopFeat(verbose = FALSE))
  expect_s3_class(result, "tbl_df")
  expect_gt(nrow(result), 0)
  expect_true("Variable" %in% names(result))
  expect_true("Importance" %in% names(result))
})

# ── get_metadata_path ────────────────────────────────────────────────────────

test_that("get_metadata_path returns NULL for nonexistent species code", {
  result <- get_metadata_path("Zzz")
  expect_null(result)
})

test_that("get_metadata_path finds demo metadata parquet", {
  extdata <- system.file("extdata", package = "amRviz")
  skip_if(!nzchar(extdata), "No extdata directory found")

  result <- get_metadata_path("Sfl")
  if (!is.null(result)) {
    expect_true(file.exists(result))
    expect_true(grepl("Sfl_metadata\\.parquet$", result))
  }
})

# ── loadDrugClassMap ─────────────────────────────────────────────────────────

test_that("loadDrugClassMap returns a tibble with expected columns", {
  result <- suppressMessages(loadDrugClassMap())
  expect_s3_class(result, "tbl_df")
  expect_true("drug.antibiotic_name" %in% names(result))
  expect_true("drug_class" %in% names(result))
  expect_gt(nrow(result), 0)
})

test_that("loadDrugClassMap returns distinct rows", {
  result <- suppressMessages(loadDrugClassMap())
  expect_equal(nrow(result), nrow(dplyr::distinct(result)))
})
