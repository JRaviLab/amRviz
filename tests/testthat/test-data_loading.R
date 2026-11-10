## Tests for data loading functions: loadMLResults, loadTopFeat,
## listResultSources, loadDrugClassMap

# ── listResultSources ──────────────────────────────────────────────────

test_that("listResultSources returns empty for NULL/invalid root", {
  expect_equal(listResultSources(NULL), character(0))
  expect_equal(listResultSources(""), character(0))
  expect_equal(
    listResultSources("/nonexistent/path"),
    character(0)
  )
})

test_that("listResultSources finds the older-format demo folders", {
  extdata <- system.file("extdata", package = "amRviz")
  skip_if(!nzchar(extdata), "No extdata directory found")

  result <- listResultSources(extdata, verbose = FALSE)
  expect_type(result, "character")
  expect_true(length(result) > 0)
  expect_true(all(nzchar(names(result))))
})

test_that("listResultSources ignores dirs without results", {
  tmp <- tempfile("test_species_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))

  empty_dir <- file.path(tmp, "Empty_species")
  dir.create(empty_dir)

  result <- listResultSources(tmp, verbose = FALSE)
  expect_equal(length(result), 0)
})

# ── loadMLResults ────────────────────────────────────────────────────────────

test_that("loadMLResults returns empty tibble when results_root given but no result_dirs", {
  tmp <- tempdir()
  result <- suppressMessages(
    loadMLResults(
      results_root = tmp,
      result_dirs = NULL,
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

test_that("loadMLResults loads from specific result_dirs", {
  extdata <- system.file("extdata", package = "amRviz")
  skip_if(!nzchar(extdata), "No extdata directory found")

  folders <- listResultSources(extdata, verbose = FALSE)
  skip_if(
    length(folders) == 0,
    "No species dirs with baseline perf in extdata"
  )

  result <- suppressMessages(
    loadMLResults(
      results_root = extdata,
      result_dirs = folders[1],
      verbose = FALSE
    )
  )
  expect_s3_class(result, "tbl_df")
  expect_gt(nrow(result), 0)
  expect_true("species_label" %in% names(result))
})

# ── loadTopFeat ──────────────────────────────────────────────────────────────

test_that("loadTopFeat returns empty tibble when results_root given but no result_dirs", {
  tmp <- tempdir()
  result <- suppressMessages(
    loadTopFeat(
      results_root = tmp,
      result_dirs = NULL,
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
