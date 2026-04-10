## Tests for internal data helper functions

# ── .normalize_results_root ──────────────────────────────────────────────────

test_that(".normalize_results_root returns NULL for invalid inputs", {
  expect_null(.normalize_results_root(NULL))
  expect_null(.normalize_results_root(NA))
  expect_null(.normalize_results_root(""))
  expect_null(.normalize_results_root(character(0)))
  # Multiple values
  expect_null(.normalize_results_root(c("/a", "/b")))
})

test_that(".normalize_results_root normalizes a valid path", {
  tmp <- tempdir()
  result <- .normalize_results_root(tmp)
  expect_type(result, "character")
  expect_equal(nchar(result) > 0, TRUE)
  # Should be an absolute path
  expect_true(startsWith(result, "/") || grepl("^[A-Z]:", result))
})

# ── .read_parquet_safe ───────────────────────────────────────────────────────

test_that(".read_parquet_safe returns empty tibble for missing file", {
  result <- suppressMessages(.read_parquet_safe("/nonexistent/path.parquet", verbose = FALSE))
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 0)
})

test_that(".read_parquet_safe returns empty tibble for NULL path", {
  result <- suppressMessages(.read_parquet_safe(NULL, verbose = FALSE))
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 0)
})

test_that(".read_parquet_safe reads a valid parquet file", {
  extdata <- system.file("extdata", package = "amRshiny")
  perf_files <- list.files(extdata, pattern = "_ML_perf\\.parquet$", recursive = TRUE, full.names = TRUE)
  skip_if(length(perf_files) == 0, "No demo parquet files found")

  result <- .read_parquet_safe(perf_files[1], verbose = FALSE)
  expect_s3_class(result, "tbl_df")
  expect_gt(nrow(result), 0)
})

# ── .prep_nmcc_data ──────────────────────────────────────────────────────────

test_that(".prep_nmcc_data returns NULL for empty/NULL input", {
  expect_null(.prep_nmcc_data(NULL))
  expect_null(.prep_nmcc_data(data.frame()))
})

test_that(".prep_nmcc_data filters to baseline non-cross rows", {
  df <- tibble::tibble(
    species = c("Sau", "Sau", "Sau"),
    strat_label = c(NA_character_, "country", NA_character_),
    cross_test = c(FALSE, FALSE, TRUE),
    nmcc = c(0.8, 0.7, 0.6),
    species_label = c("Staphylococcus_aureus", "Staphylococcus_aureus", "Staphylococcus_aureus")
  )
  result <- .prep_nmcc_data(df)
  # Should only keep rows where strat_label is NA and cross_test is FALSE

  expect_equal(nrow(result), 1)
  expect_equal(result$nmcc, 0.8)
  expect_true("species_display" %in% names(result))
  expect_equal(result$species_display, "Staphylococcus aureus")
})

test_that(".prep_nmcc_data falls back to species when species_label missing", {
  df <- tibble::tibble(
    species = c("Sau"),
    strat_label = c(NA_character_),
    cross_test = c(FALSE),
    nmcc = c(0.8)
  )
  result <- .prep_nmcc_data(df)
  expect_equal(result$species_display, "Sau")
})
