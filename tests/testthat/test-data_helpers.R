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
  result <- suppressMessages(
    .read_parquet_safe("/nonexistent/path.parquet", verbose = FALSE)
  )
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 0)
})

test_that(".read_parquet_safe returns empty tibble for NULL path", {
  result <- suppressMessages(
    .read_parquet_safe(NULL, verbose = FALSE)
  )
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 0)
})

test_that(".read_parquet_safe reads a valid parquet file", {
  extdata <- system.file("extdata", package = "amRviz")
  perf_files <- list.files(
    extdata,
    pattern = "_perf\\.parquet$",
    recursive = TRUE, full.names = TRUE
  )
  skip_if(length(perf_files) == 0, "No demo parquet files found")

  result <- .read_parquet_safe(perf_files[1], verbose = FALSE)
  expect_s3_class(result, "tbl_df")
  expect_gt(nrow(result), 0)
})
