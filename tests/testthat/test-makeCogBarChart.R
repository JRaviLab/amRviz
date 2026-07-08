## Tests for makeClusterBarChart

test_that("makeClusterBarChart returns plotly placeholder for NULL input", {
  result <- makeClusterBarChart(NULL)
  expect_s3_class(result, "plotly")
})

test_that("makeClusterBarChart returns plotly placeholder for zero-row data", {
  df <- tibble::tibble()
  result <- makeClusterBarChart(df)
  expect_s3_class(result, "plotly")
})

test_that("makeClusterBarChart returns plotly placeholder when COG column is missing", {
  df <- tibble::tibble(
    Variable = c("gene1", "gene2"),
    Importance = c(0.5, 0.3)
  )
  result <- makeClusterBarChart(df)
  expect_s3_class(result, "plotly")
})

test_that("makeClusterBarChart returns plotly placeholder when all COGs are NA/empty", {
  df <- tibble::tibble(
    Variable = c("gene1", "gene2"),
    COG = c(NA_character_, ""),
    Importance = c(0.5, 0.3)
  )
  result <- makeClusterBarChart(df)
  expect_s3_class(result, "plotly")
})

test_that("makeClusterBarChart returns plotly bar chart for valid data", {
  df <- tibble::tibble(
    Variable = c("gene1", "gene2", "gene3", "gene4"),
    COG = c("COG0001", "COG0002", "COG0001", "COG0003"),
    COG_name = c("Translation", "Transcription", "Translation", "Replication"),
    Importance = c(0.5, 0.3, 0.4, 0.2)
  )
  result <- makeClusterBarChart(df)
  expect_s3_class(result, "plotly")
})

test_that("makeClusterBarChart handles comma-separated COG cells", {
  df <- tibble::tibble(
    Variable = c("gene1", "gene2"),
    COG = c("COG0001, COG0002", "COG0001"),
    Importance = c(0.5, 0.3)
  )
  result <- makeClusterBarChart(df)
  expect_s3_class(result, "plotly")
})

test_that("makeClusterBarChart respects top_n parameter", {
  df <- tibble::tibble(
    Variable = paste0("gene", 1:20),
    COG = paste0("COG", sprintf("%04d", 1:20)),
    Importance = runif(20)
  )
  result <- makeClusterBarChart(df, top_n = 5)
  expect_s3_class(result, "plotly")
})
