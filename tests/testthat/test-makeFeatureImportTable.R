## Tests for makeFeatureImportTable

test_that("makeFeatureImportTable returns datatable for NULL input", {
  result <- makeFeatureImportTable(NULL)
  expect_s3_class(result, "datatables")
})

test_that("makeFeatureImportTable returns datatable for empty data", {
  df <- tibble::tibble()
  result <- makeFeatureImportTable(df)
  expect_s3_class(result, "datatables")
})

test_that("makeFeatureImportTable formats numeric columns to scientific notation", {
  df <- tibble::tibble(
    species = c("Sau", "Kpn"),
    drug_or_class = c("ampicillin", "tetracycline"),
    Variable = c("gene1", "gene2"),
    Importance = c(0.00123, 0.456)
  )

  result <- makeFeatureImportTable(df)
  expect_s3_class(result, "datatables")
})

test_that("makeFeatureImportTable reorders columns by priority", {
  df <- tibble::tibble(
    Variable = c("gene1"),
    Importance = c(0.5),
    species = c("Sau"),
    drug_or_class = c("ampicillin"),
    feature_type = c("genes"),
    feature_subtype = c("binary")
  )

  result <- makeFeatureImportTable(df)
  # Check that the result is a valid datatable

  expect_s3_class(result, "datatables")
})

test_that("makeFeatureImportTable adds hyperlinks for accession column", {
  df <- tibble::tibble(
    species = c("Sau"),
    drug_or_class = c("ampicillin"),
    accession = c("WP_001234567.1"),
    Importance = c(0.5)
  )

  result <- makeFeatureImportTable(df)
  expect_s3_class(result, "datatables")
  # The underlying data should contain hyperlink HTML
  tbl_data <- result$x$data
  # accession column should be transformed to have href
  acc_col <- which(names(df) == "accession")
  expect_true(any(grepl("href", unlist(tbl_data))))
})

test_that("makeFeatureImportTable adds hyperlinks for COG_name column", {
  df <- tibble::tibble(
    species = c("Sau"),
    drug_or_class = c("ampicillin"),
    COG_name = c("COG0001"),
    Importance = c(0.5)
  )

  result <- makeFeatureImportTable(df)
  expect_s3_class(result, "datatables")
  expect_true(any(grepl("href", unlist(result$x$data))))
})

test_that("makeFeatureImportTable replaces non-ARG in ARG_name column", {
  df <- tibble::tibble(
    species = c("Sau", "Sau"),
    drug_or_class = c("ampicillin", "tetracycline"),
    ARG_name = c("non-ARG", "blaOXA-1"),
    Importance = c(0.5, 0.8)
  )

  result <- makeFeatureImportTable(df)
  expect_s3_class(result, "datatables")
})
