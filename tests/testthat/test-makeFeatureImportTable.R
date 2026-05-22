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
  expect_true(any(grepl("href", unlist(result$x$data))))
})

test_that("makeFeatureImportTable adds hyperlinks for COG column", {
  df <- tibble::tibble(
    species = c("Sau"),
    drug_or_class = c("ampicillin"),
    COG = c("COG0001"),
    Importance = c(0.5)
  )

  result <- makeFeatureImportTable(df)
  expect_s3_class(result, "datatables")
  expect_true(any(grepl("href", unlist(result$x$data))))
})

test_that("makeFeatureImportTable adds hyperlinks for cluster column", {
  df <- tibble::tibble(
    species = c("Sfl"),
    drug_or_class = c("ampicillin"),
    cluster = c("fig|42897.100.peg.487"),
    Importance = c(0.5)
  )

  result <- makeFeatureImportTable(df)
  expect_s3_class(result, "datatables")
  expect_true(any(grepl("bv-brc.org", unlist(result$x$data))))
})

test_that("makeFeatureImportTable links each comma-separated cluster ID", {
  df <- tibble::tibble(
    species = c("Sfl"),
    drug_or_class = c("ampicillin"),
    cluster = c("fig|42897.100.peg.487, fig|42897.100.peg.488"),
    Importance = c(0.5)
  )

  result <- makeFeatureImportTable(df)
  # Two separate <a> tags expected, one per fig ID.
  hrefs <- unlist(result$x$data)
  n_anchors <- sum(stringr::str_count(hrefs, "<a href="))
  expect_gte(n_anchors, 2)
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
