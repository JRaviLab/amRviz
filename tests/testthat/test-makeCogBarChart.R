## Tests for makeCogBarChart

test_that("makeCogBarChart returns plotly placeholder for NULL input", {
  result <- makeCogBarChart(NULL)
  expect_s3_class(result, "plotly")
})

test_that("makeCogBarChart returns plotly placeholder for zero-row data", {
  df <- tibble::tibble()
  result <- makeCogBarChart(df)
  expect_s3_class(result, "plotly")
})

test_that("makeCogBarChart returns plotly placeholder when COG column is missing", {
  df <- tibble::tibble(
    Variable = c("gene1", "gene2"),
    Importance = c(0.5, 0.3)
  )
  result <- makeCogBarChart(df)
  expect_s3_class(result, "plotly")
})

test_that("makeCogBarChart returns plotly placeholder when all COGs are NA/empty", {
  df <- tibble::tibble(
    Variable = c("gene1", "gene2"),
    COG = c(NA_character_, ""),
    Importance = c(0.5, 0.3)
  )
  result <- makeCogBarChart(df)
  expect_s3_class(result, "plotly")
})

test_that("makeCogBarChart returns plotly bar chart for valid data", {
  df <- tibble::tibble(
    Variable = c("gene1", "gene2", "gene3", "gene4"),
    COG = c("COG0001", "COG0002", "COG0001", "COG0003"),
    COG_name = c("Translation", "Transcription", "Translation", "Replication"),
    Importance = c(0.5, 0.3, 0.4, 0.2)
  )
  result <- makeCogBarChart(df)
  expect_s3_class(result, "plotly")
})

test_that("makeCogBarChart handles comma-separated COG cells", {
  df <- tibble::tibble(
    Variable = c("gene1", "gene2"),
    COG = c("COG0001, COG0002", "COG0001"),
    Importance = c(0.5, 0.3)
  )
  result <- makeCogBarChart(df)
  expect_s3_class(result, "plotly")
})

test_that("makeCogBarChart respects top_n parameter", {
  df <- tibble::tibble(
    Variable = paste0("gene", 1:20),
    COG = paste0("COG", sprintf("%04d", 1:20)),
    Importance = runif(20)
  )
  result <- makeCogBarChart(df, top_n = 5)
  expect_s3_class(result, "plotly")
})

test_that(".clean_cog_name strips standalone NA tokens", {
  expect_equal(.clean_cog_name("Alanine racemase NA NA NA"), "Alanine racemase")
  expect_equal(
    .clean_cog_name("D-serine deaminase NA"), "D-serine deaminase"
  )
  # NA tokens in the middle collapse too
  expect_equal(.clean_cog_name("foo NA bar"), "foo bar")
})

test_that(".clean_cog_name returns NA when nothing meaningful remains", {
  expect_true(is.na(.clean_cog_name("NA")))
  expect_true(is.na(.clean_cog_name("NA NA NA")))
  expect_true(is.na(.clean_cog_name(NA_character_)))
})

test_that(".clean_cog_name strips NA tokens adjacent to semicolons", {
  expect_equal(
    .clean_cog_name("Foo NA NA NA; Bar NA NA"),
    "Foo; Bar"
  )
  expect_equal(
    .clean_cog_name("Foo NA; NA; Bar"),
    "Foo; Bar"
  )
})

test_that(".clean_cog_name leaves clean names untouched and is vectorised", {
  expect_equal(.clean_cog_name("Beta-lactamase class C"), "Beta-lactamase class C")
  expect_equal(
    .clean_cog_name(c("Helicase NA", "NA", "Transcription")),
    c("Helicase", NA, "Transcription")
  )
})
