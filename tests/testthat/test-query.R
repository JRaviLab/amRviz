## Tests for the Query Data tab (R/query_server.R).

query_rows <- function() {
  tibble::tibble(
    species_label = "Shigella sonnei", drug_or_class = "AMP", mcc = 0.9,
    metadata_path = "/Users/someone/orb/metadata.parquet"
  )
}

test_that("queryColumns hides where files sit on the computer", {
  expect_equal(queryColumns(query_rows()), c("species_label", "drug_or_class", "mcc"))
})

test_that("queryTable keeps the chosen columns that exist, in order", {
  expect_named(queryTable(query_rows(), c("mcc", "species_label")), c("mcc", "species_label"))
  expect_named(queryTable(query_rows(), c("mcc", "not_a_column")), "mcc")
  # No choice shows every column that can be shown.
  expect_named(queryTable(query_rows()), c("species_label", "drug_or_class", "mcc"))
  expect_named(queryTable(query_rows(), character(0)), c("species_label", "drug_or_class", "mcc"))
})

test_that("the downloads contain the chosen columns and no file paths", {
  core <- list(queryData = query_rows, topFeatures = query_rows)
  shiny::testServer(queryDataServer, args = list(core = core), {
    session$setInputs(query_data_columns = c("drug_or_class", "mcc"))
    csv <- utils::read.csv(output$query_data_download)
    expect_named(csv, c("drug_or_class", "mcc"))

    session$setInputs(top_features_columns = NULL)
    csv <- utils::read.csv(output$top_features_download)
    expect_false("metadata_path" %in% names(csv))
  })
})
