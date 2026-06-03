## Tests for makeCrossModelRidgePlot

test_that("makeCrossModelRidgePlot returns plotly placeholder for NULL input", {
  result <- makeCrossModelRidgePlot(NULL, bug = "Sfl", cross_model = "country")
  expect_s3_class(result, "plotly")
})

test_that("makeCrossModelRidgePlot returns plotly placeholder for zero-row data", {
  result <- makeCrossModelRidgePlot(
    tibble::tibble(),
    bug = "Sfl", cross_model = "country"
  )
  expect_s3_class(result, "plotly")
})

test_that("makeCrossModelRidgePlot returns plotly when no rows match selection", {
  df <- tibble::tibble(
    species = c("Eco", "Kpn"),
    strat_label = c("country", "year"),
    drug_label = c("drug_class", "drug_class"),
    drug_or_class = c("FLQ", "FLQ"),
    cross_test = c(TRUE, FALSE),
    bal_acc = c(0.8, 0.7)
  )
  result <- makeCrossModelRidgePlot(df, bug = "Sfl", cross_model = "country")
  expect_s3_class(result, "plotly")
})

test_that("makeCrossModelRidgePlot returns plotly for valid country selection", {
  perf <- suppressMessages(loadMLResults(verbose = FALSE))
  skip_if(nrow(perf) == 0, "No demo performance data")

  bug <- unique(perf$species)[1]
  result <- makeCrossModelRidgePlot(perf, bug = bug, cross_model = "country")
  expect_s3_class(result, "plotly")
})

test_that("makeCrossModelRidgePlot returns plotly for valid year selection", {
  perf <- suppressMessages(loadMLResults(verbose = FALSE))
  skip_if(nrow(perf) == 0, "No demo performance data")

  bug <- unique(perf$species)[1]
  result <- makeCrossModelRidgePlot(perf, bug = bug, cross_model = "year")
  expect_s3_class(result, "plotly")
})
