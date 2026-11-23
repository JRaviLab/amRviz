## Tests for plot generation functions

# ── makeModelPerformancePlot ─────────────────────────────────────────────────

test_that("makeModelPerformancePlot returns plotly for NULL data", {
  result <- makeModelPerformancePlot(
    data = NULL, bug = "Sau", model_scale = "genes",
    data_type = "binary", metrics = "mcc",
    amr_drug_class = NULL, amr_drug = NULL
  )
  expect_s3_class(result, "plotly")
})

test_that("makeModelPerformancePlot returns plotly for zero-row data", {
  result <- makeModelPerformancePlot(
    data = data.frame(), bug = "Sau", model_scale = "genes",
    data_type = "binary", metrics = "mcc",
    amr_drug_class = NULL, amr_drug = NULL
  )
  expect_s3_class(result, "plotly")
})

test_that("makeModelPerformancePlot generates valid plot with demo data", {
  perf <- suppressMessages(loadMLResults(verbose = FALSE))
  skip_if(nrow(perf) == 0, "No demo performance data")

  species <- unique(perf$species)[1]
  scales <- unique(perf$feature_type)
  subtypes <- unique(perf$feature_subtype)

  result <- makeModelPerformancePlot(
    data = perf,
    bug = species,
    model_scale = scales,
    data_type = subtypes,
    metrics = "mcc",
    amr_drug_class = "all",
    amr_drug = NULL
  )
  expect_s3_class(result, "plotly")
})

# ── makeFeatureImportancePlot ────────────────────────────────────────────────

test_that("makeFeatureImportancePlot returns NULL for empty data", {
  result <- makeFeatureImportancePlot(
    data = NULL, bug = "Sau", amr_drug = "ampicillin",
    model_scale = "genes", data_type_ = "binary",
    top_n_features = 10,
    feature_importance_tabset = "across_bug"
  )
  expect_null(result)
})

# ── makeCrossModelPerformancePlot ────────────────────────────────────────────

test_that("makeCrossModelPerformancePlot returns placeholder for empty data", {
  result <- makeCrossModelPerformancePlot(
    perf_data = NULL, bug = "Sau",
    drug = "ampicillin", cross_model = "country"
  )
  expect_s3_class(result, "plotly")
})

test_that("makeCrossModelPerformancePlot returns placeholder for zero-row data", {
  result <- makeCrossModelPerformancePlot(
    perf_data = data.frame(), bug = "Sau",
    drug = "ampicillin", cross_model = "country"
  )
  expect_s3_class(result, "plotly")
})
