## Tests for plot generation functions

# ── makeDatAvailabilityPlot ──────────────────────────────────────────────────

test_that("makeDatAvailabilityPlot returns a plotly object", {
  df <- tibble::tibble(
    genome_drug.antibiotic = c(
      "ampicillin", "ampicillin", "tetracycline"
    ),
    genome_drug.resistant_phenotype = c(
      "Resistant", "Susceptible", "Resistant"
    )
  )

  result <- makeDatAvailabilityPlot(df)
  expect_s3_class(result, "plotly")
})

# ── makeTimeSeriesAMRPlot ────────────────────────────────────────────────────

test_that("makeTimeSeriesAMRPlot returns a plotly object", {
  df <- tibble::tibble(
    genome.collection_year = c(2018, 2019, 2020, 2018, 2019, 2020),
    genome_drug.resistant_phenotype = rep(
      c("Resistant", "Susceptible"),
      each = 3
    ),
    n = c(10, 15, 20, 30, 25, 35)
  )

  result <- makeTimeSeriesAMRPlot(df, "ampicillin")
  expect_s3_class(result, "plotly")
})

test_that("makeTimeSeriesAMRPlot handles amr_drug = 'all'", {
  df <- tibble::tibble(
    genome.collection_year = c(2019, 2020),
    genome_drug.resistant_phenotype = c("Resistant", "Susceptible"),
    n = c(100, 200)
  )

  result <- makeTimeSeriesAMRPlot(df, "all")
  expect_s3_class(result, "plotly")
})

# ── makeHostIsolatePlot ──────────────────────────────────────────────────────

test_that("makeHostIsolatePlot returns a plotly object", {
  df <- tibble::tibble(
    genome_drug.antibiotic = c(
      "ampicillin", "ampicillin", "tetracycline"
    ),
    genome.host_common_name = c("Human", "Bovine", "Human"),
    genome.isolation_source = c("blood", "feces", "urine")
  )

  result <- makeHostIsolatePlot(df)
  expect_s3_class(result, "plotly")
})

# ── makeIsolationSourcesPlot ─────────────────────────────────────────────────

test_that("makeIsolationSourcesPlot returns a plotly object", {
  df <- tibble::tibble(
    genome_drug.antibiotic = c(
      "ampicillin", "tetracycline", "ampicillin"
    ),
    genome.isolation_source = c("Blood", "Urine", "Stool")
  )

  result <- makeIsolationSourcesPlot(df)
  expect_s3_class(result, "plotly")
})

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
    model_scale = "genes", data_type = "binary",
    top_n_features = 10,
    feature_importance_tabset = "across_bug"
  )
  expect_null(result)
})

# ── makeCrossModelPerformancePlot ────────────────────────────────────────────

test_that("makeCrossModelPerformancePlot returns NULL for empty data", {
  result <- makeCrossModelPerformancePlot(
    perf_data = NULL, bug = "Sau",
    drug = "ampicillin", cross_model = "country"
  )
  expect_null(result)
})

test_that("makeCrossModelPerformancePlot returns NULL for zero-row data", {
  result <- makeCrossModelPerformancePlot(
    perf_data = data.frame(), bug = "Sau",
    drug = "ampicillin", cross_model = "country"
  )
  expect_null(result)
})
