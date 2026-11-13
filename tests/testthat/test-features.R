## Tests for the Feature Importance tab (R/features_plots.R, R/features_ui.R)

test_that("feature heatmap is relative to each group's top feature", {
  top <- suppressMessages(loadTopFeat(verbose = FALSE))
  w <- makeFeatureImportancePlot(
    top, "Sso", "AMP", "genes", "binary", 100, "across_drug"
  )
  trace <- plotly::plotly_build(w)$x$data[[1]]
  z <- unlist(trace$z)
  z <- z[!is.na(z)]
  expect_equal(max(z), 1)
  # The least important feature shown is still above 0, and colors span 0-1.
  expect_gt(min(z), 0)
  expect_equal(c(trace$zmin, trace$zmax), c(0, 1))
})

test_that("feature tables scroll inside their boxes", {
  html <- as.character(featureImportanceUI("features"))
  for (id in c("across_bug_feature_importance_table", "across_drug_feature_importance_table")) {
    box <- regmatches(html, regexpr(paste0('<div[^>]*overflow-y: auto[^>]*>\\s*<div[^>]*id="features-', id), html))
    expect_length(box, 1)
  }
})
