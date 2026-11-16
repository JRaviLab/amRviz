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

test_that("across-bug heatmap keeps every species, ESKAPE species first", {
  top <- suppressMessages(loadTopFeat(verbose = FALSE))
  sso <- top[top$species == "Sso", ]
  mixed <- rbind(
    dplyr::mutate(sso, species = "Sso"),
    dplyr::mutate(sso, species = "Kpn"),
    dplyr::mutate(sso, species = "Aba")
  )
  drug <- sso$drug_or_class[sso$feature_type == "genes"][1]
  w <- makeFeatureImportancePlot(
    mixed, c("Sso", "Kpn", "Aba"), drug, "genes", "binary", 10, "across_bug"
  )
  expect_equal(as.character(plotly::plotly_build(w)$x$data[[1]]$x), c("Kpn", "Aba", "Sso"))
})

test_that("feature heatmap shows only the selected model scale", {
  top <- suppressMessages(loadTopFeat(verbose = FALSE))
  sso <- top[top$species == "Sso" & top$feature_type == "genes", ]
  drug <- sso$drug_or_class[1]
  struct <- dplyr::mutate(sso, feature_type = "struct", Variable = paste0("struct_", Variable))
  both <- rbind(sso, struct)
  genes <- makeFeatureImportancePlot(both, "Sso", drug, "genes", "binary", "all", "across_drug")
  features <- unlist(plotly::plotly_build(genes)$x$data[[1]]$y)
  expect_false(any(startsWith(features, "struct_")))
  structural <- makeFeatureImportancePlot(both, "Sso", drug, "struct", "binary", "all", "across_drug")
  expect_true(all(startsWith(unlist(plotly::plotly_build(structural)$x$data[[1]]$y), "struct_")))
})

test_that("model scale menu offers structural features", {
  html <- as.character(featureImportanceUI("features"))
  expect_match(html, '<option value="struct"', fixed = TRUE)
})
