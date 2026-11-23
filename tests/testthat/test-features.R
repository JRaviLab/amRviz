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

# ── table under the heatmap (review on #66) ─────────────────────────────────

test_that("the feature table holds exactly the features the heatmap shows", {
  top <- suppressMessages(loadTopFeat(verbose = FALSE))
  sso <- top[top$species == "Sso", ]
  drug <- "AMP"
  shown <- function(n, scale = "genes", type = "binary") {
    featureImportanceTable(sso, "Sso", drug, scale, type, n, "across_drug")
  }
  rows <- function(n, scale = "genes", type = "binary") {
    featureImportanceRows(sso, "Sso", drug, scale, type, n, "across_drug", names = list())
  }
  # The table follows the slider ...
  expect_setequal(unique(rows(5)$Variable), shown(5)$cluster_name)
  expect_lt(length(unique(rows(5)$Variable)), length(unique(rows(20)$Variable)))
  # ... the model scale and the data type.
  expect_true(all(rows(5, "domains")$feature_type == "domains"))
  expect_true(all(rows(5, "genes", "counts")$feature_subtype == "counts"))
  # And the heatmap draws the same features.
  w <- makeFeatureImportancePlot(sso, "Sso", drug, "genes", "binary", 5, "across_drug", names = list())
  expect_setequal(unlist(plotly::plotly_build(w)$x$data[[1]]$y), shown(5)$cluster_name)
})

test_that("the feature table leaves out the task and constant columns", {
  top <- suppressMessages(loadTopFeat(verbose = FALSE))
  rows <- featureImportanceRows(
    top, "Sso", "AMP", "genes", "binary", 5, "across_drug", names = list()
  )
  expect_gt(nrow(rows), 0)
  expect_false(any(c(
    "task_id", "shuffled", "strat_label", "strat_value", "strat_value_test",
    "cross_test", "species_label", "metadata_path"
  ) %in% names(rows)))
  expect_true(all(c("species", "drug_or_class", "Variable", "seed", "Importance") %in% names(rows)))
})

test_that("cross-tested models are left out of the heatmap and table", {
  top <- suppressMessages(loadTopFeat(verbose = FALSE))
  sso <- top[top$species == "Sso" & top$drug_or_class == "AMP" &
    top$feature_type == "genes" & top$feature_subtype == "binary" &
    !top$cross_test %in% TRUE & (is.na(top$strat_label) | !nzchar(top$strat_label)), ]
  cross <- dplyr::mutate(sso[1, ], cross_test = TRUE, Variable = "cross_only", Importance = 1e6)
  both <- rbind(sso, cross)
  expect_false("cross_only" %in%
    featureImportanceTable(both, "Sso", "AMP", "genes", "binary", "all", "across_drug")$cluster_name)
  expect_false("cross_only" %in%
    featureImportanceRows(both, "Sso", "AMP", "genes", "binary", "all", "across_drug", names = list())$Variable)
})

test_that("drugs and classes are named in full in the heatmap and table", {
  top <- suppressMessages(loadTopFeat(verbose = FALSE))
  sso <- top[top$species == "Sso", ]
  names <- amrNames(sso)
  expect_equal(names$drug[["AMP"]], "ampicillin")
  expect_equal(names$drug[["AMG"]], "aminoglycosides")
  w <- makeFeatureImportancePlot(sso, "Sso", c("AMP", "TET"), "genes", "binary", 5, "across_drug", names = names)
  expect_setequal(
    as.character(plotly::plotly_build(w)$x$data[[1]]$x),
    c("ampicillin (AMP)", "tetracycline (TET)")
  )
  rows <- featureImportanceRows(sso, "Sso", "AMP", "genes", "binary", 5, "across_drug", names = names)
  expect_equal(unique(rows$drug_or_class), "ampicillin (AMP)")
  expect_equal(.amr_choices(c("AMG", "XYZ"), names$drug), c(`aminoglycosides (AMG)` = "AMG", XYZ = "XYZ"))
})
