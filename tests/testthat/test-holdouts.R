## Tests for the Model holdouts tab's tables (R/holdouts_plots.R).

# Country models for one drug class: trained on AUS and tested on AUS
# ("Same") or GBR ("Different"), and trained on GBR and tested on AUS, over two
# seeds. One model has no score.
holdout_perf <- function() {
  tibble::tibble(
    species = "Sso",
    drug_label = "drug_class",
    drug_or_class = "AMG",
    strat_label = "country",
    strat_value = c("AUS", "AUS", "AUS", "AUS", "GBR", "GBR"),
    strat_value_test = c(NA, NA, "GBR", "GBR", "AUS", "AUS"),
    cross_test = c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE),
    bal_acc = c(0.9, 1.0, 0.6, 0.7, 0.8, NA)
  )
}

test_that("holdoutAccuracyTable labels each model by its test group", {
  tbl <- holdoutAccuracyTable(holdout_perf(), "Sso", "country")
  expect_named(tbl, c("drug_class", "test", "value"))
  expect_equal(sum(tbl$test == "Same"), 2)
  # The model without a score is left out.
  expect_equal(nrow(tbl), 5)
  expect_equal(nrow(holdoutAccuracyTable(holdout_perf(), "Sso", "time")), 0)
})

test_that("holdoutGridTable averages each train and test pair", {
  tbl <- holdoutGridTable(holdout_perf(), "Sso", "AMG", "country")
  expect_named(tbl, c("train", "test", "value", "n_models"))
  same <- tbl[tbl$train == "AUS" & tbl$test == "AUS", ]
  expect_equal(c(same$value, same$n_models), c(0.95, 2))
  expect_equal(tbl$value[tbl$train == "GBR"], 0.8)
  expect_equal(nrow(holdoutGridTable(holdout_perf(), "Sso", "TET", "country")), 0)
})

test_that("the grid leaves a cell blank where there is no model", {
  p <- plotly::plotly_build(makeCrossModelPerformancePlot(holdout_perf(), "Sso", "AMG", "country"))
  trace <- p$x$data[[1]]
  # GBR was never tested on itself.
  expect_true("" %in% unlist(trace$text))
})

test_that("holdoutFeatureTable scales the top features within each group", {
  top <- tibble::tibble(
    species = "Sso", drug_or_class = "AMG", strat_label = "country",
    strat_value = c("AUS", "AUS", "AUS", "GBR", "GBR"),
    cross_test = c(FALSE, FALSE, FALSE, FALSE, TRUE),
    Variable = c("tetA", "tetA", "blaTEM", "tetA", "cross_only"),
    feature_type = c("genes", "args", "genes", "genes", "genes"),
    Importance = c(0.2, 0.4, 0.6, 0.5, 0.9)
  )
  tbl <- holdoutFeatureTable(top, "Sso", "AMG", "country", top_n = 10)
  expect_named(tbl, c("feature", "feature_type", "group", "importance"))
  # The same name in two feature types stays two features.
  expect_equal(sum(tbl$feature == "tetA" & tbl$group == "AUS"), 2)
  # Cross-tested models are left out.
  expect_false("cross_only" %in% tbl$feature)
  # Relative to the group's top feature, so the least important shown stays
  # above 0.
  expect_equal(sort(tbl$importance[tbl$group == "AUS"]), c(0.2, 0.4, 0.6) / 0.6)
  expect_equal(nrow(holdoutFeatureTable(top, "Sso", "AMG", "country", top_n = 1)), 2)
})
