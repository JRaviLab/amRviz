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

# Top features of AMG country models. AUS has two genes and one ARG; GBR has
# one gene and one feature from a cross-tested model. NEG points to
# Resistant, POS to Susceptible.
holdout_top <- function() {
  tibble::tibble(
    task_id = "t1", shuffled = FALSE, seed = 123L,
    species = "Sso", species_label = "Shigella sonnei",
    drug_label = "drug_class", drug_or_class = "AMG", strat_label = "country",
    strat_value = c("AUS", "AUS", "AUS", "GBR", "GBR"),
    strat_value_test = NA_character_,
    cross_test = c(FALSE, FALSE, FALSE, FALSE, TRUE),
    Variable = c("tetA", "tetA", "blaTEM", "tetA", "cross_only"),
    feature_type = c("genes", "args", "genes", "genes", "genes"),
    Importance = c(0.2, 0.4, 0.6, 0.5, 0.9),
    Sign = c("NEG", "POS", "NEG", "NEG", "NEG"),
    metadata_path = "/some/metadata.parquet"
  )
}

test_that("holdoutFeatureTable scales the top features within each group and type", {
  top <- holdout_top()
  tbl <- holdoutFeatureTable(top, "Sso", "AMG", "country", top_n = 10)
  expect_named(tbl, c("feature", "feature_type", "group", "importance", "direction"))
  # The same name in two feature types stays two features.
  expect_equal(sum(tbl$feature == "tetA" & tbl$group == "AUS"), 2)
  # Cross-tested models are left out.
  expect_false("cross_only" %in% tbl$feature)
  # Relative to the top feature of the same group and type, so the least
  # important shown stays above 0.
  aus_genes <- tbl[tbl$group == "AUS" & tbl$feature_type == "genes", ]
  expect_equal(sort(aus_genes$importance), c(0.2, 0.6) / 0.6)
  expect_equal(tbl$importance[tbl$group == "AUS" & tbl$feature_type == "args"], 1)
  # top_n applies per group and feature type: AUS genes, AUS args, GBR genes.
  expect_equal(nrow(holdoutFeatureTable(top, "Sso", "AMG", "country", top_n = 1)), 3)
})

test_that("holdoutFeatureTable keeps each feature's direction", {
  tbl <- holdoutFeatureTable(holdout_top(), "Sso", "AMG", "country")
  expect_equal(tbl$direction[tbl$feature_type == "args"], "Susceptible")
  expect_true(all(tbl$direction[tbl$feature_type == "genes"] == "Resistant"))
  # Seeds that disagree on direction count against each other: the mean of
  # +0.6 and -0.2.
  two_seeds <- holdout_top()[c(3, 3), ]
  two_seeds$Sign <- c("NEG", "POS")
  two_seeds$Importance <- c(0.6, 0.2)
  one <- holdoutFeatureTable(two_seeds, "Sso", "AMG", "country")
  expect_equal(one$direction, "Resistant")
  # Without Sign there is no direction.
  no_sign <- holdoutFeatureTable(
    dplyr::select(holdout_top(), -"Sign"), "Sso", "AMG", "country"
  )
  expect_true(all(is.na(no_sign$direction)))
})

test_that("the feature heatmap colors by direction", {
  p <- plotly::plotly_build(makeCrossModelFeatureImportancePlot(
    holdout_top(), "Sso", "AMG", "country", 5, names = list()
  ))
  z <- unlist(p$x$data[[1]]$z)
  expect_true(any(z < 0, na.rm = TRUE))
  expect_true(any(z > 0, na.rm = TRUE))
})

test_that("holdoutFeatureRows drops the task and constant columns", {
  rows <- holdoutFeatureRows(holdout_top(), "Sso", "AMG", "country")
  expect_equal(nrow(rows), 4)
  expect_false(any(c(
    "task_id", "shuffled", "cross_test", "strat_value_test",
    "species_label", "metadata_path", "Sign"
  ) %in% names(rows)))
  expect_setequal(rows$direction, c("Resistant", "Susceptible"))
})

test_that("makeFeatureImportTable keeps whole numbers such as seeds whole", {
  dt <- makeFeatureImportTable(holdoutFeatureRows(holdout_top(), "Sso", "AMG", "country"))
  expect_equal(unique(as.character(dt$x$data$seed)), "123")
  expect_match(dt$x$data$Importance[[1]], "e-01$")
})

test_that("holdoutNames reads code-name pairs from the ORB metadata", {
  path <- tempfile(fileext = ".parquet")
  arrow::write_parquet(tibble::tibble(
    genome_drug.antibiotic = c("gentamicin", "ampicillin"),
    drug_abbr = c("GEN", "AMP"),
    drug_class = c("aminoglycosides", "penicillins"),
    class_abbr = c("AMG", "PEN"),
    genome.isolation_country = c("Australia", NA),
    country_abbr = c("AUS", NA)
  ), path)
  names <- holdoutNames(tibble::tibble(metadata_path = path))
  expect_equal(names$drug[["GEN"]], "gentamicin")
  expect_equal(names$drug[["AMG"]], "aminoglycosides")
  expect_equal(names$country, c(AUS = "Australia"))
  expect_equal(.holdout_label(c("AUS", "XYZ"), names$country), c("Australia (AUS)", "XYZ"))
  # No metadata, no names.
  expect_length(holdoutNames(holdout_perf())$drug, 0)
})

test_that("holdoutDrugMenu groups drug classes and drugs, with full names", {
  perf <- dplyr::bind_rows(
    holdout_perf(),
    dplyr::mutate(holdout_perf(), drug_label = "drug", drug_or_class = "GEN")
  )
  names <- list(drug = c(GEN = "gentamicin", AMG = "aminoglycosides"))
  menu <- holdoutDrugMenu(perf, "Sso", names)
  expect_named(menu, c("Drug classes", "Drugs"))
  expect_equal(menu$`Drug classes`, list(`aminoglycosides (AMG)` = "AMG"))
  expect_equal(menu$Drugs, list(`gentamicin (GEN)` = "GEN"))
  expect_equal(holdoutDrugMenu(perf, "Kpn", names), character(0))
})

test_that("the grid hover text names each country", {
  names <- list(country = c(AUS = "Australia", GBR = "United Kingdom"))
  p <- plotly::plotly_build(makeCrossModelPerformancePlot(
    holdout_perf(), "Sso", "AMG", "country", names = names
  ))
  hover <- unlist(p$x$data[[1]]$hovertext)
  expect_true(any(grepl("Australia (AUS)", hover, fixed = TRUE)))
})
