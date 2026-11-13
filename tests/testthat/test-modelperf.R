## Tests for the Model Performance tab's tables (R/modelperf_plots.R). The
## tables are what the plots draw, so these pin down their columns.

# Baseline models for two species, one ESKAPE and one not, plus a stratified
# and a cross-tested row that the tables leave out.
perf_rows <- function() {
  tibble::tibble(
    species = c("Sfl", "Sfl", "Sau", "Sau", "Sfl", "Sfl"),
    species_label = c(
      "Shigella flexneri", "Shigella flexneri", "Staphylococcus aureus",
      "Staphylococcus aureus", "Shigella flexneri", "Shigella flexneri"
    ),
    feature_type = c("genes", "DefenseCas", "genes", "domains", "genes", "genes"),
    feature_subtype = "binary",
    drug_label = c("drug", "drug_class", "drug_class", "drug", "drug", "drug"),
    drug_or_class = c("AMP", "BLA", "BLA", "OXA", "AMP", "AMP"),
    strat_label = c(NA, NA, NA, NA, "country", NA),
    cross_test = c(FALSE, NA, FALSE, FALSE, FALSE, TRUE),
    mcc = c(0.9, 0.5, 0.8, 0.7, 0.1, 0.2)
  )
}

test_that("modelPerfTable keeps baseline models and marks the selected drug", {
  tbl <- modelPerfTable(
    perf_rows(),
    bug = c("Sfl", "Sau"), model_scale = c("genes", "domains", "DefenseCas"),
    data_type = "binary", metrics = "mcc", amr_drug = "AMP"
  )
  expect_named(tbl, c(
    "species", "species_display", "feature_type", "feature_subtype",
    "drug_label", "drug_or_class", "value", "highlighted"
  ))
  expect_equal(nrow(tbl), 4)
  expect_equal(tbl$value[tbl$highlighted], 0.9)
  # ESKAPE species first, and other species kept rather than dropped.
  expect_equal(levels(tbl$species_display), c("Staphylococcus aureus", "Shigella flexneri"))
  expect_false(anyNA(tbl$species_display))
  expect_true("DefenseCas" %in% tbl$feature_type)
})

test_that("modelPerfTable returns an empty table when nothing matches", {
  tbl <- modelPerfTable(perf_rows(), "Kpn", "genes", "binary", "mcc")
  expect_equal(nrow(tbl), 0)
  expect_named(tbl, names(modelPerfTable(perf_rows(), "Sfl", "genes", "binary", "mcc")))
  expect_equal(nrow(modelPerfTable(NULL, "Sfl", "genes", "binary", "mcc")), 0)
})

test_that("mccStripTable labels every scale, including Defense/Cas", {
  tbl <- mccStripTable(perf_rows(), selected_drug_class = "BLA")
  expect_equal(nrow(tbl), 4)
  expect_setequal(as.character(tbl$scale_label), c("Gene", "Defense/Cas", "Domain"))
  expect_equal(tbl$mcc[tbl$highlighted], c(0.5, 0.8))
  expect_false(any(mccStripTable(perf_rows(), "all", "AMP")$highlighted))
})

test_that("mccHeatmapTable summarizes drug-class models by each panel", {
  tbl <- mccHeatmapTable(perf_rows())
  expect_named(tbl, c("panel", "column", "drug_or_class", "median_mcc", "n_models"))
  expect_setequal(unique(tbl$panel), c("Species", "Molecular scale", "Data type"))
  species <- tbl[tbl$panel == "Species", ]
  expect_equal(species$median_mcc[species$column == "Shigella flexneri"], 0.5)
  binary <- tbl[tbl$panel == "Data type", ]
  expect_equal(binary$column, "Binary")
  expect_equal(binary$n_models, 2)
  expect_null(mccHeatmapTable(perf_rows()[perf_rows()$drug_label == "drug", ]))
})

test_that("a one-cell heatmap keeps its color instead of turning gray", {
  one_class <- perf_rows()[perf_rows()$species == "Sfl", ]
  p <- .fix_single_value_heatmap(plotly::ggplotly(
    ggplot2::ggplot(
      data.frame(x = "a", y = "b", z = 0.7),
      ggplot2::aes(x, y, fill = z)
    ) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient(limits = c(0, 1))
  ))
  expect_length(p$x$data[[1]]$colorscale, 2)
  expect_s3_class(makeMCCHeatmap(one_class), "plotly")
})

test_that("the scale selector offers only the scales in the results", {
  expect_equal(
    .scale_choices(c("genes", "DefenseCas")),
    c(Gene = "genes", "Defense/Cas" = "DefenseCas")
  )
  expect_true(all(names(SCALE_COLORS) %in% .scale_choices()))
})

test_that("the tab's plots draw the demo data", {
  perf <- suppressMessages(loadMLResults(verbose = FALSE))
  expect_s3_class(makeMCCStripPlot(perf, "all"), "plotly")
  expect_s3_class(makeMCCHeatmap(perf, unique(perf$drug_or_class)[1]), "plotly")
  expect_s3_class(
    makeModelPerformancePlot(
      perf, unique(perf$species), unique(perf$feature_type), "binary", "mcc",
      "all", unique(perf$drug_or_class)[1]
    ),
    "plotly"
  )
})
