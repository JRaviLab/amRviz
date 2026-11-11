## Tests for the Network tab: makeDrugFeatureNetwork, makeFeatureEgoNetwork,
## and networkUI.

# ── makeDrugFeatureNetwork ──────────────────────────────────────────────────

test_that("makeDrugFeatureNetwork returns NULL for NULL input", {
  expect_null(makeDrugFeatureNetwork(NULL, bug = "Sfl"))
})

test_that("makeDrugFeatureNetwork returns NULL for non-data.frame input", {
  expect_null(makeDrugFeatureNetwork("not a data frame", bug = "Sfl"))
})

test_that("makeDrugFeatureNetwork returns NULL for zero-row data", {
  expect_null(makeDrugFeatureNetwork(tibble::tibble(), bug = "Sfl"))
})

test_that("makeDrugFeatureNetwork returns a forceNetwork for demo data", {
  skip_if_not_installed("networkD3")

  top_data <- suppressMessages(loadTopFeat(verbose = FALSE))
  skip_if(nrow(top_data) == 0, "No demo top-features data")

  bug <- unique(top_data$species)[1]
  result <- makeDrugFeatureNetwork(top_data, bug = bug, top_n = 5)
  expect_s3_class(result, "htmlwidget")
})

# ── networkTables ──────────────────────────────────────────────────────────

# Baseline models for two drugs, plus a stratified and a cross-tested row the
# network leaves out. One feature is named like a drug.
network_rows <- function() {
  tibble::tibble(
    species = "Sso",
    drug_or_class = c("AMP", "AMP", "AMP", "TET", "TET", "AMP", "AMP"),
    feature_type = c("domains", "domains", "domains", "genes", "genes", "genes", "genes"),
    Variable = c("DDE_Tnp_Tn3", "DDE_Tnp_1", "PF00001_IPR000001", "tetA", "AMP", "strat", "cross"),
    Importance = c(0.9, 0.8, 0.7, 0.6, 0.5, 2, 3),
    strat_label = c(NA, NA, NA, NA, NA, "country", NA),
    cross_test = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE)
  )
}

test_that("networkTables links drugs to their baseline features", {
  net <- networkTables(network_rows(), "Sso", top_n = 10)
  expect_named(net$nodes, c("id", "label", "group"))
  expect_named(net$links, c("from", "to", "importance"))
  # Stratified and cross-tested models are left out.
  expect_false(any(c("strat", "cross") %in% net$nodes$label))
  expect_true(all(c(net$links$from, net$links$to) %in% net$nodes$id))
})

test_that("networkTables keeps current Pfam names whole", {
  net <- networkTables(network_rows(), "Sso", top_n = 10)
  domains <- net$nodes$label[net$nodes$group == "domains"]
  expect_setequal(domains, c("DDE_Tnp_Tn3", "DDE_Tnp_1", "PF00001"))
})

test_that("networkTables keeps a feature named like a drug as its own node", {
  net <- networkTables(network_rows(), "Sso", top_n = 10)
  expect_equal(sum(net$nodes$label == "AMP"), 2)
  expect_equal(nrow(networkTables(network_rows(), "Sso", top_n = 1)$links), 2)
  expect_null(networkTables(network_rows(), "Sfl"))
})

# ── makeFeatureEgoNetwork ───────────────────────────────────────────────────

test_that("makeFeatureEgoNetwork returns NULL for NULL input", {
  expect_null(makeFeatureEgoNetwork(NULL, variable = "x"))
})

test_that("makeFeatureEgoNetwork returns NULL for zero-row data", {
  expect_null(makeFeatureEgoNetwork(tibble::tibble(), variable = "x"))
})

test_that("makeFeatureEgoNetwork returns NULL when variable is NULL or empty", {
  df <- tibble::tibble(Variable = "x", cluster = "c1", COG = "COG0001")
  expect_null(makeFeatureEgoNetwork(df, variable = NULL))
  expect_null(makeFeatureEgoNetwork(df, variable = ""))
})

test_that("makeFeatureEgoNetwork returns NULL when variable not in table", {
  df <- tibble::tibble(Variable = "x", cluster = "c1", COG = "COG0001")
  expect_null(makeFeatureEgoNetwork(df, variable = "not-here"))
})

test_that("makeFeatureEgoNetwork returns NULL when no clusters or COGs", {
  df <- tibble::tibble(
    Variable = "x",
    cluster = NA_character_,
    COG = NA_character_
  )
  expect_null(makeFeatureEgoNetwork(df, variable = "x"))
})

test_that("makeFeatureEgoNetwork returns forceNetwork for valid annotated row", {
  skip_if_not_installed("networkD3")

  df <- tibble::tibble(
    Variable = "x",
    cluster = "fig|42897.100.peg.487",
    COG = "COG0001, COG0002"
  )
  result <- makeFeatureEgoNetwork(df, variable = "x")
  expect_s3_class(result, "htmlwidget")
})

# ── networkUI ───────────────────────────────────────────────────────────────

test_that("networkUI returns a shiny tabPanel", {
  ui <- networkUI("network")
  expect_s3_class(ui, "shiny.tag")
})

# ── shared font + node-shape styling ─────────────────────────────────────────

test_that("network visualisations carry the shared font + node-shape hook", {
  skip_if_not_installed("networkD3")

  top_data <- suppressMessages(loadTopFeat(verbose = FALSE))
  skip_if(nrow(top_data) == 0, "No demo top-features data")

  net <- makeDrugFeatureNetwork(
    top_data,
    bug = unique(top_data$species)[1], top_n = 5
  )
  # Font matches the plotly/ggplot panels.
  expect_equal(net$x$options$fontFamily, "Arial, sans-serif")
  # Per-group node shapes are applied via an onRender hook.
  expect_false(is.null(net$jsHooks$render))
  expect_match(net$jsHooks$render[[1]]$code, "d3.symbol")
  # Shapes are mirrored in both the nodes and the legend.
  expect_match(net$jsHooks$render[[1]]$code, "node-shape")
  expect_match(net$jsHooks$render[[1]]$code, "legend-shape")
})
