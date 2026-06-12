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
  ui <- networkUI()
  expect_s3_class(ui, "shiny.tag")
})
