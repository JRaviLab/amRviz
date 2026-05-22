## Tests for UI module functions

test_that("queryDataUI returns a shiny tabPanel", {
  result <- queryDataUI()
  expect_s3_class(result, "shiny.tag")
  expect_true("shiny.tag" %in% class(result))
})

test_that("metadataUI returns a shiny tabPanel", {
  result <- metadataUI()
  expect_s3_class(result, "shiny.tag")
})

test_that("modelPerfUI returns a shiny tabPanel", {
  result <- modelPerfUI()
  expect_s3_class(result, "shiny.tag")
})

test_that("featureImportanceUI returns a shiny tabPanel", {
  result <- featureImportanceUI()
  expect_s3_class(result, "shiny.tag")
})

test_that("crossModelComparisonUI returns a shiny tabPanel", {
  result <- crossModelComparisonUI()
  expect_s3_class(result, "shiny.tag")
})

# ── UI helper functions ──────────────────────────────────────────────────────

test_that("quickStatBox returns a shiny tag", {
  result <- quickStatBox("Test Title", "42")
  expect_s3_class(result, "shiny.tag")
})

test_that("quickStatBox accepts custom colors", {
  result <- quickStatBox(
    "Title", "100",
    bg_color = "#fff", text_color = "black"
  )
  expect_s3_class(result, "shiny.tag")
})

test_that("amr_select returns a shiny tag with selectInput", {
  result <- amr_select(
    "test_id", "Test Label",
    choices = c("A", "B", "C")
  )
  expect_s3_class(result, "shiny.tag")
})

test_that("amr_select respects multiple parameter", {
  result_multi <- amr_select(
    "id1", "Label", c("A", "B"),
    multiple = TRUE
  )
  result_single <- amr_select(
    "id2", "Label", c("A", "B"),
    multiple = FALSE
  )
  expect_s3_class(result_multi, "shiny.tag")
  expect_s3_class(result_single, "shiny.tag")
})

test_that("amr_button returns a shiny tag with actionButton", {
  result <- amr_button(
    "btn_id", "Click Me", "download", "btn-primary"
  )
  expect_s3_class(result, "shiny.tag")
})

test_that("launchAMRDashboard returns a shiny.appobj", {
  app <- launchAMRDashboard()
  expect_s3_class(app, "shiny.appobj")
})
