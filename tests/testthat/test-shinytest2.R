# shinytest2 smoke tests
#
# These drive the real dashboard in a headless browser to confirm that the app
# launches and that every tab renders against the bundled Shigella flexneri
# demo data without throwing a Shiny output error. They complement the unit
# tests in the other test files, which exercise the internal data/plot helpers
# in isolation.
#
# The tests are skipped automatically when shinytest2 or a Chrome/Chromium
# binary is unavailable (e.g. on some build machines), so they never turn into
# spurious errors during R CMD check / BiocCheck.

skip_if_no_shinytest2 <- function() {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("shinytest2")
  chrome <- tryCatch(chromote::find_chrome(), error = function(e) NA_character_)
  if (length(chrome) != 1 || is.na(chrome) || !nzchar(chrome)) {
    testthat::skip("No Chrome/Chromium binary available for shinytest2")
  }
}

# Build an AppDriver for the demo app with generous load/idle timeouts, since
# the first render reads several Parquet files and draws plotly/networkD3
# widgets. A fixed seed keeps any stochastic layout reproducible.
new_demo_app <- function(name) {
  shinytest2::AppDriver$new(
    launchAMRDashboard(),
    name = name,
    seed = 1234,
    timeout = 60 * 1000,
    load_timeout = 60 * 1000
  )
}

# Assert that a given Shiny output rendered and is not displaying an error.
expect_output_ok <- function(app, output_id) {
  html <- app$get_html(sprintf("#%s", output_id))
  testthat::expect_true(
    !is.null(html) && nzchar(html),
    info = sprintf("output '%s' is missing from the DOM", output_id)
  )
  testthat::expect_false(
    grepl("shiny-output-error", html, fixed = TRUE),
    info = sprintf("output '%s' rendered a Shiny error", output_id)
  )
}

test_that("dashboard launches on the home tab with no output errors", {
  skip_if_no_shinytest2()

  app <- new_demo_app("amRviz-launch")
  withr::defer(app$stop())
  app$wait_for_idle(timeout = 60 * 1000)

  # The app starts on the Home tab.
  expect_identical(app$get_value(input = "tabselected"), "home")

  # Nothing on the initial page should be in an errored state.
  errors <- app$get_html(".shiny-output-error")
  expect_true(is.null(errors) || !nzchar(errors))
})

test_that("each dashboard tab renders its primary output", {
  skip_if_no_shinytest2()

  app <- new_demo_app("amRviz-tabs")
  withr::defer(app$stop())
  app$wait_for_idle(timeout = 60 * 1000)

  # tab value (input$tabselected) -> a primary output that must render on it.
  tabs <- list(
    list(value = "Metadata", output = "quick_metadata_stats"),
    list(value = "model_perf_tab", output = "model_perfomance_plot"),
    list(
      value = "Bug/Drug feature comparison",
      output = "across_bug_feature_importance_plot"
    ),
    list(value = "Model holdouts", output = "cross_model_ridge_country"),
    list(value = "Network", output = "drug_feature_network"),
    list(value = "query_datatable", output = "queryDataTable")
  )

  for (tab in tabs) {
    app$set_inputs(tabselected = tab$value)
    app$wait_for_idle(timeout = 60 * 1000)
    expect_identical(app$get_value(input = "tabselected"), tab$value)
    expect_output_ok(app, tab$output)
  }
})
