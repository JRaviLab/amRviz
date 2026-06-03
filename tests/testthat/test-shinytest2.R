# shinytest2 app tests
#
# These open the real dashboard in a hidden browser and check that the app
# starts up and that every tab shows its plots using the built-in Shigella
# flexneri demo data, with no errors. The other test files check the helper
# functions on their own; these check the whole app running together.
#
# These tests only run on a developer's machine. They are skipped on CRAN, on
# continuous integration (GitHub Actions), and whenever shinytest2 or a Chrome
# browser isn't installed. Skipping on CI matters because the browser leaves
# temporary files behind that make R CMD check fail with a "detritus in the
# temp directory" warning.

skip_if_no_shinytest2 <- function() {
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  testthat::skip_if_not_installed("shinytest2")
  chrome <- tryCatch(chromote::find_chrome(), error = function(e) NA_character_)
  if (length(chrome) != 1 || is.na(chrome) || !nzchar(chrome)) {
    testthat::skip("No Chrome/Chromium binary available for shinytest2")
  }
}

# Start the demo app for testing. We allow long wait times because the app
# loads several data files and draws interactive plots when it first opens.
# Setting a fixed seed makes the results the same each time the test is run.
new_demo_app <- function(name) {
  shinytest2::AppDriver$new(
    launchAMRDashboard(),
    name = name,
    seed = 1234,
    timeout = 60 * 1000,
    load_timeout = 60 * 1000
  )
}

# Check that a given plot or table appeared on the page and is not showing
# an error message.
expect_output_ok <- function(app, output_id) {
  html <- app$get_html(sprintf("#%s", output_id))
  testthat::expect_true(
    !is.null(html) && nzchar(html),
    info = sprintf("output '%s' did not appear on the page", output_id)
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

  # Each entry pairs a tab with one plot or table that should appear on it.
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
