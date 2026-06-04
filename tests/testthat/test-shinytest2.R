# shinytest2 app tests
#
# These open the real dashboard in a hidden browser and check that the app
# starts up and that every tab shows its plots using the built-in Shigella
# flexneri demo data, with no errors. The other test files check the helper
# functions on their own; these check the whole app running together.
#
# By default these run only on a developer's machine. They are skipped on CRAN,
# and whenever shinytest2 or a Chrome browser isn't installed.
#
# On CI they are skipped too, UNLESS the job sets the environment variable
# RUN_APP_TESTS=true. This lets a dedicated app-test job run them (using
# devtools::test() on their own), while the regular R CMD check workflow leaves
# the variable unset and skips them -- otherwise the browser leaves temporary
# files behind that make R CMD check fail with a "detritus in the temp
# directory" warning.

skip_if_no_shinytest2 <- function() {
  testthat::skip_on_cran()
  if (!identical(Sys.getenv("RUN_APP_TESTS"), "true")) {
    testthat::skip_on_ci()
  }
  testthat::skip_if_not_installed("shinytest2")
  chrome <- tryCatch(chromote::find_chrome(), error = function(e) NA_character_)
  if (length(chrome) != 1 || is.na(chrome) || !nzchar(chrome)) {
    testthat::skip("No Chrome/Chromium binary available for shinytest2")
  }
}

# Start the demo app for testing. The app loads several data files
# and draws interactive plots when it first opens.
# Setting a fixed seed makes the results the same each time the test is run.
#
# When it starts the app, shinytest2 scans the server function for global
# variables and warns that it cannot resolve the bare column names used in
# dplyr (non-standard evaluation). This is a harmless false positive -- the app
# runs correctly -- so we quietly drop just that one warning and let any other
# warning through.
new_demo_app <- function(name) {
  withCallingHandlers(
    shinytest2::AppDriver$new(
      launchAMRDashboard(),
      name = name,
      seed = 1234,
      timeout = 60 * 1000,
      load_timeout = 60 * 1000
    ),
    warning = function(w) {
      if (grepl("locate globals|globalsByName", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
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

# Give the browser its own temp directory for the duration of this test file, so
# the scratch files Chrome creates land there and are deleted when the file
# finishes, instead of lingering in the shared session temp dir (which would
# clutter a reviewer's session and can trip R CMD check's "detritus in the temp
# directory" check). Chrome chooses where to write from these OS temp env vars,
# and withr both creates the directory and removes it at teardown.
#
# NOTE: on Windows, deleting these files can fail if the (shared) Chrome process
# still holds them open when teardown runs. Confirming that across Win/Mac/Ubuntu
# is the main purpose of the dedicated app-tests workflow.
.chrome_tmp <- withr::local_tempdir(
  "amrviz-chrome-",
  .local_envir = testthat::teardown_env()
)
withr::local_envvar(
  c(TMPDIR = .chrome_tmp, TMP = .chrome_tmp, TEMP = .chrome_tmp),
  .local_envir = testthat::teardown_env()
)

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
