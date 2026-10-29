# Each dashboard tab run on its own with shiny::testServer(), using the
# packaged Shigella demo data. The Metadata tab is tested in test-metadata.R.

# Shared data in the shape setupServerCore() returns, built from the demo data
# without starting the app.
demo_core <- function() {
  perf <- suppressMessages(loadMLResults(verbose = FALSE))
  top <- suppressMessages(loadTopFeat(verbose = FALSE))
  list(
    queryData = function() perf,
    topFeatures = function() top,
    available_species = function() c("Shigella flexneri" = "Sfl"),
    filtered_top_features = function() {
      top[is.na(top$strat_label) & !top$cross_test, ]
    },
    loadDrugClassMapRec = function() tibble::tibble()
  )
}

test_that("setupServerCore loads the demo results and lists their species", {
  shiny::testServer(
    function(id) {
      shiny::moduleServer(id, function(input, output, session) {
        core <- setupServerCore(NULL, shiny::reactive(NULL))
      })
    },
    {
      session$flushReact()
      expect_gt(nrow(core$queryData()), 0)
      expect_true(all(c("Sfl", "Sso") %in% core$available_species()))
    }
  )
})

test_that("resultsSourceServer returns no selection in demo mode", {
  shiny::testServer(resultsSourceServer, args = list(results_root = NULL), {
    expect_null(session$getReturned()())
  })
})

test_that("resultsSourceServer offers the species folders under a results root", {
  root <- system.file("extdata", package = "amRviz")
  shiny::testServer(resultsSourceServer, args = list(results_root = root), {
    expect_match(
      output$results_species_selector_ui$html, "Shigella_flexneri"
    )
  })
})

test_that("modelPerfServer renders the performance plots for a species", {
  shiny::testServer(
    modelPerfServer,
    args = list(core = demo_core()),
    {
      session$setInputs(
        bug_perf_id = "Sfl", model_scale = "genes", data_type = "binary",
        model_metrics = "mcc", drug_class_perf_id = "all",
        drug_perf_id = "GEN"
      )
      expect_true(nzchar(output$model_perfomance_plot))
      expect_true(nzchar(output$mcc_heatmap))
    }
  )
})

test_that("featureImportanceServer renders the across-bug plot and table", {
  shiny::testServer(
    featureImportanceServer,
    args = list(core = demo_core()),
    {
      session$setInputs(
        across_bug_id = "drug", bug_search_amr_across_bug = "Sfl",
        bug_drug_comp_model_scale = "genes", feature_data_type = "binary",
        top_n_features = 10, feature_importance_tabset = "across_bug",
        amr_drug_across_bug = "GEN"
      )
      expect_identical(bug_norm_input(), "Sfl")
      expect_true(nzchar(output$across_bug_feature_importance_plot))
      expect_gt(nrow(enriched_across_bug()), 0)
    }
  )
})

test_that("crossModelComparisonServer renders the holdout plots", {
  shiny::testServer(
    crossModelComparisonServer,
    args = list(core = demo_core()),
    {
      session$setInputs(bug_cross_model_comparison_id = "Sfl")
      expect_true(nzchar(output$cross_model_ridge_country))
    }
  )
})

test_that("networkServer renders the drug-feature network", {
  shiny::testServer(
    networkServer,
    args = list(core = demo_core()),
    {
      session$setInputs(network_bug_id = "Sfl", network_top_n = 5)
      expect_true(nzchar(output$drug_feature_network))
    }
  )
})

test_that("queryDataServer shows the chosen performance columns", {
  shiny::testServer(
    queryDataServer,
    args = list(core = demo_core()),
    {
      session$setInputs(query_data_columns = c("species", "mcc"))
      expect_true(nzchar(output$queryDataTable))
    }
  )
})
