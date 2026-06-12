## Tests for plots ported from amRml PR #8:
## makeTopFeatsVIPlot, makeMDRPerformancePlot, and the loadMDRResults loader.

# ── makeTopFeatsVIPlot ───────────────────────────────────────────────────────

make_topfeat <- function() {
  tibble::tibble(
    species = "Sfl",
    drug_or_class = "AMP",
    feature_type = "genes",
    feature_subtype = "binary",
    strat_label = NA_character_,
    Variable = c("geneA", "geneB", "geneC", "geneD"),
    Importance = c(0.9, 0.7, 0.5, 0.3),
    Sign = c("POS", "NEG", "POS", "NEG")
  )
}

test_that("makeTopFeatsVIPlot returns a plotly object for a valid selection", {
  result <- makeTopFeatsVIPlot(
    make_topfeat(),
    bug = "Sfl", amr_drug = "AMP",
    model_scale = "genes", data_type_ = "binary",
    top_n_features = 3
  )
  expect_s3_class(result, "plotly")
})

test_that("makeTopFeatsVIPlot returns NULL when nothing matches the filter", {
  expect_null(makeTopFeatsVIPlot(
    make_topfeat(),
    bug = "Zzz", amr_drug = "nope",
    model_scale = "genes", data_type_ = "binary"
  ))
})

test_that("makeTopFeatsVIPlot handles a group with all-NA Sign without error", {
  df <- make_topfeat()
  df$Sign <- NA_character_
  result <- makeTopFeatsVIPlot(
    df,
    bug = "Sfl", amr_drug = "AMP",
    model_scale = "genes", data_type_ = "binary"
  )
  expect_s3_class(result, "plotly")
})

test_that("makeTopFeatsVIPlot returns NULL for empty or invalid input", {
  expect_null(makeTopFeatsVIPlot(
    NULL, "Sfl", "AMP", "genes", "binary"
  ))
  expect_null(makeTopFeatsVIPlot(
    tibble::tibble(a = 1), "Sfl", "AMP", "genes", "binary"
  ))
})

# ── makeMDRPerformancePlot ───────────────────────────────────────────────────

make_mdr <- function() {
  tibble::tibble(
    feature_type = c("genes", "genes", "domains", "domains"),
    feature_subtype = c("binary", "counts", "binary", "counts"),
    nmcc = c(0.8, 0.6, 0.7, 0.5)
  )
}

test_that("makeMDRPerformancePlot returns a plotly object", {
  expect_s3_class(makeMDRPerformancePlot(make_mdr()), "plotly")
})

test_that("makeMDRPerformancePlot returns NULL for empty or invalid input", {
  expect_null(makeMDRPerformancePlot(NULL))
  expect_null(makeMDRPerformancePlot(data.frame()))
  expect_null(makeMDRPerformancePlot(tibble::tibble(a = 1)))
})

test_that("makeMDRPerformancePlot returns NULL when all nmcc are NA", {
  df <- make_mdr()
  df$nmcc <- NA_real_
  expect_null(makeMDRPerformancePlot(df))
})

# ── loadMDRResults + end-to-end with bundled demo data ───────────────────────

test_that("loadMDRResults loads bundled MDR demo data", {
  mdr <- loadMDRResults(verbose = FALSE)
  skip_if(nrow(mdr) == 0, "No bundled MDR demo data")
  expect_true("nmcc" %in% names(mdr))
  expect_true("feature_type" %in% names(mdr))
  expect_s3_class(makeMDRPerformancePlot(mdr), "plotly")
})

test_that("makeTopFeatsVIPlot renders from bundled demo top features", {
  tf <- loadTopFeat(verbose = FALSE)
  skip_if(nrow(tf) == 0, "No bundled top-feature demo data")

  bug <- normalize_species(unique(tf$species))[1]
  drug <- tf$drug_or_class[!is.na(tf$drug_or_class)][1]
  result <- makeTopFeatsVIPlot(
    tf,
    bug = bug, amr_drug = drug,
    model_scale = "genes", data_type_ = "binary",
    top_n_features = 10
  )
  # Either a plotly object (matches found) or NULL (no genes/binary rows for
  # this drug) — both are valid; assert it does not error.
  expect_true(is.null(result) || inherits(result, "plotly"))
})
