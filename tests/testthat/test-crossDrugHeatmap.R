## Tests for the cross-drug generalization heatmap functions
## (.buildCrossDrugMatrix, makeCrossDrugHeatmap, makeCrossDrugHeatmapPlotly)

# Minimal synthetic performance tibble with the columns the cross-drug
# functions require: cross rows (cross_test == TRUE, ref_drug/test_drug set)
# plus baseline diagonal rows (cross_test == FALSE, strat_label NA).
make_cross_perf <- function() {
  # Off-diagonal cross-testing only (a drug is not cross-tested against itself);
  # the matrix diagonal is supplied by the baseline rows below.
  cross <- tibble::tibble(
    cross_test = TRUE,
    drug_label = "drug",
    strat_label = NA_character_,
    ref_drug = c("AMP", "TET"),
    test_drug = c("TET", "AMP"),
    drug_or_class = c("AMP", "TET"),
    nmcc = c(0.40, 0.30),
    species = "cross"
  )
  diag <- tibble::tibble(
    cross_test = FALSE,
    drug_label = "drug",
    strat_label = NA_character_,
    ref_drug = NA_character_,
    test_drug = NA_character_,
    drug_or_class = c("AMP", "TET"),
    nmcc = c(0.95, 0.85),
    species = "Escherichia_coli"
  )
  dplyr::bind_rows(cross, diag)
}

make_cross_meta <- function() {
  tibble::tibble(
    drug_abbr = c("AMP", "TET"),
    class_abbr = c("Penicillins", "Tetracyclines")
  )
}

# ── .buildCrossDrugMatrix ────────────────────────────────────────────────────

test_that(".buildCrossDrugMatrix returns an ordered square matrix", {
  prep <- .buildCrossDrugMatrix(make_cross_perf(), meta = make_cross_meta())

  expect_type(prep, "list")
  expect_true(is.matrix(prep$mat))
  expect_equal(dim(prep$mat), c(2L, 2L))
  expect_setequal(rownames(prep$mat), c("AMP", "TET"))
  expect_setequal(colnames(prep$mat), c("AMP", "TET"))
  # Diagonal comes from the baseline self-evaluation rows.
  expect_equal(prep$mat["AMP", "AMP"], 0.95)
  expect_equal(prep$mat["TET", "TET"], 0.85)
  # Class lookup is populated from the metadata.
  expect_equal(unname(prep$class_of["AMP"]), "Penicillins")
})

test_that(".buildCrossDrugMatrix works without metadata", {
  prep <- .buildCrossDrugMatrix(make_cross_perf(), meta = NULL)

  expect_true(is.matrix(prep$mat))
  expect_true(all(is.na(prep$class_of)))
})

test_that(".buildCrossDrugMatrix returns NULL for empty or invalid input", {
  expect_null(.buildCrossDrugMatrix(NULL))
  expect_null(.buildCrossDrugMatrix(data.frame()))
  # Missing the required cross-drug columns.
  expect_null(.buildCrossDrugMatrix(tibble::tibble(a = 1, b = 2)))
})

test_that(".buildCrossDrugMatrix returns NULL with no cross-drug rows", {
  baseline_only <- make_cross_perf() |>
    dplyr::filter(.data$cross_test %in% FALSE)
  expect_null(.buildCrossDrugMatrix(baseline_only))
})

# ── makeCrossDrugHeatmap (ComplexHeatmap, static export) ─────────────────────

test_that("makeCrossDrugHeatmap returns a ComplexHeatmap Heatmap object", {
  hm <- makeCrossDrugHeatmap(make_cross_perf(), meta = make_cross_meta())
  expect_s4_class(hm, "Heatmap")
})

test_that("makeCrossDrugHeatmap returns NULL for empty input", {
  expect_null(makeCrossDrugHeatmap(NULL))
  expect_null(makeCrossDrugHeatmap(data.frame()))
})

test_that("makeCrossDrugHeatmap renders to a non-empty PDF", {
  hm <- makeCrossDrugHeatmap(make_cross_perf(), meta = make_cross_meta())
  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp)
  ComplexHeatmap::draw(hm)
  grDevices::dev.off()
  expect_gt(file.info(tmp)$size, 0)
})

# ── makeCrossDrugHeatmapPlotly (interactive dashboard panel) ─────────────────

test_that("makeCrossDrugHeatmapPlotly returns a plotly object", {
  result <- makeCrossDrugHeatmapPlotly(
    make_cross_perf(),
    meta = make_cross_meta()
  )
  expect_s3_class(result, "plotly")
})

test_that("makeCrossDrugHeatmapPlotly returns NULL for empty input", {
  expect_null(makeCrossDrugHeatmapPlotly(NULL))
  expect_null(makeCrossDrugHeatmapPlotly(data.frame()))
})

# ── End-to-end with bundled demo data ────────────────────────────────────────

test_that("makeCrossDrugHeatmap renders from bundled Shigella demo data", {
  fp <- system.file(
    "extdata", "Shigella_flexneri", "Sfl_cross_ML_perf.parquet",
    package = "amRviz"
  )
  skip_if(!nzchar(fp) || !file.exists(fp), "No bundled cross-drug demo data")

  perf <- arrow::read_parquet(fp)
  mfp <- system.file(
    "extdata", "Shigella_flexneri", "Sfl_metadata.parquet",
    package = "amRviz"
  )
  meta <- if (nzchar(mfp) && file.exists(mfp)) {
    arrow::read_parquet(mfp)
  } else {
    NULL
  }

  hm <- makeCrossDrugHeatmap(perf, meta = meta)
  expect_s4_class(hm, "Heatmap")

  ply <- makeCrossDrugHeatmapPlotly(perf, meta = meta)
  expect_s3_class(ply, "plotly")
})
