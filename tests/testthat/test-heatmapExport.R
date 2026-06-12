## Tests for the generic heatmap PDF export helpers:
## .plotlyHeatmapMatrix, makeHeatmapExport, .writeHeatmapPdf

make_heatmap_plotly <- function() {
  m <- matrix(
    c(0.1, 0.9, 0.5, 0.3),
    nrow = 2,
    dimnames = list(c("r1", "r2"), c("c1", "c2"))
  )
  plotly::plot_ly(x = colnames(m), y = rownames(m), z = m, type = "heatmap")
}

make_bar_plotly <- function() {
  plotly::plot_ly(x = c("a", "b"), y = c(1, 2), type = "bar")
}

# ── .plotlyHeatmapMatrix ─────────────────────────────────────────────────────

test_that(".plotlyHeatmapMatrix extracts a labelled matrix from a heatmap", {
  mat <- .plotlyHeatmapMatrix(make_heatmap_plotly())
  expect_true(is.matrix(mat))
  expect_equal(dim(mat), c(2L, 2L))
  expect_setequal(rownames(mat), c("r1", "r2"))
  expect_setequal(colnames(mat), c("c1", "c2"))
})

test_that(".plotlyHeatmapMatrix returns NULL for non-heatmap / NULL input", {
  expect_null(.plotlyHeatmapMatrix(make_bar_plotly()))
  expect_null(.plotlyHeatmapMatrix(NULL))
})

# ── makeHeatmapExport ────────────────────────────────────────────────────────

test_that("makeHeatmapExport returns a ComplexHeatmap from a heatmap plotly", {
  hm <- makeHeatmapExport(make_heatmap_plotly(), name = "v")
  expect_s4_class(hm, "Heatmap")
})

test_that("makeHeatmapExport returns NULL for non-heatmap / NULL input", {
  expect_null(makeHeatmapExport(make_bar_plotly()))
  expect_null(makeHeatmapExport(NULL))
})

# ── .writeHeatmapPdf ─────────────────────────────────────────────────────────

test_that(".writeHeatmapPdf renders a non-empty PDF for a heatmap", {
  hm <- makeHeatmapExport(make_heatmap_plotly())
  tmp <- tempfile(fileext = ".pdf")
  .writeHeatmapPdf(tmp, hm)
  expect_gt(file.info(tmp)$size, 0)
})

test_that(".writeHeatmapPdf renders a placeholder PDF when hm is NULL", {
  tmp <- tempfile(fileext = ".pdf")
  .writeHeatmapPdf(tmp, NULL)
  expect_gt(file.info(tmp)$size, 0)
})

# ── End-to-end: dashboard heatmaps export to ComplexHeatmap ──────────────────

test_that("dashboard heatmaps convert to ComplexHeatmap exports", {
  tf <- loadTopFeat(verbose = FALSE)
  perf <- loadMLResults(verbose = FALSE)
  skip_if(nrow(tf) == 0 || nrow(perf) == 0, "No demo data")

  drugs <- unique(tf$drug_or_class[!is.na(tf$drug_or_class)])

  fi <- makeFeatureImportancePlot(
    tf, "Sfl", drugs[seq_len(min(3, length(drugs)))],
    "genes", "binary", 10, "across_drug"
  )
  expect_s4_class(makeHeatmapExport(fi, name = "Importance"), "Heatmap")

  perf_hm <- makeCrossModelPerformancePlot(perf, "Sfl", drugs[1], "country")
  skip_if(is.null(perf_hm), "No cross-model performance demo data")
  expect_s4_class(makeHeatmapExport(perf_hm, name = "nMCC"), "Heatmap")
})
