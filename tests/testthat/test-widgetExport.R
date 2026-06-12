## Tests for the static (PDF) export of interactive widgets:
## .writeWidgetStatic and the .exportBtn UI helper.

# ── .exportBtn ───────────────────────────────────────────────────────────────

test_that(".exportBtn builds a download button with the paired id", {
  btn <- .exportBtn("model_perfomance_plot")
  html <- as.character(btn)
  expect_match(html, "model_perfomance_plot_static_download")
  expect_match(html, "Export \\(PDF\\)")
  expect_match(html, "btn-export")
})

# ── .writeWidgetStatic ───────────────────────────────────────────────────────

test_that(".writeWidgetStatic writes a placeholder PDF for a NULL widget", {
  tmp <- tempfile(fileext = ".pdf")
  .writeWidgetStatic(tmp, NULL)
  expect_gt(file.info(tmp)$size, 0)
})

test_that(".writeWidgetStatic snapshots a plotly widget to a non-empty PDF", {
  skip_if_not_installed("webshot2")
  skip_if_not_installed("chromote")
  skip_if(
    !nzchar(tryCatch(chromote::find_chrome(), error = function(e) "")),
    "No headless Chrome available"
  )

  p <- plotly::plot_ly(x = c("a", "b"), y = c(1, 2), type = "bar")
  tmp <- tempfile(fileext = ".pdf")
  .writeWidgetStatic(tmp, p)
  expect_gt(file.info(tmp)$size, 0)
})
