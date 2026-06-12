# Heatmap / widget export helpers (PDF + static PNG).


# .plotlyHeatmapMatrix: pull the z matrix (with x/y as column/row names) out of a
# plotly heatmap object built by any of the dashboard's plot_ly(type="heatmap")
# renderers. Returns a labelled numeric matrix, or NULL if the object is not a
# heatmap. This lets a single ComplexHeatmap exporter serve every heatmap panel
# without re-deriving each one's matrix.
.plotlyHeatmapMatrix <- function(p) {
  if (is.null(p) || !inherits(p, c("plotly", "htmlwidget"))) {
    return(NULL)
  }
  b <- tryCatch(plotly::plotly_build(p), error = function(e) NULL)
  if (is.null(b) || is.null(b$x$data)) {
    return(NULL)
  }
  traces <- b$x$data
  idx <- which(vapply(
    traces, function(t) identical(t$type, "heatmap"), logical(1)
  ))
  if (!length(idx)) {
    return(NULL)
  }
  d <- traces[[idx[1]]]
  if (is.null(d$z)) {
    return(NULL)
  }
  mat <- suppressWarnings(as.matrix(d$z))
  if (!is.numeric(mat) || !length(mat)) {
    return(NULL)
  }
  rows <- unlist(d$y)
  cols <- unlist(d$x)
  if (length(rows) == nrow(mat)) rownames(mat) <- rows
  if (length(cols) == ncol(mat)) colnames(mat) <- cols
  mat
}

# makeHeatmapExport: convert any dashboard heatmap (a plotly object) into a static
# ComplexHeatmap::Heatmap for publication-quality PDF export. Colours use an RdBu
# diverging ramp mapped to the data range. Returns a Heatmap object, or NULL when
# the input is not a heatmap or has no finite values.
makeHeatmapExport <- function(p, name = "value",
                              column_title = NULL, row_title = NULL) {
  mat <- .plotlyHeatmapMatrix(p)
  if (is.null(mat)) {
    return(NULL)
  }
  finite_vals <- mat[is.finite(mat)]
  if (!length(finite_vals)) {
    return(NULL)
  }
  rng <- range(finite_vals)
  if (rng[1] == rng[2]) {
    rng <- rng + c(-0.5, 0.5)
  }
  heat_colors <- grDevices::colorRampPalette(
    RColorBrewer::brewer.pal(11, "RdBu")
  )(100)
  col_fun <- circlize::colorRamp2(
    seq(rng[1], rng[2], length.out = length(heat_colors)),
    heat_colors
  )

  tryCatch(
    ComplexHeatmap::Heatmap(
      mat,
      name = name,
      col = col_fun,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      column_title = column_title,
      row_title = row_title,
      na_col = "grey95",
      rect_gp = grid::gpar(col = NA),
      row_names_gp = grid::gpar(fontsize = 8),
      column_names_gp = grid::gpar(fontsize = 10),
      show_heatmap_legend = TRUE
    ),
    error = function(e) NULL
  )
}

# .writeHeatmapPdf: render a ComplexHeatmap object to a PDF file for a Shiny
# downloadHandler, or a placeholder message when there is nothing to draw.
.writeHeatmapPdf <- function(file, hm, width = 9, height = 7,
                             empty_msg = "No data available for this selection.") {
  grDevices::pdf(file, width = width, height = height)
  on.exit(grDevices::dev.off())
  drawn <- if (!is.null(hm)) {
    isTRUE(tryCatch(
      {
        ComplexHeatmap::draw(hm)
        TRUE
      },
      error = function(e) FALSE
    ))
  } else {
    FALSE
  }
  if (!drawn) {
    grid::grid.newpage()
    grid::grid.text(empty_msg, gp = grid::gpar(fontsize = 12, col = "grey40"))
  }
}

# .exportBtn: small right-aligned "Export (PDF)" button for an interactive
# panel, paired with the server's "<id>_static_download" downloadHandler.
.exportBtn <- function(id) {
  shiny::div(
    style = "text-align: right; padding: 2px 6px;",
    shiny::downloadButton(
      paste0(id, "_static_download"), "Export (PDF)",
      class = "btn btn-export"
    )
  )
}

# .writeWidgetStatic: render any interactive htmlwidget (plotly, networkD3
# force/sankey) to a static PDF for a download handler, via a headless-browser
# snapshot (webshot2 + chromote). Used for the dashboard panels that are not
# plain heatmaps (which use the vector ComplexHeatmap path instead). Falls back
# to a placeholder PDF when the widget is NULL or the browser tooling/Chrome is
# unavailable, so a download never hard-errors.
.writeWidgetStatic <- function(file, widget, vwidth = 1000, vheight = 700,
                               empty_msg = "No data available for this selection.") {
  ok <- !is.null(widget) &&
    requireNamespace("htmlwidgets", quietly = TRUE) &&
    requireNamespace("webshot2", quietly = TRUE)
  if (ok) {
    ok <- tryCatch(
      {
        html <- tempfile(fileext = ".html")
        out <- tempfile(fileext = ".pdf")
        on.exit(unlink(c(html, out)), add = TRUE)
        htmlwidgets::saveWidget(widget, html, selfcontained = TRUE)
        webshot2::webshot(
          html, out,
          delay = 0.8, vwidth = vwidth, vheight = vheight
        )
        file.copy(out, file, overwrite = TRUE)
      },
      error = function(e) FALSE
    )
  }
  if (!isTRUE(ok)) {
    grDevices::pdf(file)
    on.exit(grDevices::dev.off(), add = TRUE)
    grid::grid.newpage()
    grid::grid.text(empty_msg, gp = grid::gpar(fontsize = 12, col = "grey40"))
  }
}
