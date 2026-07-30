# Headless batch export of every dashboard visualization to static image
# files. Lets a user render the full set of amRviz figures from the packaged
# demo data (or their own amRml results) without ever launching the Shiny app.


#' Snapshot one htmlwidget to static image file(s)
#'
#' Saves the widget to a temp HTML bundle, renders it in a headless Chrome,
#' then writes each requested format. Raster formats (png/jpg/pdf) all derive
#' from a single webshot PNG so they share crop and dimensions; svg is
#' extracted from the rendered DOM via chromote.
#'
#' @param widget An htmlwidget (plotly / networkD3), or NULL.
#' @param path_base Output path without extension; each format is written as
#'   `path_base.<ext>`.
#' @param formats Character vector of extensions among png/pdf/jpg/jpeg/svg.
#' @param width,height Snapshot viewport size in pixels.
#' @param scale Device-pixel multiplier for the underlying PNG render.
#' @param trim When TRUE, `magick::image_trim()` crops the PNG before other
#'   raster formats are derived from it. Meant for figures with known excess
#'   whitespace (currently the drug-feature network only).
#' @param delay Seconds to let the widget's JavaScript render before capture.
#' @param verbose Whether to message on per-format failures.
#' @return Character vector of files actually written.
#' @keywords internal
#' @noRd
.exportWidgetFile <- function(widget, path_base, formats,
                              width = 1200, height = 800, scale = 2,
                              trim = FALSE,
                              delay = 1.5, verbose = TRUE) {
  if (is.null(widget)) return(character(0))
  written <- character(0)
  raster <- intersect(formats, c("png", "pdf", "jpg", "jpeg"))
  want_svg <- "svg" %in% formats
  if (!length(raster) && !want_svg) return(written)

  tmpdir <- tempfile("amrviz_widget_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)
  html <- file.path(tmpdir, "widget.html")
  # selfcontained = TRUE inlines JS/CSS so the extracted SVG stands alone.
  saved <- tryCatch(
    {
      htmlwidgets::saveWidget(widget, html, selfcontained = TRUE)
      TRUE
    },
    error = function(e) {
      if (verbose) message("    saveWidget failed: ", conditionMessage(e))
      FALSE
    }
  )
  if (!saved) return(written)

  if (length(raster)) {
    tmp_png <- file.path(tmpdir, "render.png")
    rendered <- tryCatch(
      {
        webshot2::webshot(
          html, tmp_png,
          vwidth = width, vheight = height,
          zoom = scale, delay = delay, quiet = TRUE
        )
        file.exists(tmp_png) && file.info(tmp_png)$size > 0
      },
      error = function(e) {
        if (verbose) message("    render failed: ", conditionMessage(e))
        FALSE
      }
    )
    if (rendered) {
      if (trim) .trim_raster(tmp_png, verbose = verbose)
      for (ext in raster) {
        out <- paste0(path_base, ".", ext)
        ok <- if (ext == "png") {
          file.copy(tmp_png, out, overwrite = TRUE)
        } else {
          .convert_raster(tmp_png, out, ext, verbose = verbose)
        }
        if (isTRUE(ok)) written <- c(written, out)
      }
    }
  }

  if (want_svg) {
    out <- paste0(path_base, ".svg")
    ok <- tryCatch(
      .save_widget_svg(html, out, width, height, delay),
      error = function(e) {
        if (verbose) message("    svg failed: ", conditionMessage(e))
        FALSE
      }
    )
    if (isTRUE(ok)) written <- c(written, out)
  }

  written
}


#' Re-encode a PNG into another raster/PDF format via magick
#'
#' PDF is written as a rasterised image wrapped in a PDF page (not vector).
#'
#' @param src PNG source path.
#' @param out Destination path.
#' @param ext Target extension ("jpg", "jpeg", "pdf").
#' @param verbose Whether to message on failure.
#' @return TRUE on success, FALSE otherwise.
#' @keywords internal
#' @noRd
.convert_raster <- function(src, out, ext, verbose = TRUE) {
  fmt <- if (ext %in% c("jpg", "jpeg")) "jpeg" else ext
  tryCatch(
    {
      magick::image_write(magick::image_read(src), out, format = fmt)
      file.exists(out) && file.info(out)$size > 0
    },
    error = function(e) {
      if (verbose) message("    ", ext, " failed: ", conditionMessage(e))
      FALSE
    }
  )
}


#' Extract an htmlwidget's rendered SVG element into a standalone .svg file
#'
#' Plotly widgets go through `Plotly.toImage` so title/axis/legend overlays
#' end up inlined as SVG text; other widgets use the raw SVG outerHTML.
#'
#' @param html Path to the saved widget HTML.
#' @param out Path to write the SVG to.
#' @param width,height Viewport size the widget renders into.
#' @param delay Seconds to wait after navigation before extracting.
#' @return TRUE on success, FALSE if the widget has no SVG or writes fail.
#' @keywords internal
#' @noRd
.save_widget_svg <- function(html, out, width, height, delay) {
  b <- chromote::ChromoteSession$new()
  on.exit(b$close(), add = TRUE)
  b$Emulation$setDeviceMetricsOverride(
    width = as.integer(width), height = as.integer(height),
    deviceScaleFactor = 1, mobile = FALSE
  )
  b$Page$navigate(paste0("file://", html))
  Sys.sleep(delay)
  js <- sprintf(
    "(async function() {
       var gd = document.querySelector('.js-plotly-plot');
       if (gd && window.Plotly) {
         var url = await Plotly.toImage(gd, {format: 'svg', width: %d, height: %d});
         return decodeURIComponent(url.replace(/^data:image\\/svg\\+xml,/, ''));
       }
       var s = document.querySelector('svg');
       if (!s) return null;
       if (!s.getAttribute('xmlns')) s.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
       return s.outerHTML;
     })()",
    as.integer(width), as.integer(height)
  )
  res <- b$Runtime$evaluate(js, returnByValue = TRUE, awaitPromise = TRUE)
  svg_str <- res$result$value
  if (is.null(svg_str) || !nzchar(svg_str)) return(FALSE)
  writeLines(svg_str, out)
  file.exists(out) && file.info(out)$size > 0
}


#' Baseline (non-stratified) drug labels available for a species
#' @keywords internal
#' @noRd
.export_ml_drugs <- function(perf_data, code) {
  if (is.null(perf_data) || !nrow(perf_data)) {
    return(character(0))
  }
  perf_data |>
    dplyr::filter(
      normalize_species(.data$species) %in% normalize_species(code)
    ) |>
    dplyr::filter(is.na(.data$strat_label) | !nzchar(.data$strat_label)) |>
    dplyr::filter(.data$drug_label == "drug") |>
    dplyr::pull(.data$drug_or_class) |>
    unique() |>
    sort()
}


#' Baseline top-feature drug labels for a species (feature-importance defaults)
#' @keywords internal
#' @noRd
.export_tf_drugs <- function(top_features, code,
                             scale = "genes",
                             subtype = c("counts", "binary")) {
  if (is.null(top_features) || !nrow(top_features)) {
    return(character(0))
  }
  top_features |>
    dplyr::filter(
      normalize_species(.data$species) %in% normalize_species(code)
    ) |>
    dplyr::filter(.data$feature_type %in% scale) |>
    dplyr::filter(.data$feature_subtype %in% subtype) |>
    dplyr::filter(is.na(.data$strat_label) | !nzchar(.data$strat_label)) |>
    dplyr::filter(.data$drug_label == "drug") |>
    dplyr::pull(.data$drug_or_class) |>
    unique() |>
    sort()
}


#' Crop uniform-white margins from a raster file in place
#'
#' No-op when `magick` isn't installed or the trim errors.
#'
#' @param path Raster file to overwrite in place.
#' @param verbose Whether to message on failure.
#' @keywords internal
#' @noRd
.trim_raster <- function(path, verbose = TRUE) {
  tryCatch(
    magick::image_write(magick::image_trim(magick::image_read(path)), path),
    error = function(e) {
      if (verbose) message("    trim failed: ", conditionMessage(e))
    }
  )
  invisible()
}


#' Fit a networkD3 forceNetwork to its container after simulation settles
#'
#' Waits for the force simulation, then makes the SVG fill its container and
#' rewrites the viewBox to the graph's bbox. Export-only.
#'
#' @param widget A `networkD3::forceNetwork` htmlwidget.
#' @return The widget with the onRender hook attached.
#' @keywords internal
#' @noRd
.fit_network_to_content <- function(widget) {
  htmlwidgets::onRender(widget, "
    function(el, x) {
      setTimeout(function() {
        try {
          var svg = el.querySelector('svg');
          svg.setAttribute('width', '100%');
          svg.setAttribute('height', '100%');
          var bb = svg.getBBox();
          var p = 40;
          svg.setAttribute('viewBox',
            (bb.x - p) + ' ' + (bb.y - p) + ' ' +
            (bb.width + 2*p) + ' ' + (bb.height + 2*p));
        } catch (e) {}
      }, 1300);
    }
  ")
}


#' Assemble the list of figures to export
#'
#' Returns a list of specs; each is `list(group, name, build)` where `build()`
#' lazily produces the widget. Grouping mirrors the dashboard tabs. Global
#' overviews live under `_overview` / `_across_species`; everything else is
#' filed under its species folder. Selections replicate the dashboard defaults.
#'
#' @keywords internal
#' @noRd
.exportPlanSpecs <- function(perf_data, top_features,
                             ml_species, meta_species,
                             results_root, amrdata_root,
                             top_n_features = 10, network_top_n = 5) {
  specs <- list()
  add <- function(group, name, build,
                  width = NULL, height = NULL, trim = FALSE) {
    specs[[length(specs) + 1]] <<- list(
      group = group, name = name, build = build,
      width = width, height = height, trim = trim
    )
  }

  model_scale_all <- c("genes", "proteins", "domains", "cogs", "args")
  fi_scale <- "genes"
  fi_subtype <- c("counts", "binary")

  # Defensive: drop species with no code. Blank/NA codes share a species_label
  # folder with a real species, so their per-species specs would collide with
  # (and overwrite) the real ones. Guaranteeing unique output paths here keeps
  # the plan robust regardless of how the caller derived `ml_species`.
  valid <- !is.na(ml_species$code) & nzchar(ml_species$code)
  ml_species <- list(
    code = ml_species$code[valid], label = ml_species$label[valid]
  )

  # ---- Global performance overviews (all species at once) ----
  if (nrow(perf_data)) {
    add("_overview", "mcc_strip_by_species_scale", function() {
      makeMCCStripPlot(perf_data, selected_drug_class = "all")
    })
    add("_overview", "mcc_heatmap_overview", function() {
      makeMCCHeatmap(perf_data, selected_drug_class = "all")
    })
  }

  # ---- Cross-species feature importance (across_bug view) ----
  if (nrow(top_features) && length(ml_species$code)) {
    pooled_drugs <- .export_tf_drugs(
      top_features, ml_species$code, fi_scale, fi_subtype
    )
    if (length(pooled_drugs)) {
      drug <- if ("GEN" %in% pooled_drugs) "GEN" else pooled_drugs[1]
      add("_across_species", paste0("feature_importance_", drug), function() {
        makeFeatureImportancePlot(
          top_features, ml_species$code, drug,
          fi_scale, fi_subtype, top_n_features, "across_bug",
          amrdata_root = amrdata_root, results_root = results_root
        )
      })
    }
  }

  # ---- Per ML species: performance, feature importance, holdouts, network ----
  for (i in seq_along(ml_species$code)) {
    code <- ml_species$code[i]
    folder <- ml_species$label[i]

    # Model performance box/point plot.
    drugs <- .export_ml_drugs(perf_data, code)
    drug_sel <- if ("GEN" %in% drugs) {
      "GEN"
    } else if (length(drugs)) {
      drugs[1]
    } else {
      NULL
    }
    local({
      code <- code
      drug_sel <- drug_sel
      add(folder, "model_performance", function() {
        makeModelPerformancePlot(
          perf_data, code, model_scale_all, "binary", "mcc",
          "all", drug_sel
        )
      })
    })

    # Feature importance across drugs for this species.
    tf_drugs <- .export_tf_drugs(top_features, code, fi_scale, fi_subtype)
    fi_drugs <- intersect(c("OXA", "PEN", "MET"), tf_drugs)
    if (!length(fi_drugs)) {
      fi_drugs <- utils::head(tf_drugs, min(3, length(tf_drugs)))
    }
    if (length(fi_drugs)) {
      local({
        code <- code
        fi_drugs <- fi_drugs
        add(folder, "feature_importance_across_drugs", function() {
          makeFeatureImportancePlot(
            top_features, code, fi_drugs,
            fi_scale, fi_subtype, top_n_features, "across_drug",
            amrdata_root = amrdata_root, results_root = results_root
          )
        })
        # COG category barplot from annotation-enriched features.
        add(folder, "cog_categories", function() {
          tf <- top_features |>
            dplyr::filter(
              is.na(.data$strat_label) | !nzchar(.data$strat_label)
            ) |>
            dplyr::filter(!isTRUE(.data$cross_test)) |>
            dplyr::filter(
              normalize_species(.data$species) %in% normalize_species(code)
            ) |>
            dplyr::filter(.data$drug_or_class %in% fi_drugs)
          if (!nrow(tf)) {
            return(makeCogBarChart(NULL))
          }
          enriched <- dplyr::bind_rows(lapply(unique(tf$species), function(sp) {
            enrich_with_annotations(
              tf[tf$species == sp, ],
              species_code = sp, results_root = results_root
            )
          }))
          makeCogBarChart(enriched)
        })
      })
    }

    # Cross-model holdouts (country + time strata).
    holdout_drugs <- getHoldoutsDrugChoices(perf_data, code)
    holdout_drug <- if (length(holdout_drugs)) holdout_drugs[1] else NULL
    local({
      code <- code
      holdout_drug <- holdout_drug
      add(folder, "holdout_ridge_country", function() {
        makeCrossModelRidgePlot(perf_data, code, "country")
      })
      add(folder, "holdout_ridge_time", function() {
        makeCrossModelRidgePlot(perf_data, code, "time")
      })
      if (!is.null(holdout_drug)) {
        add(
          folder, paste0("holdout_performance_country_", holdout_drug),
          function() {
            makeCrossModelPerformancePlot(
              perf_data, code, holdout_drug, "country"
            )
          }
        )
        add(
          folder, paste0("holdout_performance_time_", holdout_drug),
          function() {
            makeCrossModelPerformancePlot(
              perf_data, code, holdout_drug, "time"
            )
          }
        )
        add(
          folder, paste0("holdout_feature_importance_country_", holdout_drug),
          function() {
            makeCrossModelFeatureImportancePlot(
              top_features, code, holdout_drug, "country", top_n_features
            )
          }
        )
      }
    })

    # Drug-feature network. Square viewport since the graph settles roughly
    # square; .fit_network_to_content() then fits the SVG to the actual bbox.
    local({
      code <- code
      add(folder, "drug_feature_network", function() {
        makeDrugFeatureNetwork(
          top_features, code,
          top_n = network_top_n,
          include_clusters = FALSE, include_cogs = FALSE,
          results_root = results_root
        ) |> .fit_network_to_content()
      }, width = 1600, height = 1600, trim = TRUE)
    })
  }

  # ---- Per metadata species: distributions + sankey ----
  for (sp in meta_species) {
    folder <- sp
    local({
      sp <- sp
      folder <- folder
      meta_raw <- function() .read_metadata_for_bug(sp, results_root)

      add(folder, "metadata_data_availability", function() {
        makeDatAvailabilityPlot(meta_raw())
      })
      add(folder, "metadata_geographic", function() {
        data <- meta_raw() |>
          dplyr::filter(.data$genome.isolation_country != "") |>
          .add_evidence_column() |>
          dplyr::group_by(
            .data$genome.isolation_country, .data$genome_drug.antibiotic
          ) |>
          dplyr::summarize(count = dplyr::n(), .groups = "drop") |>
          dplyr::group_by(.data$genome.isolation_country) |>
          dplyr::summarise(count = sum(.data$count), .groups = "drop")
        makeGeoChloroPlot(data)
      })
      add(folder, "metadata_resistance_over_time", function() {
        data <- meta_raw() |>
          .add_evidence_column() |>
          dplyr::filter(!is.na(.data$genome.collection_year)) |>
          dplyr::group_by(
            .data$genome_drug.antibiotic,
            .data$genome_drug.resistant_phenotype,
            .data$genome.isolation_country,
            .data$genome.collection_year
          ) |>
          dplyr::summarize(n = dplyr::n(), .groups = "drop")
        makeTimeSeriesAMRPlot(data, "all")
      })
      add(folder, "metadata_hosts", function() {
        data <- meta_raw() |>
          dplyr::filter(.data$genome.host_common_name != "") |>
          .add_evidence_column()
        makeHostIsolatePlot(data)
      })
      add(folder, "metadata_isolation_sources", function() {
        data <- meta_raw() |>
          dplyr::filter(.data$genome.host_common_name != "") |>
          .add_evidence_column()
        makeIsolationSourcesPlot(data)
      })
      add(folder, "metadata_resistance_sankey", function() {
        meta <- meta_raw()
        classes <- if (nrow(meta) && "drug_class" %in% names(meta)) {
          meta |>
            dplyr::filter(!is.na(.data$drug_class)) |>
            dplyr::count(.data$drug_class, name = "n") |>
            dplyr::arrange(dplyr::desc(.data$n)) |>
            dplyr::slice_head(n = 3) |>
            dplyr::pull(.data$drug_class)
        } else {
          NULL
        }
        makeMetadataSankey(meta, drug_classes = classes)
      })
    })
  }

  specs
}


#' Export all amRviz dashboard visualizations to static image files
#'
#' Renders every amRviz figure (metadata, model performance, feature
#' importance, cross-model holdouts, drug-feature networks) to files, without
#' launching the Shiny app. One figure set per species, with global overviews
#' under `_overview` / `_across_species`. Selections mirror the dashboard's
#' defaults; `top_n_features` and `network_top_n` override the corresponding
#' sliders.
#'
#' Widgets are driven by a headless Chrome: `png` is a webshot2 snapshot;
#' `jpg` and `pdf` are re-encodes of that PNG via \pkg{magick} (identical
#' crop and dimensions); `svg` is extracted from the DOM via \pkg{chromote},
#' using `Plotly.toImage` for plotly so titles and labels survive.
#'
#' @param output_dir Directory to write figures into; created if needed. Files
#'   are organised as `output_dir/<species>/<panel>.<ext>`, with global
#'   overviews under `output_dir/_overview` and `output_dir/_across_species`.
#' @param formats Any of `"png"`, `"pdf"`, `"jpg"`, `"svg"`.
#' @param results_root Directory of amRml model outputs (per-species subdirs
#'   of `*_perf.parquet` / `*_top_features.parquet` / `metadata.parquet`).
#'   `NULL` uses the packaged demo data.
#' @param amrdata_root Directory of amRdata annotation parquets (for COG
#'   enrichment). `NULL` tries `~/amRdata/data`, else falls back to unenriched.
#' @param species Optional character vector restricting which species folders
#'   to export.
#' @param top_n_features,network_top_n Feature counts for the
#'   feature-importance panels and the drug-feature network.
#' @param width,height Viewport size in pixels.
#' @param scale Device-pixel multiplier for the PNG render; higher = sharper.
#' @param delay Seconds to wait for each widget's JavaScript to settle.
#' @param verbose Per-figure progress messages.
#'
#' @return Invisibly, a data frame with one row per figure: `group`, `name`,
#'   `written` (file count), `ok`.
#' @export
#' @examples
#' if (interactive()) {
#'   # Export the packaged demo figures as PNG + PDF into ./amRviz_exports
#'   exportAMRVisualizations()
#'
#'   # Your own results, PNG only, one species, more features per panel
#'   exportAMRVisualizations(
#'     output_dir = "figs",
#'     formats = "png",
#'     results_root = "~/my_amRml_results",
#'     species = "Shigella_flexneri",
#'     top_n_features = 25,
#'     network_top_n = 10
#'   )
#' }
exportAMRVisualizations <- function(output_dir = "amRviz_exports",
                                    formats = c("png", "pdf"),
                                    results_root = NULL,
                                    amrdata_root = NULL,
                                    species = NULL,
                                    top_n_features = 10,
                                    network_top_n = 5,
                                    width = 1200,
                                    height = 800,
                                    scale = 2,
                                    delay = 1.5,
                                    verbose = TRUE) {
  # Validate dependencies up front with actionable messages.
  for (pkg in c("htmlwidgets", "webshot2", "chromote", "magick")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(
        "Package '", pkg, "' is required for exportAMRVisualizations(). ",
        "Install it with install.packages('", pkg, "').",
        call. = FALSE
      )
    }
  }
  chrome_ok <- tryCatch(
    nzchar(chromote::find_chrome()),
    error = function(e) FALSE
  )
  if (!isTRUE(chrome_ok)) {
    stop(
      "A Chrome/Chromium browser is required to render the figures but none ",
      "was found. Install Google Chrome or Chromium (webshot2/chromote use it ",
      "to photograph the interactive plots).",
      call. = FALSE
    )
  }

  formats <- tolower(formats)
  formats[formats == "jpeg"] <- "jpg"
  valid <- c("png", "pdf", "jpg", "svg")
  bad <- setdiff(formats, valid)
  if (length(bad)) {
    stop(
      "Unsupported format(s): ", paste(bad, collapse = ", "),
      ". Choose from: ", paste(valid, collapse = ", "), ".",
      call. = FALSE
    )
  }

  if (!is.numeric(scale) || length(scale) != 1 || is.na(scale) || scale <= 0) {
    stop("`scale` must be a single positive number.", call. = FALSE)
  }

  for (nm in c("top_n_features", "network_top_n")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 1 || is.na(v) || v < 1) {
      stop("`", nm, "` must be a single positive number.", call. = FALSE)
    }
  }

  # Default amrdata_root: ~/amRdata/data when present (mirrors the dashboard).
  if (is.null(amrdata_root)) {
    default_amrdata <- file.path(path.expand("~"), "amRdata", "data")
    if (dir.exists(default_amrdata)) amrdata_root <- default_amrdata
  }

  # Load model results (user results_root, else packaged demo data).
  species_dirs <- if (!is.null(results_root) && nzchar(results_root)) {
    unname(listAmRmlSpeciesFolders(results_root))
  } else {
    NULL
  }
  perf_data <- loadMLResults(results_root, species_dirs, verbose = verbose)
  top_features <- loadTopFeat(results_root, species_dirs, verbose = verbose)

  # ML species: (code, folder label) pairs, excluding cross/MDR pseudo-species.
  ml_species <- list(code = character(0), label = character(0))
  if (nrow(perf_data) &&
    all(c("species", "species_label") %in% names(perf_data))) {
    # Drop rows with no species code: some parquets carry NA/blank species,
    # which would otherwise yield a phantom species whose folder collides with
    # (and overwrites) a real one. The dashboard drops these implicitly via
    # sort(); we do it explicitly.
    pairs <- perf_data |>
      dplyr::filter(!is.na(.data$species) & nzchar(.data$species)) |>
      dplyr::filter(!(.data$species %in% c("cross", "MDR"))) |>
      dplyr::distinct(.data$species, .data$species_label) |>
      dplyr::arrange(.data$species_label, .data$species)
    ml_species <- list(
      code = as.character(pairs$species),
      label = as.character(pairs$species_label)
    )
  }

  # Metadata species: subdirectories that hold a metadata.parquet.
  scan_root <- if (!is.null(results_root) && nzchar(results_root)) {
    results_root
  } else {
    system.file("extdata", package = "amRviz")
  }
  meta_species <- character(0)
  if (nzchar(scan_root)) {
    for (d in list.dirs(scan_root, full.names = TRUE, recursive = FALSE)) {
      if (file.exists(file.path(d, "metadata.parquet"))) {
        meta_species <- c(meta_species, basename(d))
      }
    }
    meta_species <- sort(unique(meta_species))
  }

  # Optional species filter (matched against folder labels).
  if (!is.null(species)) {
    keep <- ml_species$label %in% species | ml_species$code %in% species
    ml_species <- list(
      code = ml_species$code[keep], label = ml_species$label[keep]
    )
    meta_species <- meta_species[meta_species %in% species]
  }

  if (!length(ml_species$code) && !length(meta_species)) {
    stop(
      "No species found to export. Check `results_root` / `species`, or omit ",
      "them to use the packaged demo data.",
      call. = FALSE
    )
  }

  specs <- .exportPlanSpecs(
    perf_data, top_features, ml_species, meta_species,
    results_root, amrdata_root,
    top_n_features = top_n_features, network_top_n = network_top_n
  )

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  if (verbose) {
    message(
      "Exporting ", length(specs), " figures (",
      paste(formats, collapse = ", "), ") to ",
      normalizePath(output_dir, mustWork = FALSE)
    )
  }

  results <- vector("list", length(specs))
  for (i in seq_along(specs)) {
    spec <- specs[[i]]
    group_dir <- file.path(output_dir, spec$group)
    if (!dir.exists(group_dir)) dir.create(group_dir, recursive = TRUE)
    path_base <- file.path(group_dir, spec$name)

    if (verbose) {
      message(
        "  [", i, "/", length(specs), "] ",
        spec$group, "/", spec$name
      )
    }

    written <- tryCatch(
      {
        widget <- spec$build()
        .exportWidgetFile(
          widget, path_base, formats,
          width = spec$width %||% width,
          height = spec$height %||% height,
          scale = scale, trim = isTRUE(spec$trim),
          delay = delay, verbose = verbose
        )
      },
      error = function(e) {
        if (verbose) message("    build failed: ", conditionMessage(e))
        character(0)
      }
    )

    results[[i]] <- data.frame(
      group = spec$group, name = spec$name,
      written = length(written), ok = length(written) > 0,
      stringsAsFactors = FALSE
    )
  }

  summary_df <- do.call(rbind, results)
  if (verbose) {
    n_ok <- sum(summary_df$ok)
    n_files <- sum(summary_df$written)
    message(
      "Done: ", n_ok, "/", nrow(summary_df),
      " figures rendered, ", n_files, " files written."
    )
    failed <- summary_df[!summary_df$ok, , drop = FALSE]
    if (nrow(failed)) {
      message(
        "  No output for: ",
        paste(
          paste0(failed$group, "/", failed$name),
          collapse = ", "
        )
      )
    }
  }

  invisible(summary_df)
}
