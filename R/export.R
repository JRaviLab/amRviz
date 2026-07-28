# Headless batch export of every dashboard visualization to static image
# files. Lets a user render the full set of amRviz figures from the packaged
# demo data (or their own amRml results) without ever launching the Shiny app.


#' Snapshot one htmlwidget to static image file(s)
#'
#' Every amRviz plot returns a plotly or networkD3 htmlwidget. There is no
#' server-side raster renderer for these, so we save the widget to a local HTML
#' bundle and photograph it with a headless Chrome via webshot2. png/pdf/jpg are
#' all handled this way (the format follows the file extension). svg is
#' best-effort: it needs the plotly image engine (kaleido) and only applies to
#' plotly charts, so it is skipped silently when unavailable.
#'
#' @param widget An htmlwidget (plotly / networkD3), or NULL.
#' @param path_base Output path without extension; one file per requested format
#'   is written as `path_base.<ext>`.
#' @param formats Character vector of extensions among png/pdf/jpg/jpeg/svg.
#' @param width,height Snapshot viewport size in pixels.
#' @param scale Device-pixel multiplier for raster (png/jpg) output: the saved
#'   image is `width * scale` by `height * scale` pixels. Higher values give
#'   sharper, higher-resolution figures (`scale = 2` ~ 340 dpi at a 7 in width;
#'   `scale = 4` ~ 680 dpi). Does not affect the vector PDF or SVG output.
#' @param delay Seconds to let the widget's JavaScript render before the shot.
#' @param verbose Whether to message on per-format failures.
#' @return Character vector of files actually written (possibly empty).
#' @keywords internal
#' @noRd
.exportWidgetFile <- function(widget, path_base, formats,
                              width = 1200, height = 800, scale = 2,
                              delay = 1.5, verbose = TRUE) {
  if (is.null(widget)) {
    return(character(0))
  }
  written <- character(0)
  raster <- intersect(formats, c("png", "pdf", "jpg", "jpeg"))
  want_svg <- "svg" %in% formats

  # png/pdf/jpg: save the widget once, then photograph it in each format.
  if (length(raster)) {
    tmpdir <- tempfile("amrviz_widget_")
    dir.create(tmpdir)
    on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)
    html <- file.path(tmpdir, "widget.html")
    saved <- tryCatch(
      {
        htmlwidgets::saveWidget(widget, html, selfcontained = FALSE)
        TRUE
      },
      error = function(e) {
        if (verbose) message("    saveWidget failed: ", conditionMessage(e))
        FALSE
      }
    )
    if (isTRUE(saved)) {
      for (ext in raster) {
        out <- paste0(path_base, ".", ext)
        # Supersample raster output for resolution; PDF is vector, so leave its
        # zoom at 1 (zooming would only rescale the page, not sharpen it).
        zoom <- if (ext == "pdf") 1 else scale
        ok <- tryCatch(
          {
            webshot2::webshot(
              html, out,
              vwidth = width, vheight = height,
              zoom = zoom, delay = delay, quiet = TRUE
            )
            file.exists(out) && file.info(out)$size > 0
          },
          error = function(e) {
            if (verbose) {
              message("    ", ext, " failed: ", conditionMessage(e))
            }
            FALSE
          }
        )
        if (isTRUE(ok)) written <- c(written, out)
      }
    }
  }

  # svg: plotly-only, via the kaleido engine when present.
  if (isTRUE(want_svg)) {
    out <- paste0(path_base, ".svg")
    ok <- tryCatch(
      {
        if (inherits(widget, "plotly")) {
          plotly::save_image(widget, out, width = width, height = height)
          file.exists(out) && file.info(out)$size > 0
        } else {
          FALSE
        }
      },
      error = function(e) FALSE
    )
    if (isTRUE(ok)) {
      written <- c(written, out)
    } else if (verbose) {
      message("    svg skipped (needs plotly + kaleido image engine)")
    }
  }

  written
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
  add <- function(group, name, build) {
    specs[[length(specs) + 1]] <<- list(
      group = group, name = name, build = build
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

    # Drug-feature network.
    local({
      code <- code
      add(folder, "drug_feature_network", function() {
        makeDrugFeatureNetwork(
          top_features, code,
          top_n = network_top_n,
          include_clusters = FALSE, include_cogs = FALSE,
          results_root = results_root
        )
      })
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
#' Renders the complete set of amRviz figures - metadata distributions, model
#' performance, feature importance, cross-model holdouts, and drug-feature
#' networks - to static image files, without launching the interactive Shiny
#' dashboard. This lets a user install the package, point it at model results
#' (or use the packaged demo data), and obtain figures for every panel in one
#' call.
#'
#' One figure set is produced per species using the same default selections the
#' dashboard opens with (e.g. all molecular scales, binary encoding, gentamicin
#' where present). The number of features shown is adjustable via
#' `top_n_features` (feature-importance panels) and `network_top_n` (the
#' drug-feature network). Species-agnostic overviews (the performance heatmaps
#' and the cross-species feature-importance panel) are written once under
#' `_overview` and `_across_species`.
#'
#' Because every dashboard plot is an interactive htmlwidget (plotly or
#' networkD3), static export photographs each widget with a headless Chrome via
#' the \pkg{webshot2} package. `png`, `pdf`, and `jpg` are fully supported.
#' `svg` is best-effort: it requires the plotly image engine (kaleido) and is
#' silently skipped for widgets or environments where that is unavailable.
#'
#' On output quality: `pdf` is written as a true vector figure (Chrome's
#' Skia PDF backend over the plots' underlying SVG), so it is resolution
#' independent and the best choice for publication. Raster formats (`png`,
#' `jpg`) are screenshots whose resolution is `width * scale` by
#' `height * scale` pixels; raise `scale` for high-DPI raster figures.
#'
#' @param output_dir Directory to write figures into; created if needed. Files
#'   are organised as `output_dir/<species>/<panel>.<ext>`, with global
#'   overviews under `output_dir/_overview` and `output_dir/_across_species`.
#' @param formats Character vector of output formats, any of `"png"`, `"pdf"`,
#'   `"jpg"`, `"svg"`. Defaults to `c("png", "pdf")`.
#' @param results_root Path to a directory of amRml model outputs (per-species
#'   subdirectories of `*_perf.parquet` / `*_top_features.parquet` /
#'   `metadata.parquet`). When `NULL` (default), the packaged demo data bundled
#'   with amRviz is used.
#' @param amrdata_root Path to amRdata annotation parquets used to enrich
#'   feature-importance panels (COG categories, etc.). When `NULL` (default),
#'   `~/amRdata/data` is used if present; otherwise annotation-based panels fall
#'   back to unenriched output.
#' @param species Optional character vector restricting which species are
#'   exported, matched against the species folder names. `NULL` (default)
#'   exports every species found.
#' @param top_n_features Number of top features to show in each
#'   feature-importance panel (per model), matching the dashboard's "Top
#'   features" control. Defaults to `10`.
#' @param network_top_n Number of top features per drug to include in the
#'   drug-feature network, matching the dashboard's network slider. Defaults to
#'   `5`.
#' @param width,height Snapshot viewport size in pixels.
#' @param scale Device-pixel multiplier for raster (`png`/`jpg`) output; the
#'   saved image is `width * scale` by `height * scale` pixels. Defaults to `2`
#'   (~340 dpi at a 7 in figure width); use `3`-`4` for ~500-680 dpi. Ignored
#'   for the vector `pdf` and `svg` output.
#' @param delay Seconds to wait for each widget's JavaScript to render before
#'   the screenshot is taken. Increase if figures come out partially rendered.
#' @param verbose Whether to print per-figure progress.
#'
#' @return Invisibly, a data frame with one row per attempted figure: its
#'   `group`, `name`, the number of files `written`, and whether it `ok`.
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
  for (pkg in c("htmlwidgets", "webshot2")) {
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
          width = width, height = height, scale = scale,
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
