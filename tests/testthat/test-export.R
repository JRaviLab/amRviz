## Tests for the headless visualization exporter:
## exportAMRVisualizations() and its internal plan/snapshot helpers.

# ── .exportPlanSpecs (no browser needed) ─────────────────────────────────────

test_that(".exportPlanSpecs builds figure specs from the bundled demo data", {
  perf <- loadMLResults(verbose = FALSE)
  top <- loadTopFeat(verbose = FALSE)
  skip_if(!nrow(perf), "No demo performance data available")

  pairs <- perf |>
    dplyr::filter(!(species %in% c("cross", "MDR"))) |>
    dplyr::distinct(species, species_label)
  ml_species <- list(
    code = as.character(pairs$species),
    label = as.character(pairs$species_label)
  )
  extdata <- system.file("extdata", package = "amRviz")
  meta_species <- basename(
    list.dirs(extdata, full.names = TRUE, recursive = FALSE)
  )
  meta_species <- meta_species[vapply(
    file.path(extdata, meta_species, "metadata.parquet"),
    file.exists, logical(1)
  )]

  specs <- .exportPlanSpecs(
    perf, top, ml_species, meta_species,
    results_root = NULL, amrdata_root = NULL
  )

  expect_gt(length(specs), 0)
  # Every spec is addressable and lazily builds a widget.
  for (s in specs) {
    expect_true(all(c("group", "name", "build") %in% names(s)))
    expect_true(is.function(s$build))
    expect_true(nzchar(s$group) && nzchar(s$name))
  }
  # Global overviews are present exactly once each.
  groups_names <- vapply(
    specs, function(s) paste(s$group, s$name), character(1)
  )
  expect_true(any(grepl("_overview mcc_strip", groups_names)))
  expect_true(any(grepl("_overview mcc_heatmap", groups_names)))

  # No duplicate group/name specs: a phantom NA/blank species code shares a
  # species_label folder with a real one and would otherwise overwrite its
  # figures (regression guard).
  expect_equal(anyDuplicated(groups_names), 0L)
})

test_that("exportAMRVisualizations drops species with no code (no phantom folder)", {
  perf <- loadMLResults(verbose = FALSE)
  skip_if(!nrow(perf), "No demo performance data available")
  # The demo Shigella_sonnei parquets carry some rows with a blank/NA species
  # code; those must not produce their own export iteration.
  skip_if(!any(is.na(perf$species) | !nzchar(perf$species)),
    "Demo data has no blank species codes to exercise the guard"
  )
  top <- loadTopFeat(verbose = FALSE)
  # Deliberately build ml_species the naive way, INCLUDING the phantom
  # NA/blank-code entry, to prove .exportPlanSpecs drops it defensively.
  pairs <- perf |>
    dplyr::filter(!(species %in% c("cross", "MDR"))) |>
    dplyr::distinct(species, species_label)
  ml_species <- list(
    code = as.character(pairs$species),
    label = as.character(pairs$species_label)
  )
  expect_true(any(is.na(ml_species$code) | !nzchar(ml_species$code)))

  specs <- .exportPlanSpecs(
    perf, top, ml_species, character(0),
    results_root = NULL, amrdata_root = NULL
  )
  groups_names <- vapply(
    specs, function(s) paste(s$group, s$name), character(1)
  )
  expect_equal(anyDuplicated(groups_names), 0L)
})

# ── format validation (no browser needed when it errors early) ───────────────

test_that("exportAMRVisualizations rejects unsupported formats", {
  skip_if_not_installed("webshot2")
  skip_if_not_installed("htmlwidgets")
  skip_if(
    tryCatch(!nzchar(chromote::find_chrome()), error = function(e) TRUE),
    "No Chrome/Chromium available"
  )
  expect_error(
    exportAMRVisualizations(
      output_dir = tempfile("amrviz_bad"),
      formats = "tiff", verbose = FALSE
    ),
    "Unsupported format"
  )
})

# ── end-to-end export (needs a headless browser) ─────────────────────────────

test_that("exportAMRVisualizations writes files for the demo data", {
  skip_on_cran()
  skip_if_not_installed("webshot2")
  skip_if_not_installed("htmlwidgets")
  skip_if(
    tryCatch(!nzchar(chromote::find_chrome()), error = function(e) TRUE),
    "No Chrome/Chromium available"
  )

  out <- tempfile("amrviz_export")
  on.exit(unlink(out, recursive = TRUE), add = TRUE)

  res <- exportAMRVisualizations(
    output_dir = out,
    formats = "png",
    species = "Shigella_flexneri",
    verbose = FALSE
  )

  expect_s3_class(res, "data.frame")
  expect_true(all(c("group", "name", "written", "ok") %in% names(res)))
  expect_gt(nrow(res), 0)
  # At least the core per-species figures should render.
  expect_true(any(res$ok))
  # No two figures share a group/name (which would mean colliding output paths).
  expect_equal(anyDuplicated(paste(res$group, res$name)), 0L)
  png_files <- list.files(out, pattern = "\\.png$", recursive = TRUE)
  expect_gt(length(png_files), 0)
})
