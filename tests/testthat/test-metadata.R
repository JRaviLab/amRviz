## Tests for the Metadata tab: data functions and the tab's server module

# The packaged Shigella sonnei metadata, as listMetadataSpecies() offers it.
sonnei <- function() listMetadataSpecies()["Shigella sonnei"]

test_that("listMetadataSpecies finds the demo species with readable labels", {
  choices <- listMetadataSpecies()
  expect_identical(names(choices), "Shigella sonnei")
  expect_true(all(file.exists(choices)))
  expect_identical(basename(sonnei()), "metadata.parquet")
})

test_that("listMetadataSpecies returns empty for a folder with no metadata", {
  expect_length(listMetadataSpecies(withr::local_tempdir()), 0)
})

test_that("readSpeciesMetadata tags rows with the species label", {
  meta <- readSpeciesMetadata(sonnei())
  expect_gt(nrow(meta), 0)
  expect_true(all(meta$species == "Shigella sonnei"))
})

test_that("readSpeciesMetadata returns an empty tibble for a missing file", {
  missing <- c("Not a species" = file.path(withr::local_tempdir(), "x.parquet"))
  expect_identical(nrow(readSpeciesMetadata(missing)), 0L)
})

test_that("the demo metadata builds every table the tab draws", {
  meta <- readSpeciesMetadata(sonnei())
  expect_equal(sum(metadataByCountry(meta)$records), nrow(filter_valid_country(meta)))
  expect_lte(max(metadataByCountry(meta)$genomes), length(unique(meta$genome.genome_id)))
  expect_false(anyNA(metadataByYear(meta)$year))
  expect_equal(sum(metadataBySource(meta)$records), nrow(meta))
  expect_gt(nrow(metadataSankey(meta)$links), 0)
})

test_that("topDrugClasses returns the most-recorded classes first", {
  meta <- tibble::tibble(drug_class = c("A", "B", "B", "C", "C", "C", NA))
  expect_identical(topDrugClasses(meta, n = 2), c("C", "B"))
  expect_null(topDrugClasses(tibble::tibble(x = 1)))
})

test_that("metadataServer reads the selected species and renders its outputs", {
  shiny::testServer(metadataServer, args = list(results_root = NULL), {
    session$setInputs(bug_metadata_id = unname(sonnei()))
    expect_true(all(meta()$species == "Shigella sonnei"))
    expect_match(
      output$quick_metadata_stats$html, "Data summary for Shigella sonnei"
    )
    expect_true(nzchar(output$geo_isolate_plot))
  })
})
