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

test_that("metadataByCountry sums records per country, excluding blanks", {
  meta <- readSpeciesMetadata(sonnei())
  by_country <- metadataByCountry(meta)
  expect_named(by_country, c("genome.isolation_country", "count"))
  expect_false(any(by_country$genome.isolation_country == ""))
  expect_identical(
    sum(by_country$count),
    sum(meta$genome.isolation_country != "", na.rm = TRUE)
  )
})

test_that("metadataByYear drops rows without a collection year", {
  by_year <- metadataByYear(readSpeciesMetadata(sonnei()))
  expect_true("n" %in% names(by_year))
  expect_false(anyNA(by_year$genome.collection_year))
})

test_that("metadataWithHost drops blank hosts and adds evidence", {
  with_host <- metadataWithHost(readSpeciesMetadata(sonnei()))
  expect_false(any(with_host$genome.host_common_name == ""))
  expect_true("genome_drug.evidence" %in% names(with_host))
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
