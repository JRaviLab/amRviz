## Tests for the Metadata tab's tables (R/metadata_plots.R). Each plot draws
## one table, so these pin down what each table counts.

# Two genomes tested against two drugs each, plus a third genome tested
# against one: 5 genome-drug records. Free text is typed inconsistently, as in
# BV-BRC.
meta_rows <- function() {
  tibble::tibble(
    species = "Shigella sonnei",
    genome.genome_id = c("G1", "G1", "G2", "G2", "G3"),
    genome_drug.antibiotic = c("ampicillin", "tetracycline", "ampicillin", "tetracycline", "ampicillin"),
    drug_class = c("penicillins", "tetracyclines", "penicillins", "tetracyclines", "penicillins"),
    genome_drug.resistant_phenotype = c("Resistant", "Susceptible", "Resistant", "Resistant", "Susceptible"),
    genome.isolation_country = c("France", "France", "France", "France", ""),
    genome.collection_year = c(2015, 2015, 2015, 2015, NA),
    genome.host_common_name = c("Human", "Human", "human ", "human ", ""),
    genome.isolation_source = c("Feces", "Feces", "feces", "feces", NA)
  )
}

test_that("metadataSummary counts records, genomes and phenotypes", {
  s <- metadataSummary(meta_rows())
  expect_equal(s$records, 5)
  expect_equal(s$genomes, 3)
  expect_equal(s$drugs, 2)
  expect_equal(c(s$resistant_records, s$susceptible_records), c(3, 2))
  expect_equal(s$top_countries, "France")
  expect_equal(s$species, "Shigella sonnei")
})

test_that("metadataAvailability counts genomes per drug and phenotype", {
  tbl <- metadataAvailability(meta_rows())
  expect_named(tbl, c("drug", "phenotype", "genomes"))
  amp <- tbl[tbl$drug == "ampicillin", ]
  expect_equal(amp$genomes[amp$phenotype == "Resistant"], 2)
})

test_that("metadataByCountry counts genomes as well as records", {
  tbl <- metadataByCountry(meta_rows())
  expect_named(tbl, c("country", "iso3", "genomes", "records"))
  # The blank country is left out.
  expect_equal(tbl$country, "France")
  expect_equal(c(tbl$genomes, tbl$records), c(2, 4))
  expect_equal(tbl$iso3, "FRA")
  expect_true(is.na(metadataByCountry(
    dplyr::mutate(meta_rows(), genome.isolation_country = "Southeast Asia")
  )$iso3))
})

test_that("metadataByYear has one row per year and phenotype", {
  tbl <- metadataByYear(meta_rows())
  expect_named(tbl, c("year", "phenotype", "records"))
  # Drugs and countries are summed, so each line has one point per year.
  expect_equal(nrow(tbl), 2)
  expect_equal(tbl$records[tbl$phenotype == "Resistant"], 3)
})

test_that("hosts and sources are compared regardless of case and spacing", {
  hosts <- metadataByHost(meta_rows())
  expect_equal(unique(hosts$host), "human")
  expect_equal(sum(hosts$records), 4)

  sources <- metadataBySource(meta_rows())
  expect_equal(levels(sources$source), c("feces", "not recorded"))
  # Every record is counted, including those without a host.
  expect_equal(sum(sources$records), 5)
})

test_that("metadataBySource groups rare sources as other, apart from missing ones", {
  meta <- meta_rows()
  meta$genome.isolation_source <- c("a", "a", "b", "c", NA)
  sources <- metadataBySource(meta, top_n = 1)
  expect_equal(levels(sources$source), c("a", "other", "not recorded"))
})

test_that("metadataSankey keeps a value that appears in two tiers as two nodes", {
  meta <- meta_rows()
  meta$genome.isolation_source <- meta$genome.host_common_name
  flow <- metadataSankey(meta)
  expect_named(flow$nodes, c("id", "label", "tier"))
  expect_equal(sum(flow$nodes$label == "human"), 2)
  expect_true(all(c(flow$links$from, flow$links$to) %in% flow$nodes$id))
  expect_null(metadataSankey(meta[0, ]))
})

test_that("the tab's plots draw their tables", {
  meta <- meta_rows()
  expect_s3_class(makeDataAvailabilityPlot(metadataAvailability(meta)), "plotly")
  expect_s3_class(makeGeoChloroPlot(metadataByCountry(meta)), "plotly")
  expect_s3_class(makeTimeSeriesAMRPlot(metadataByYear(meta), "all"), "plotly")
  expect_s3_class(makeHostIsolatePlot(metadataByHost(meta)), "plotly")
  expect_s3_class(makeIsolationSourcesPlot(metadataBySource(meta)), "plotly")
  expect_s3_class(makeMetadataSankey(meta), "sankeyNetwork")
  # Empty tables show a message instead of an empty chart.
  expect_s3_class(makeHostIsolatePlot(metadataByHost(meta[0, ])), "plotly")
})
