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
  expect_named(tbl, c("drug", "drug_abbr", "phenotype", "genomes"))
  amp <- tbl[tbl$drug == "ampicillin", ]
  expect_equal(amp$genomes[amp$phenotype == "Resistant"], 2)
})

test_that("per-drug tables carry the drug abbreviation for axis labels", {
  meta <- dplyr::mutate(
    meta_rows(),
    drug_abbr = ifelse(.data$genome_drug.antibiotic == "ampicillin", "AMP", NA)
  )
  for (tbl in list(metadataAvailability(meta), metadataByHost(meta), metadataBySource(meta))) {
    expect_equal(unique(tbl$drug_abbr[tbl$drug == "ampicillin"]), "AMP")
    # A missing abbreviation falls back to the drug name.
    expect_equal(unique(tbl$drug_abbr[tbl$drug == "tetracycline"]), "tetracycline")
  }
  # Without the column, the drug name is used.
  expect_equal(metadataAvailability(meta_rows())$drug_abbr, metadataAvailability(meta_rows())$drug)
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


# ---- Resistance tables ----

# Three genomes. G1: resistant to ampicillin and gentamicin, susceptible to
# amoxicillin (same class as ampicillin). G2: susceptible to ampicillin,
# intermediate to amoxicillin. G3: resistant to three classes, no year.
resistance_rows <- function() {
  tibble::tibble(
    species = "Shigella sonnei",
    genome.genome_id = c("G1", "G1", "G1", "G2", "G2", "G3"),
    genome_drug.antibiotic = c("ampicillin", "amoxicillin", "gentamicin", "ampicillin", "amoxicillin", "ampicillin"),
    drug_class = c("penicillins", "penicillins", "aminoglycosides", "penicillins", "penicillins", "penicillins"),
    genome_drug.resistant_phenotype = c("Resistant", "Susceptible", "Resistant", "Susceptible", "Intermediate", "Resistant"),
    genome.isolation_country = c("France", "France", "France", "Spain", "Spain", ""),
    genome.collection_year = c(2015, 2015, 2015, 2016, 2016, NA),
    num_resistant_classes = c(2, 2, 2, 0, 0, 3),
    resistant_classes = c("AMG_PEN", "AMG_PEN", "AMG_PEN", "Susceptible", "Susceptible", "AMG_PEN_TET")
  )
}

test_that("metadataResistance counts records per drug", {
  tbl <- metadataResistance(resistance_rows(), "drug")
  expect_named(tbl, c("label", "resistant", "intermediate", "susceptible", "total", "pct_resistant"))
  amp <- tbl[tbl$label == "ampicillin", ]
  expect_equal(c(amp$resistant, amp$susceptible, amp$total), c(2, 1, 3))
  # Most resistant first.
  expect_equal(tbl$label[1], "gentamicin")
})

test_that("metadataResistance counts each genome once per drug class", {
  tbl <- metadataResistance(resistance_rows(), "class")
  pen <- tbl[tbl$label == "penicillins", ]
  # G1 resistant (ampicillin), G2 intermediate (amoxicillin), G3 resistant.
  expect_equal(c(pen$resistant, pen$intermediate, pen$susceptible, pen$total), c(2, 1, 0, 3))
  # No drug_class column: an empty table, not an error.
  expect_equal(nrow(metadataResistance(resistance_rows()[, -4], "class")), 0)
})

test_that("metadataDrugYear leaves out undated records and rare drugs", {
  tbl <- metadataDrugYear(resistance_rows(), min_records = 1)
  expect_named(tbl, c("drug", "year", "resistant", "total", "pct_resistant"))
  expect_false(anyNA(tbl$year))
  expect_equal(sum(tbl$total), 5)
  expect_equal(nrow(metadataDrugYear(resistance_rows(), min_records = 3)), 0)
})

test_that("metadataCountryYear counts genomes, not records", {
  tbl <- metadataCountryYear(resistance_rows())
  expect_equal(tbl$genomes[tbl$country == "France"], 1)
  expect_equal(sum(tbl$genomes), 2)
})

test_that("metadataMdr counts genomes by resistant classes", {
  tbl <- metadataMdr(resistance_rows())
  expect_equal(tbl$classes, 0:3)
  expect_equal(tbl$genomes, c(1, 0, 1, 1))
  expect_equal(tbl$mdr, c(FALSE, FALSE, FALSE, TRUE))
  expect_equal(nrow(metadataMdr(resistance_rows()[, 1:7])), 0)
})

test_that("metadataProfiles lists resistant-class combinations", {
  tbl <- metadataProfiles(resistance_rows())
  expect_equal(tbl$profile, c("AMG + PEN", "AMG + PEN + TET"))
  expect_equal(tbl$genomes, c(1, 1))
})

test_that("the resistance plots draw, and draw a message when empty", {
  m <- resistance_rows()
  empty <- m[0, ]
  for (p in list(
    makeResistanceSharePlot(metadataResistance(m)),
    makeDrugYearHeatmap(metadataDrugYear(m, min_records = 1)),
    makeCountryYearHeatmap(metadataCountryYear(m)),
    makeMdrPlot(metadataMdr(m)),
    makeProfilesPlot(metadataProfiles(m)),
    makeResistanceSharePlot(metadataResistance(empty)),
    makeDrugYearHeatmap(metadataDrugYear(empty)),
    makeMdrPlot(metadataMdr(empty)),
    makeProfilesPlot(metadataProfiles(empty))
  )) {
    expect_s3_class(p, "plotly")
  }
})
