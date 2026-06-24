## Tests for makeQuickStats

test_that("makeQuickStats returns a tagList with expected structure", {
  df <- tibble::tibble(
    genome.genome_id = c("G1", "G1", "G2", "G3"),
    genome_drug.antibiotic = c(
      "ampicillin", "tetracycline",
      "ampicillin", "ciprofloxacin"
    ),
    genome_drug.resistant_phenotype = c(
      "Resistant", "Susceptible",
      "Resistant", "Susceptible"
    ),
    drug_class = c(
      "penicillins", "tetracyclines",
      "penicillins", "fluoroquinolones"
    ),
    genome.isolation_country = c("USA", "UK", "USA", "Germany"),
    species = c(
      "staphylococcus aureus", "staphylococcus aureus",
      "staphylococcus aureus", "staphylococcus aureus"
    )
  )

  result <- makeQuickStats(df)
  expect_s3_class(result, "shiny.tag.list")
})

test_that("makeQuickStats handles single-row data", {
  df <- tibble::tibble(
    genome.genome_id = "G1",
    genome_drug.antibiotic = "ampicillin",
    genome_drug.resistant_phenotype = "Resistant",
    drug_class = "penicillins",
    genome.isolation_country = "USA",
    species = "staphylococcus aureus"
  )

  result <- makeQuickStats(df)
  expect_s3_class(result, "shiny.tag.list")
})
