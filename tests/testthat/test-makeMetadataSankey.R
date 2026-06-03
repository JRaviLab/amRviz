## Tests for makeMetadataSankey

test_that("makeMetadataSankey returns NULL for NULL input", {
  expect_null(makeMetadataSankey(NULL))
})

test_that("makeMetadataSankey returns NULL for non-data.frame input", {
  expect_null(makeMetadataSankey("not a data frame"))
  expect_null(makeMetadataSankey(list(a = 1)))
})

test_that("makeMetadataSankey returns NULL for zero-row data", {
  df <- tibble::tibble()
  expect_null(makeMetadataSankey(df))
})

test_that("makeMetadataSankey returns NULL when required columns are missing", {
  df <- tibble::tibble(
    genome_drug.antibiotic = "ampicillin",
    drug_class = "beta-lactams"
  )
  expect_null(makeMetadataSankey(df))
})

test_that("makeMetadataSankey returns NULL when all rows are filtered out", {
  df <- tibble::tibble(
    genome_drug.resistant_phenotype = c("Indeterminate", "Intermediate"),
    genome_drug.antibiotic = c("ampicillin", "tetracycline"),
    genome.isolation_country = c("USA", "Canada"),
    genome.host_common_name = c("Human", "Bovine"),
    genome.isolation_source = c("blood", "feces"),
    drug_class = c("beta-lactams", "tetracyclines")
  )
  expect_null(makeMetadataSankey(df))
})

test_that("makeMetadataSankey returns a Sankey htmlwidget for valid data", {
  skip_if_not_installed("networkD3")

  df <- tibble::tibble(
    genome_drug.resistant_phenotype = rep(c("Resistant", "Susceptible"), 3),
    genome_drug.antibiotic = rep(c("ampicillin", "tetracycline"), 3),
    genome.isolation_country = rep(c("USA", "Canada", "UK"), 2),
    genome.host_common_name = rep(c("Human", "Bovine"), 3),
    genome.isolation_source = rep(c("blood", "feces", "urine"), 2),
    drug_class = rep(c("beta-lactams", "tetracyclines"), 3)
  )

  result <- makeMetadataSankey(df)
  expect_s3_class(result, "htmlwidget")
})

test_that("makeMetadataSankey respects max_classes parameter", {
  skip_if_not_installed("networkD3")

  df <- tibble::tibble(
    genome_drug.resistant_phenotype = rep("Resistant", 6),
    genome_drug.antibiotic = rep("ampicillin", 6),
    genome.isolation_country = rep("USA", 6),
    genome.host_common_name = rep("Human", 6),
    genome.isolation_source = rep("blood", 6),
    drug_class = c("a", "a", "b", "b", "c", "d")
  )

  result <- makeMetadataSankey(df, max_classes = 2)
  expect_s3_class(result, "htmlwidget")
})
