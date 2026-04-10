## Tests for getHoldoutsDrugChoices

test_that("getHoldoutsDrugChoices returns empty for NULL input", {
  expect_equal(getHoldoutsDrugChoices(NULL), character(0))
})

test_that("getHoldoutsDrugChoices returns empty for empty data frame", {
  expect_equal(getHoldoutsDrugChoices(data.frame()), character(0))
})

test_that("getHoldoutsDrugChoices returns empty for zero-row tibble", {
  df <- tibble::tibble(
    strat_label = character(0),
    drug_or_class = character(0),
    species = character(0)
  )
  expect_equal(getHoldoutsDrugChoices(df), character(0))
})

test_that("getHoldoutsDrugChoices filters to stratified rows only", {
  df <- tibble::tibble(
    strat_label = c(NA_character_, "", "country", "year"),
    drug_or_class = c("amoxicillin", "ampicillin", "ciprofloxacin", "tetracycline"),
    species = c("Sau", "Sau", "Sau", "Sau")
  )
  result <- getHoldoutsDrugChoices(df)
  # Only rows with non-empty strat_label should be included
  expect_true("ciprofloxacin" %in% result)
  expect_true("tetracycline" %in% result)
  expect_false("amoxicillin" %in% result)
  expect_false("ampicillin" %in% result)
})

test_that("getHoldoutsDrugChoices returns sorted unique values", {
  df <- tibble::tibble(
    strat_label = c("country", "country", "year"),
    drug_or_class = c("tetracycline", "ampicillin", "tetracycline"),
    species = c("Sau", "Sau", "Sau")
  )
  result <- getHoldoutsDrugChoices(df)
  expect_equal(result, c("ampicillin", "tetracycline"))
})

test_that("getHoldoutsDrugChoices filters by bug when provided", {
  df <- tibble::tibble(
    strat_label = c("country", "country", "country"),
    drug_or_class = c("ampicillin", "ciprofloxacin", "tetracycline"),
    species = c("Sau", "Kpn", "Sau")
  )
  result <- getHoldoutsDrugChoices(df, bug = "Sau")
  expect_true("ampicillin" %in% result)
  expect_true("tetracycline" %in% result)
  expect_false("ciprofloxacin" %in% result)
})

test_that("getHoldoutsDrugChoices normalizes species for comparison", {
  df <- tibble::tibble(
    strat_label = c("country", "country"),
    drug_or_class = c("ampicillin", "ciprofloxacin"),
    species = c("Esp.", "Esp")
  )
  # Both "Esp." and "Esp" should match when bug = "Esp."

  result <- getHoldoutsDrugChoices(df, bug = "Esp.")
  expect_length(result, 2)
})

test_that("getHoldoutsDrugChoices excludes NA and blank drug_or_class", {
  df <- tibble::tibble(
    strat_label = c("country", "country", "country"),
    drug_or_class = c("ampicillin", NA_character_, "  "),
    species = c("Sau", "Sau", "Sau")
  )
  result <- getHoldoutsDrugChoices(df)
  expect_equal(result, "ampicillin")
})
