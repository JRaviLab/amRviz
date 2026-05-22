test_that("normalize_species removes trailing dots", {
  expect_equal(normalize_species("Esp."), "Esp")
  expect_equal(normalize_species("Sau"), "Sau")
  expect_equal(normalize_species("Kpn."), "Kpn")
})

test_that("normalize_species handles vectors", {
  input <- c("Esp.", "Sau", "Kpn.", "Aba")
  expected <- c("Esp", "Sau", "Kpn", "Aba")
  expect_equal(normalize_species(input), expected)
})

test_that("normalize_species preserves NA", {
  expect_equal(normalize_species(NA), NA_character_)
  result <- normalize_species(c("Sau", NA, "Esp."))
  expect_equal(result, c("Sau", NA, "Esp"))
})

test_that("normalize_species handles empty strings", {
  expect_equal(normalize_species(""), "")
  expect_equal(normalize_species(character(0)), character(0))
})

test_that("normalize_species only removes trailing dot", {
  # Dot in the middle should remain
  expect_equal(normalize_species("E.coli"), "E.coli")
  # Multiple trailing dots: only last one removed
  expect_equal(normalize_species("Esp.."), "Esp.")
})
