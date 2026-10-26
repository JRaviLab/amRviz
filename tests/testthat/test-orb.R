## Tests for ORB discovery (R/utils_orb.R)

# Build {tmp}/{species}/orb/ with a manifest and (optionally) metadata.parquet.
make_orb <- function(root, species = "Haemophilus_influenzae", user_bacs = NULL,
                     metadata = TRUE, manifest_type = "amR_dataset",
                     schema_version = 1L, dataset_id = "Hin") {
  orb_dir <- file.path(root, species, "orb")
  dir.create(orb_dir, recursive = TRUE)
  if (metadata) {
    arrow::write_parquet(
      tibble::tibble(genome.genome_id = "1.1", drug_class = "x"),
      file.path(orb_dir, "metadata.parquet")
    )
  }
  manifest <- list(
    schema_version = schema_version,
    manifest_type = manifest_type,
    dataset_id = dataset_id,
    dataset = list(selection = list(user_bacs = as.list(user_bacs)))
  )
  mp <- file.path(orb_dir, "manifest_run_20260101T000000000_pid1.json")
  jsonlite::write_json(manifest, mp, auto_unbox = TRUE, null = "null")
  mp
}

test_that(".list_available_orbs describes a valid ORB", {
  root <- withr::local_tempdir()
  mp <- make_orb(root, user_bacs = "Haemophilus influenzae")

  orbs <- .list_available_orbs(mp)
  expect_equal(nrow(orbs), 1)
  expect_equal(orbs$label, "Haemophilus influenzae")
  expect_equal(orbs$species_dir, "Haemophilus_influenzae")
  expect_equal(orbs$dataset_id, "Hin")
  expect_true(file.exists(orbs$metadata_path))
  expect_equal(basename(orbs$orb_dir), "orb")
})

test_that(".list_available_orbs falls back to the folder name for the label", {
  root <- withr::local_tempdir()
  mp <- make_orb(root, species = "Shigella_sonnei")
  expect_equal(.list_available_orbs(mp)$label, "Shigella sonnei")
})

test_that(".list_available_orbs drops stale and invalid registry entries", {
  root <- withr::local_tempdir()
  good <- make_orb(root, species = "Good_bug")
  no_meta <- make_orb(root, species = "No_metadata", metadata = FALSE)
  wrong_type <- make_orb(root, species = "Wrong_type", manifest_type = "other")
  wrong_ver <- make_orb(root, species = "Wrong_ver", schema_version = 2L)
  garbage <- file.path(root, "garbage.json")
  writeLines("not json", garbage)
  flat <- file.path(root, "Flat", "manifest_run_x.json")
  dir.create(dirname(flat))
  file.copy(good, flat)

  orbs <- .list_available_orbs(c(
    good, no_meta, wrong_type, wrong_ver, garbage, flat,
    file.path(root, "missing.json"), NA_character_, ""
  ))
  expect_equal(orbs$species_dir, "Good_bug")
})

test_that(".list_available_orbs returns an empty tibble when nothing is found", {
  orbs <- .list_available_orbs(character(0))
  expect_s3_class(orbs, "tbl_df")
  expect_equal(nrow(orbs), 0)
  expect_named(orbs, names(.empty_orbs()))
})

test_that(".list_available_orbs dedupes manifests pointing at the same ORB", {
  root <- withr::local_tempdir()
  mp <- make_orb(root)
  newer <- file.path(dirname(mp), "manifest_run_20260202T000000000_pid1.json")
  file.copy(mp, newer)
  Sys.setFileTime(mp, Sys.time() - 3600)

  orbs <- .list_available_orbs(c(mp, newer))
  expect_equal(nrow(orbs), 1)
  expect_equal(basename(orbs$manifest_path), basename(newer))
})

test_that(".orb_choices disambiguates identical labels", {
  root_a <- withr::local_tempdir()
  root_b <- withr::local_tempdir()
  a <- make_orb(root_a, user_bacs = "Haemophilus influenzae")
  b <- make_orb(root_b, user_bacs = "Haemophilus influenzae")
  c3 <- make_orb(root_a, species = "Shigella_sonnei", user_bacs = "Shigella sonnei")

  choices <- .orb_choices(.list_available_orbs(c(a, b, c3)))
  expect_length(choices, 3)
  expect_false(anyDuplicated(names(choices)) > 0)
  expect_true(all(vapply(choices, .is_orb_id, logical(1))))
  expect_true("Shigella sonnei" %in% names(choices))
  expect_equal(.orb_choices(.empty_orbs()), character(0))
})

test_that("ORB ids resolve to metadata path and species name", {
  root <- withr::local_tempdir()
  mp <- make_orb(root)
  orbs <- .list_available_orbs(mp)
  id <- .orb_id(orbs$orb_dir)

  expect_true(.is_orb_id(id))
  expect_false(.is_orb_id("Shigella_flexneri"))
  expect_false(.is_orb_id(NULL))
  expect_equal(.orb_metadata_path(id), orbs$metadata_path)
  expect_equal(.metadata_species_name(id), "Haemophilus_influenzae")
  expect_equal(.metadata_species_name("Shigella_flexneri"), "Shigella_flexneri")
  expect_null(.orb_metadata_path(.orb_id(file.path(root, "nope"))))
  expect_null(.orb_metadata_path("Shigella_flexneri"))
})

test_that(".species_dir_has_metadata sees flat and orb/ layouts", {
  root <- withr::local_tempdir()
  make_orb(root, species = "Haemophilus_influenzae")
  flat <- file.path(root, "Flat")
  dir.create(flat)
  file.create(file.path(flat, "metadata.parquet"))

  expect_true(.species_dir_has_metadata(file.path(root, "Haemophilus_influenzae")))
  expect_true(.species_dir_has_metadata(flat))
  expect_false(.species_dir_has_metadata(file.path(root, "Nope")))
})
