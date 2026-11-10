## Tests for ORB discovery (R/utils_orb.R)

test_that(".list_available_orbs describes a valid ORB", {
  root <- withr::local_tempdir()
  mp <- make_orb(root, user_bacs = "Haemophilus influenzae")

  orbs <- .list_available_orbs(mp)
  expect_equal(nrow(orbs), 1)
  expect_equal(orbs$label, "Haemophilus influenzae")
  expect_equal(orbs$dataset_id, "Hin")
  expect_true(file.exists(orbs$metadata_path))
})

test_that(".list_available_orbs reads the species from the manifest, not the folder", {
  root <- withr::local_tempdir()
  mp <- make_orb(root, folder = "Shigella_sonnei", dataset_id = "Sfl")
  # No user_bacs: the label falls back to the dataset id, not "Shigella sonnei".
  expect_equal(.list_available_orbs(mp)$label, "Sfl")
})

test_that(".list_available_orbs finds the metadata file the manifest names", {
  root <- withr::local_tempdir()
  mp <- make_orb(root, metadata_file = "genome_metadata.parquet")
  expect_equal(basename(.list_available_orbs(mp)$metadata_path), "genome_metadata.parquet")
})

test_that(".list_available_orbs resolves Windows paths against the manifest", {
  root <- withr::local_tempdir()
  mp <- make_orb(root)
  manifest <- jsonlite::read_json(mp)
  manifest$artifacts$amRml_input$metadata_parquet <- "C:\\Users\\x\\orb\\metadata.parquet"
  jsonlite::write_json(manifest, mp, auto_unbox = TRUE)
  expect_true(file.exists(.list_available_orbs(mp)$metadata_path))
})

test_that(".list_available_orbs drops stale and invalid registry entries", {
  root <- withr::local_tempdir()
  good <- make_orb(root, folder = "Good_bug", dataset_id = "Good")
  no_meta <- make_orb(root, folder = "No_metadata", metadata = FALSE)
  wrong_type <- make_orb(root, folder = "Wrong_type", manifest_type = "other")
  wrong_ver <- make_orb(root, folder = "Wrong_ver", schema_version = 2L)
  garbage <- file.path(root, "garbage.json")
  writeLines("not json", garbage)

  orbs <- .list_available_orbs(c(
    good, no_meta, wrong_type, wrong_ver, garbage,
    file.path(root, "missing.json"), NA_character_, ""
  ))
  expect_equal(orbs$dataset_id, "Good")
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

test_that(".find_orb_manifests finds ORBs at, inside and two levels below a root", {
  root <- withr::local_tempdir()
  deep <- make_orb(root, folder = "Haemophilus_influenzae")
  expect_equal(normalizePath(.find_orb_manifests(root)), normalizePath(deep))
  expect_length(.find_orb_manifests(dirname(deep)), 1)
  expect_length(.find_orb_manifests(file.path(root, "nope")), 0)
})

test_that(".orb_choices disambiguates identical labels", {
  root_a <- withr::local_tempdir()
  root_b <- withr::local_tempdir()
  a <- make_orb(root_a, user_bacs = "Haemophilus influenzae")
  b <- make_orb(root_b, user_bacs = "Haemophilus influenzae")
  c3 <- make_orb(root_a, folder = "Shigella_sonnei", user_bacs = "Shigella sonnei")

  choices <- .orb_choices(.list_available_orbs(c(a, b, c3)))
  expect_length(choices, 3)
  expect_false(anyDuplicated(names(choices)) > 0)
  expect_true(all(file.exists(choices)))
  expect_true("Shigella sonnei" %in% names(choices))
  expect_equal(.orb_choices(.empty_orbs()), character(0))
})

test_that("listMetadataSpecies lists ORBs by their manifest label", {
  root <- withr::local_tempdir()
  make_orb(root, folder = "Misleading_name", user_bacs = "Haemophilus influenzae")
  choices <- listMetadataSpecies(root)
  expect_named(choices, "Haemophilus influenzae")
  expect_true(file.exists(choices))
})
