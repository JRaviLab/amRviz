## Tests for reading ML results: amRml runs (R/utils_amrml.R) and the
## packaged demo. Folder and file names in these fixtures are
## deliberately misleading, so passing means nothing is read from them.

# ── amRml runs ───────────────────────────────────────────────────────────────

test_that("a run's model details come from its tasks and matrices", {
  root <- withr::local_tempdir()
  mp <- make_orb(root, folder = "Shigella_sonnei", user_bacs = "Haemophilus influenzae")
  run <- make_amrml_run(file.path(dirname(mp), "amRml", "cross_drug_LOO_MDR"))

  perf <- .load_amrml_perf(run, verbose = FALSE)
  expect_equal(nrow(perf), 4)
  vanilla <- perf[perf$task_id %in% c("t1", "t2"), ]
  expect_equal(vanilla$shuffled, c(FALSE, TRUE))
  expect_true(all(is.na(vanilla$strat_label)))
  expect_false(any(vanilla$cross_test))
  expect_equal(unique(vanilla$feature_type), "domains")
  expect_equal(unique(vanilla$drug_or_class), "AMP")

  cross <- perf[perf$task_id == "t3", ]
  expect_true(cross$cross_test)
  expect_equal(cross$strat_label, "country")
  expect_equal(c(cross$strat_value, cross$strat_value_test), c("USA", "GBR"))
  expect_equal(cross$drug_label, "drug_class")
  expect_equal(cross$feature_type, "genes")
  expect_equal(cross$feature_subtype, "counts")
})

test_that("the loaders leave out models fitted on shuffled labels", {
  root <- withr::local_tempdir()
  mp <- make_orb(root)
  run <- make_amrml_run(file.path(dirname(mp), "amRml", "2026-10-18_130824"))

  perf <- loadMLResults(root, run, verbose = FALSE)
  expect_equal(perf$task_id, c("t1", "t3"))
  expect_false(any(perf$shuffled))

  # Older results leave `shuffled` NA on some rows; those are kept.
  expect_equal(nrow(.real_labels(tibble::tibble(shuffled = c(NA, TRUE, FALSE)))), 2)
  expect_equal(nrow(.real_labels(tibble::tibble(mcc = 1))), 1)
})

test_that("a run's species and metadata come from its ORB's manifest", {
  root <- withr::local_tempdir()
  mp <- make_orb(root, folder = "Shigella_sonnei", user_bacs = "Haemophilus influenzae")
  run <- make_amrml_run(file.path(dirname(mp), "amRml", "2026-10-18_130824"))

  perf <- loadMLResults(root, run, verbose = FALSE)
  expect_equal(unique(perf$species), "Hin")
  expect_equal(unique(perf$species_label), "Haemophilus influenzae")
  expect_equal(
    normalizePath(unique(perf$metadata_path)),
    normalizePath(file.path(dirname(mp), "metadata.parquet"))
  )

  top <- loadTopFeat(root, run, verbose = FALSE)
  expect_equal(top$Variable, c("PF00001", "blaTEM"))
  expect_equal(top$cross_test, c(FALSE, TRUE))
  expect_equal(unique(top$species), "Hin")
})

test_that("a run's tables are found through its manifest", {
  root <- withr::local_tempdir()
  mp <- make_orb(root, user_bacs = "Haemophilus influenzae")
  run <- make_amrml_run(file.path(dirname(mp), "amRml", "run"))
  # Move the performance table and say so in the manifest.
  dir.create(file.path(run, "elsewhere"))
  file.rename(
    file.path(run, "results", "performance.parquet"),
    file.path(run, "elsewhere", "perf.parquet")
  )
  manifest <- jsonlite::read_json(file.path(run, "manifest.json"))
  manifest$artifacts$amRviz_input$results$performance <- "elsewhere/perf.parquet"
  jsonlite::write_json(manifest, file.path(run, "manifest.json"), auto_unbox = TRUE)

  expect_equal(nrow(loadMLResults(root, run, verbose = FALSE)), 2)
})

test_that("only folders with an amRml manifest inside an ORB are runs", {
  root <- withr::local_tempdir()
  mp <- make_orb(root, user_bacs = "Haemophilus influenzae")
  run <- make_amrml_run(file.path(dirname(mp), "amRml", "run"))

  # A run outside any ORB's amRml/ folder isn't found.
  make_amrml_run(file.path(root, "elsewhere", "run"))
  expect_equal(normalizePath(unname(listResultSources(root))), normalizePath(run))

  # A manifest of another type, or none at all, isn't a run.
  manifest <- jsonlite::read_json(file.path(run, "manifest.json"))
  manifest$manifest_type <- "amR_dataset"
  jsonlite::write_json(manifest, file.path(run, "manifest.json"), auto_unbox = TRUE)
  expect_length(listResultSources(root), 0)
  unlink(file.path(run, "manifest.json"))
  expect_length(listResultSources(root), 0)
})

test_that("feature importance keeps the most important features of each model", {
  root <- withr::local_tempdir()
  mp <- make_orb(root)
  importance <- tibble::tibble(
    task_id = rep(c("t1", "t3"), each = 150),
    Variable = paste0("f", seq_len(300)),
    Importance = c(seq(-1.5, 0, length.out = 150), seq_len(150) / 100),
    Sign = "POS"
  )
  run <- make_amrml_run(file.path(dirname(mp), "amRml", "run"),
    importance = importance
  )

  top <- loadTopFeat(root, run, verbose = FALSE)
  expect_equal(as.vector(table(top$task_id)), c(100, 100))
  # Ranked by absolute importance, so large negative values are kept.
  expect_true(all(abs(top$Importance[top$task_id == "t1"]) >= 0.5))
  expect_equal(min(top$Importance[top$task_id == "t3"]), 0.51)
})

test_that("leave-one-out modes, MDR targets and unknown modes are left out", {
  root <- withr::local_tempdir()
  mp <- make_orb(root)
  run <- make_amrml_run(file.path(dirname(mp), "amRml", "run"))
  matrices_path <- file.path(run, "metadata", "matrices.parquet")
  matrices <- arrow::read_parquet(matrices_path)
  matrices$mode_id <- c("logo", "a_future_mode")
  arrow::write_parquet(matrices, matrices_path)
  expect_message(
    perf <- loadMLResults(root, run, verbose = TRUE),
    "a_future_mode"
  )
  expect_equal(nrow(perf), 0)

  matrices$mode_id <- "vanilla"
  matrices$unit <- c("mdr", "drug")
  arrow::write_parquet(matrices, matrices_path)
  expect_equal(unique(loadMLResults(root, run, verbose = FALSE)$drug_label), "drug")
})

test_that("listResultSources finds a run from its ORB or a folder above it", {
  root <- withr::local_tempdir()
  mp <- make_orb(root, user_bacs = "Haemophilus influenzae")
  run <- make_amrml_run(file.path(dirname(mp), "amRml", "2026-10-18_130824"))

  for (r in c(root, dirname(dirname(mp)), dirname(mp))) {
    sources <- listResultSources(r)
    expect_equal(normalizePath(unname(sources)), normalizePath(run))
    expect_equal(names(sources), "Haemophilus influenzae (amRml run 2026-11-09 12:48:10)")
  }
})

test_that("a results folder with no runs loads nothing, not the demo data", {
  root <- withr::local_tempdir()
  make_orb(root)
  sources <- listResultSources(root)
  expect_length(sources, 0)
  expect_equal(nrow(loadMLResults(root, unname(sources), verbose = FALSE)), 0)
  expect_equal(nrow(loadTopFeat(root, unname(sources), verbose = FALSE)), 0)
})

# ── packaged demo data ───────────────────────────────────────────────────────

test_that("the demo data is an ORB with an amRml run", {
  extdata <- system.file("extdata", package = "amRviz")
  expect_length(listResultSources(extdata), 1)

  perf <- loadMLResults(verbose = FALSE)
  expect_identical(unique(perf$species_label), "Shigella sonnei")
  expect_identical(unique(perf$species), "Sso")
  expect_true(any(perf$cross_test))
  expect_true(all(file.exists(unique(perf$metadata_path))))
  expect_false(any(perf$shuffled))
  expect_gt(nrow(loadTopFeat(verbose = FALSE)), 0)
})
