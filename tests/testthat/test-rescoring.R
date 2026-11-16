## Tests for amRml's feature rescoring: reading it from a run (R/utils_rescoring.R)
## and drawing it in the network.

# Rescored features of two drugs. AMP has four features in its core group; DDE
# is found in both encodings. TET has one. "mdr" rows are left out.
rescoring_features <- function() {
  row <- function(drug, type, variable, subtype, score, contr, prior, cum,
                  sign = "NEG", label = "drug") {
    tibble::tibble(
      species = "Hin", drug_label = label, drug_or_class = drug,
      feature_type = type, feature_subtype = subtype, variable = variable,
      seed_ratio = 1, median_rank_score = score, median_contribution = contr,
      median_prior_cum_contrib = prior, median_cum_contrib = cum, sign = sign
    )
  }
  dplyr::bind_rows(
    row("AMP", "Pfam", "DDE", "binary", 1.00, 0.40, 0.00, 0.40),
    row("AMP", "Pfam", "DDE", "counts", 0.95, 0.30, 0.00, 0.30, "POS"),
    row("AMP", "gene", "g1", "binary", 0.90, 0.20, 0.40, 0.60),
    row("AMP", "gene", "g2", "binary", 0.80, 0.10, 0.60, 0.70),
    row("AMP", "gene", "g3", "counts", 0.70, 0.05, 0.70, 0.75),
    row("TET", "AMRFinder", "tetA", "binary", 1.00, 0.70, 0.00, 0.70, "POS", "drug_class"),
    row("MDR", "gene", "gx", "binary", 1.00, 1.00, 0.00, 1.00, "NEG", "mdr")
  )
}

# The model fits those features were scored from: the encodings of each group,
# with MCC and the shuffled-label MCC.
rescoring_models <- function() {
  fit <- function(drug, type, subtype, mcc, shuffled, label = "drug") {
    tibble::tibble(
      species = "Hin", drug_label = label, drug_or_class = drug,
      feature_type = type, feature_subtype = subtype, seed = 1, model = "LR",
      nonshuffled_MCC = mcc, shuffled_MCC = shuffled, MCC_diff = mcc - shuffled
    )
  }
  dplyr::bind_rows(
    fit("AMP", "Pfam", "binary", 0.90, 0.10),
    fit("AMP", "Pfam", "counts", 0.70, 0.10),
    fit("AMP", "gene", "binary", 0.80, 0.30),
    fit("AMP", "gene", "counts", 0.55, 0.05),
    fit("TET", "AMRFinder", "binary", 0.90, 0.00, "drug_class"),
    fit("MDR", "gene", "binary", 0.90, 0.00, "mdr")
  )
}

# Dyads reached by those features, as "protein|gene" with prefixed features.
rescoring_dyads <- function() {
  tibble::tibble(
    species = "Hin", drug_label = c("drug", "drug", "drug_class"),
    drug_or_class = c("AMP", "AMP", "TET"),
    dyad = c("fig|1.peg.1|g1", "fig|9.peg.9|g9...g8", "fig|2.peg.2|g4"),
    variables_csv = c("pfam:DDE,gene:g1,protein:fig|1.peg.1", "gene:g2", "amr:tetA"),
    dyad_median_rank_score = c(0.99, 0.80, 0.95),
    dyad_median_contribution = c(0.30, 0.10, 0.70),
    frequency = c(3, 1, 1)
  )
}

# A run with rescoring inside an ORB; returns the ORB's results root.
make_rescored_run <- function(root, ready = TRUE) {
  mp <- make_orb(root, folder = "Hin", user_bacs = "Haemophilus influenzae")
  run <- make_amrml_run(
    file.path(dirname(mp), "amRml", "run"),
    manifest = list(artifacts = list(feature_rescoring = list(
      status = if (ready) "ready" else "running",
      directory = "results/rescoring", path_mode = "relative_to_manifest"
    )))
  )
  dir.create(file.path(run, "results", "rescoring"))
  arrow::write_parquet(
    rescoring_features(), file.path(run, "results", "rescoring", "features_all_seeds.parquet")
  )
  arrow::write_parquet(
    rescoring_dyads(), file.path(run, "results", "rescoring", "dyads_top.parquet")
  )
  arrow::write_parquet(
    rescoring_models(), file.path(run, "results", "rescoring", "models_qualified.parquet")
  )
  list(root = root, run = run, orb = dirname(mp))
}

load_rescoring_features <- function() {
  .load_rescoring_source(make_rescored_run(withr::local_tempdir(), TRUE)$run,
    verbose = FALSE
  )
}

# ── reading ─────────────────────────────────────────────────────────────────

test_that("a run's rescoring is found through its manifest", {
  made <- make_rescored_run(withr::local_tempdir())
  expect_equal(
    .amrml_rescoring_dir(made$run),
    file.path(made$run, "results/rescoring")
  )
  # A run without the artifact, or one not marked ready, has none.
  expect_true(is.na(.amrml_rescoring_dir(make_rescored_run(withr::local_tempdir(), FALSE)$run)))
  bare <- make_amrml_run(file.path(withr::local_tempdir(), "amRml", "run"))
  expect_true(is.na(.amrml_rescoring_dir(bare)))
  expect_null(.load_rescoring_source(bare))
})

test_that("loadRescoring reads drugs and classes under the tabs' names", {
  made <- make_rescored_run(withr::local_tempdir())
  resc <- loadRescoring(made$root, made$run, verbose = FALSE)
  expect_named(resc, c("features", "dyads", "models"))
  # MDR fits are left out of the models too.
  expect_equal(nrow(resc$models), 5)
  expect_setequal(resc$models$feature_type, c("domains", "genes", "args"))
  # MDR targets are left out; feature types are named as the tabs name them.
  expect_equal(nrow(resc$features), 6)
  expect_setequal(resc$features$feature_type, c("domains", "genes", "args"))
  # The species and metadata come from the run's ORB, not the table.
  expect_equal(unique(resc$features$species), "Hin")
  expect_true(all(file.exists(resc$features$metadata_path)))
  expect_equal(nrow(resc$dyads), 3)
  # A selected run with no rescoring gives NULL, and nothing selected too.
  bare <- make_amrml_run(file.path(withr::local_tempdir(), "amRml", "run"))
  expect_null(loadRescoring(made$root, bare, verbose = FALSE))
  expect_null(loadRescoring(made$root, NULL, verbose = FALSE))
})

# ── thresholds ──────────────────────────────────────────────────────────────

amp <- function(x) x[x$drug_or_class == "AMP", ]

test_that("rescoredFeatures keeps a drug's features at its rank-score quantile", {
  resc <- load_rescoring_features()
  at <- function(...) {
    amp(rescoredFeatures(resc$features, "Hin", utils::modifyList(rescoringDefaults(), list(...))))
  }
  # The quantile is over every row of the drug, so 0 keeps all, and the
  # encodings of one feature collapse to one row.
  expect_setequal(at(quantile = 0)$Variable, c("DDE", "g1", "g2", "g3"))
  expect_equal(nrow(at(quantile = 0)), 4)
  expect_setequal(at(quantile = 0.75)$Variable, "DDE")
  # The best encoding gives the feature's scores and sign.
  best <- at(quantile = 0)
  expect_equal(best$rank_score[best$Variable == "DDE"], 1)
  expect_equal(best$sign[best$Variable == "DDE"], "NEG")
})

test_that("rescoredFeatures applies the core, cumulative and subtype thresholds", {
  resc <- load_rescoring_features()
  at <- function(...) {
    amp(rescoredFeatures(
      resc$features, "Hin",
      utils::modifyList(rescoringDefaults(), c(list(quantile = 0), list(...)))
    ))
  }
  # In the core while the share before it is below the threshold.
  expect_setequal(at(core = 0.5)$Variable, c("DDE", "g1"))
  expect_setequal(at(core = 0.3)$Variable, "DDE")
  # And its own cumulative share stays at or below the limit.
  expect_setequal(at(cumulative = 0.6)$Variable, c("DDE", "g1"))
  # Only DDE is found in both encodings.
  expect_setequal(at(both_subtypes = TRUE)$Variable, "DDE")
  expect_equal(nrow(rescoredFeatures(NULL, "Hin")), 0)
  expect_equal(nrow(rescoredFeatures(resc$features, "Sfl")), 0)
})

test_that("rescoredModels keeps the fits that pass the performance thresholds", {
  resc <- load_rescoring_features()
  at <- function(...) {
    rescoredModels(resc$models, "Hin", utils::modifyList(rescoringDefaults(), list(...)))
  }
  expect_equal(nrow(at()), 5)
  expect_equal(nrow(at(mcc = 0.85)), 2)
  # MCC above the shuffled fit: 0.8, 0.6, 0.5, 0.5, 0.9.
  expect_equal(nrow(at(delta_mcc = 0.55)), 3)
  expect_equal(nrow(rescoredModels(NULL, "Hin")), 0)
  expect_equal(nrow(rescoredModels(resc$models, "Sfl")), 0)
  scores <- rescoredDrugScores(at())
  amp <- scores[scores$drug_or_class == "AMP", ]
  expect_equal(c(amp$mcc, amp$delta_mcc, amp$n_fits), c(0.75, 0.55, 4))
  expect_equal(nrow(rescoredDrugScores(at(mcc = 1))), 0)
})

test_that("rescoredFeatures drops features whose fits no longer pass", {
  resc <- load_rescoring_features()
  th <- utils::modifyList(rescoringDefaults(), list(quantile = 0))
  at <- function(...) {
    th <- utils::modifyList(th, list(...))
    amp(rescoredFeatures(
      resc$features, "Hin", th, models = rescoredModels(resc$models, "Hin", th)
    ))
  }
  expect_setequal(at()$Variable, c("DDE", "g1", "g2", "g3"))
  # Only the Pfam binary fit reaches 0.85, so only its feature is left.
  expect_setequal(at(mcc = 0.85)$Variable, "DDE")
  # The genes' fits are only 0.5 above their shuffled fits, the Pfam ones 0.8
  # and 0.6.
  expect_setequal(at(delta_mcc = 0.6)$Variable, "DDE")
  expect_equal(nrow(at(mcc = 1)), 0)
})

test_that("rescoredDyads keeps the dyads that reach the features shown", {
  resc <- load_rescoring_features()
  shown <- rescoredFeatures(
    resc$features, "Hin", utils::modifyList(rescoringDefaults(), list(quantile = 0))
  )
  reached <- rescoredDyads(resc$dyads, shown)
  # The first dyad is reached by two shown features; "protein:fig|1.peg.1"
  # isn't a rescored feature of AMP, so it makes no link.
  expect_equal(sort(reached$Variable[reached$dyad == "fig|1.peg.1|g1"]), c("DDE", "g1"))
  expect_equal(reached$Variable[reached$dyad == "fig|2.peg.2|g4"], "tetA")
  expect_equal(reached$feature_type[reached$Variable == "tetA"], "args")
  # A dyad below the rank-score threshold is left out.
  strict <- rescoredDyads(
    resc$dyads, shown, utils::modifyList(rescoringDefaults(), list(dyad_rank_score = 0.9))
  )
  expect_false("fig|9.peg.9|g9...g8" %in% strict$dyad)
  expect_true("fig|9.peg.9|g9...g8" %in% reached$dyad)
  expect_equal(nrow(rescoredDyads(resc$dyads, shown[0, ])), 0)
})

# ── the network ─────────────────────────────────────────────────────────────

test_that("networkTablesRescored draws the top features of each type per drug", {
  resc <- load_rescoring_features()
  everything <- utils::modifyList(rescoringDefaults(), list(quantile = 0))
  net <- networkTablesRescored(resc, "Hin", top_n = 1, thresholds = everything)
  expect_named(net$nodes, c("id", "label", "group", "key", "info"))
  expect_true(all(c(net$links$from, net$links$to) %in% net$nodes$id))
  # One Pfam and one gene for AMP, by rank score; and TET's one ARG.
  expect_setequal(net$nodes$id[net$nodes$group != "drug"], c("domains:DDE", "genes:g1", "args:tetA"))
  expect_setequal(net$nodes$id[net$nodes$group == "drug"], c("drug:AMP", "drug_class:TET"))
  # Links carry the sign of the feature's best encoding.
  expect_equal(net$links$sign[net$links$to == "args:tetA"], "POS")
  expect_equal(net$links$sign[net$links$to == "domains:DDE"], "NEG")
  expect_null(networkTablesRescored(resc, "Sfl"))
  expect_null(networkTablesRescored(NULL, "Hin"))
})

test_that("networkTablesRescored adds the rescored dyads", {
  resc <- load_rescoring_features()
  everything <- utils::modifyList(rescoringDefaults(), list(quantile = 0))
  net <- networkTablesRescored(
    resc, "Hin",
    top_n = 10, include_dyads = TRUE, thresholds = everything
  )
  dyads <- net$nodes[net$nodes$group == "dyad", ]
  # Dyads that aren't in the ORB are named as rescored, with a fused gene
  # spelled as the ORB spells it.
  expect_setequal(dyads$key, c("fig|1.peg.1|g1", "fig|9.peg.9|g9~~~g8", "fig|2.peg.2|g4"))
  expect_true(all(c(net$links$from, net$links$to) %in% net$nodes$id))
  into <- net$links[net$links$to == "dyad:fig|1.peg.1|g1", ]
  expect_setequal(into$from, c("domains:DDE", "genes:g1"))
  # The hover box gives the protein and gene read from the dyad's name.
  info <- dyads$info[dyads$key == "fig|1.peg.1|g1"]
  expect_match(info, "<b>Protein ID:</b> fig|1.peg.1", fixed = TRUE)
  expect_match(info, "<b>Gene ID:</b> g1")
  expect_false("dyad" %in% networkTablesRescored(resc, "Hin", 10, FALSE, everything)$nodes$group)
})

test_that("the hover text of a rescored network gives the scores", {
  resc <- load_rescoring_features()
  everything <- utils::modifyList(rescoringDefaults(), list(quantile = 0))
  net <- networkTablesRescored(resc, "Hin", 10, TRUE, everything)
  info <- stats::setNames(net$nodes$info, net$nodes$id)
  # A drug: the median MCC and margin over shuffled of its passing fits.
  expect_match(info[["drug:AMP"]], "<b>MCC (median of passing fits):</b> 0.75", fixed = TRUE)
  expect_match(info[["drug:AMP"]], "<b>MCC above shuffled (median):</b> 0.55", fixed = TRUE)
  expect_match(info[["drug:AMP"]], "<b>Passing fits:</b> 4", fixed = TRUE)
  expect_match(info[["drug_class:TET"]], "<b>Passing fits:</b> 1", fixed = TRUE)
  # A feature: its scores for each drug it is shown for.
  expect_match(info[["domains:DDE"]], "<b>Scores by drug:</b><br>AMP: rank score 1.00, share 0.400, seed ratio 1.00", fixed = TRUE)
  # A dyad: its scores, and how many features reach it.
  expect_match(info[["dyad:fig|1.peg.1|g1"]], "AMP: rank score 0.99, share 0.300, 3 features", fixed = TRUE)
  expect_match(info[["dyad:fig|2.peg.2|g4"]], "TET (class): rank score 0.95, share 0.700, 1 feature", fixed = TRUE)
  # A dyad that isn't in the ORB has no type to show.
  expect_match(info[["dyad:fig|2.peg.2|g4"]], "<b>Type:</b> Dyad<br>", fixed = TRUE)
  # The baseline network has no scores.
  base <- networkTables(tibble::tibble(
    species = "Sso", drug_or_class = "AMP", feature_type = "genes", Variable = "g",
    Importance = 1, strat_label = NA, cross_test = FALSE
  ), "Sso")
  expect_false(any(grepl("Scores by drug|MCC", base$nodes$info)))
})

test_that("a long list of drugs in the hover text is cut short", {
  nodes <- tibble::tibble(id = "genes:g", group = "genes", key = "g")
  scores <- list(features = tibble::tibble(
    drug_label = "drug", drug_or_class = paste0("D", 1:9), feature_type = "genes",
    Variable = "g", rank_score = seq(1, 0.92, length.out = 9), contribution = 0.1,
    seed_ratio = 1
  ))
  lines <- .score_lines(nodes, scores, max_rows = 6)
  expect_match(lines, "D1: rank score 1.00", fixed = TRUE)
  expect_match(lines, "+ 3 more", fixed = TRUE)
  expect_no_match(lines, "D7")
})

test_that("makeDrugFeatureNetwork draws the rescoring when given it", {
  skip_if_not_installed("networkD3")
  resc <- load_rescoring_features()
  net <- makeDrugFeatureNetwork(NULL, "Hin", 2, TRUE, rescoring = resc)
  expect_s3_class(net, "htmlwidget")
  expect_match(net$jsHooks$render[[2]]$code, "network-tip")
  expect_null(makeDrugFeatureNetwork(NULL, "Sfl", 2, TRUE, rescoring = resc))
})

test_that("the Network tab draws rescored features only when asked and available", {
  made <- make_rescored_run(withr::local_tempdir())
  resc <- loadRescoring(made$root, made$run, verbose = FALSE)
  core <- list(
    available_species = function() "Hin",
    topFeatures = function() tibble::tibble(),
    rescoring = function() resc
  )
  shiny::testServer(networkServer, args = list(core = core), {
    session$setInputs(
      network_bug_id = "Hin", network_top_n = 2, network_include_dyads = FALSE,
      network_rescored = "no"
    )
    expect_false(use_rescoring())
    session$setInputs(
      network_rescored = "yes", network_resc_mcc = 0.5, network_resc_delta = 0,
      network_resc_quantile = 0, network_resc_core = 0.9,
      network_resc_cumulative = 0.75, network_resc_dyad = 0, network_resc_both = FALSE
    )
    expect_true(use_rescoring())
    expect_true(nzchar(output$drug_feature_network))
  })
  none <- utils::modifyList(core, list(rescoring = function() NULL))
  shiny::testServer(networkServer, args = list(core = none), {
    session$setInputs(network_rescored = "yes")
    expect_false(use_rescoring())
    expect_match(output$network_rescoring_note$html, "no rescoring")
  })
})

test_that("networkUI has the rescoring radio and its threshold sliders", {
  html <- as.character(networkUI("network"))
  expect_match(html, "network-network_rescored", fixed = TRUE)
  expect_match(html, "network-network_resc_quantile", fixed = TRUE)
  # Grouped by model performance, features and dyads; performance sliders grey
  # and dyad sliders green.
  expect_match(html, "network-network_resc_mcc", fixed = TRUE)
  expect_match(html, "network-sliders-grey", fixed = TRUE)
  expect_match(html, "network-sliders-green", fixed = TRUE)
  expect_match(html, "Model performance", fixed = TRUE)
  expect_match(html, .NETWORK_DYAD_COLOR, fixed = TRUE)
  expect_match(html, "input.network_rescored == &#39;yes&#39;", fixed = TRUE)
})
