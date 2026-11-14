## Tests for the Network tab: makeDrugFeatureNetwork, makeFeatureEgoNetwork,
## and networkUI.

# ── makeDrugFeatureNetwork ──────────────────────────────────────────────────

test_that("makeDrugFeatureNetwork returns NULL for NULL input", {
  expect_null(makeDrugFeatureNetwork(NULL, bug = "Sfl"))
})

test_that("makeDrugFeatureNetwork returns NULL for non-data.frame input", {
  expect_null(makeDrugFeatureNetwork("not a data frame", bug = "Sfl"))
})

test_that("makeDrugFeatureNetwork returns NULL for zero-row data", {
  expect_null(makeDrugFeatureNetwork(tibble::tibble(), bug = "Sfl"))
})

test_that("makeDrugFeatureNetwork returns a forceNetwork for demo data", {
  skip_if_not_installed("networkD3")

  top_data <- suppressMessages(loadTopFeat(verbose = FALSE))
  skip_if(nrow(top_data) == 0, "No demo top-features data")

  bug <- unique(top_data$species)[1]
  result <- makeDrugFeatureNetwork(top_data, bug = bug, top_n = 5)
  expect_s3_class(result, "htmlwidget")
})

# ── networkTables ──────────────────────────────────────────────────────────

# Baseline models for two drugs, plus a stratified and a cross-tested row the
# network leaves out. One feature is named like a drug.
network_rows <- function() {
  tibble::tibble(
    species = "Sso",
    drug_or_class = c("AMP", "AMP", "AMP", "TET", "TET", "AMP", "AMP"),
    feature_type = c("domains", "domains", "domains", "genes", "genes", "genes", "genes"),
    Variable = c("DDE_Tnp_Tn3", "DDE_Tnp_1", "PF00001_IPR000001", "tetA", "AMP", "strat", "cross"),
    Importance = c(0.9, 0.8, 0.7, 0.6, 0.5, 2, 3),
    strat_label = c(NA, NA, NA, NA, NA, "country", NA),
    cross_test = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE)
  )
}

test_that("networkTables links drugs to their baseline features", {
  net <- networkTables(network_rows(), "Sso", top_n = 10)
  expect_named(net$nodes, c("id", "label", "group", "key", "info"))
  expect_named(net$links, c("from", "to", "importance", "sign"))
  # Stratified and cross-tested models are left out.
  expect_false(any(c("strat", "cross") %in% net$nodes$label))
  expect_true(all(c(net$links$from, net$links$to) %in% net$nodes$id))
})

test_that("networkTables keeps the top features of each type per drug", {
  rows <- network_rows()
  # Two more AMP genes, so AMP has three domains and three genes.
  rows <- dplyr::bind_rows(rows, tibble::tibble(
    species = "Sso", drug_or_class = "AMP", feature_type = "genes",
    Variable = c("g1", "g2"), Importance = c(0.4, 0.3),
    strat_label = NA, cross_test = FALSE
  ))
  net <- networkTables(rows, "Sso", top_n = 2)
  amp <- net$links[net$links$from == "drug:AMP", ]
  expect_equal(sort(sub(":.*", "", amp$to)), c("domains", "domains", "genes", "genes"))
  expect_setequal(
    net$nodes$label[net$nodes$id %in% amp$to & net$nodes$group == "genes"],
    c("g1", "g2")
  )
})

test_that("networkTables gives each link the sign of the feature's best model", {
  rows <- network_rows()[1:2, ]
  rows$Sign <- c("POS", "NEG")
  net <- networkTables(rows, "Sso")
  expect_equal(net$links$sign, c("POS", "NEG"))
  # Without a Sign column the links have none.
  expect_true(all(is.na(networkTables(network_rows(), "Sso")$links$sign)))
  # The sign comes from the model with the highest importance.
  two <- rbind(rows[1, ], rows[1, ])
  two$Importance <- c(0.2, 0.9)
  two$Sign <- c("POS", "NEG")
  expect_equal(networkTables(two, "Sso")$links$sign, "NEG")
})

test_that("makeDrugFeatureNetwork colours links by sign", {
  skip_if_not_installed("networkD3")
  rows <- network_rows()[1:2, ]
  rows$Sign <- c("POS", "NEG")
  net <- makeDrugFeatureNetwork(rows, "Sso")
  expect_equal(net$x$links$colour, c("#c0392b", "#2c7bb6"))
})

test_that("networkTables keeps current Pfam names whole", {
  net <- networkTables(network_rows(), "Sso", top_n = 10)
  domains <- net$nodes$label[net$nodes$group == "domains"]
  expect_setequal(domains, c("DDE_Tnp_Tn3", "DDE_Tnp_1", "PF00001"))
})

test_that("networkTables keeps a feature named like a drug as its own node", {
  net <- networkTables(network_rows(), "Sso", top_n = 10)
  expect_equal(sum(net$nodes$label == "AMP"), 2)
  expect_equal(nrow(networkTables(network_rows(), "Sso", top_n = 1)$links), 2)
  expect_null(networkTables(network_rows(), "Sfl"))
})

# An ORB with dyads: protein p1 is in dyads 1 and 2, gene g1 in dyad 1, and
# the domain "DDE_Tnp_Tn3" and the AMR hit "tetA-NCBIFAM" sit on p1 and p2.
make_dyad_orb <- function(dir) {
  arrow::write_parquet(
    tibble::tibble(
      dyad_id = 1:3, protein = c("fig|1.peg.1", "fig|1.peg.1", "fig|1.peg.2"),
      gene = c("g1", "f1...f2", NA), dyad_type = "paired"
    ),
    file.path(dir, "dyads.parquet")
  )
  hit <- function(protein, name) {
    tibble::tibble(protein = protein, query_name = name)
  }
  arrow::write_parquet(
    hit(c("fig|1.peg.1", "fig|1.peg.2"), "DDE_Tnp_Tn3"),
    file.path(dir, "protein_Pfam.parquet")
  )
  arrow::write_parquet(
    hit("fig|1.peg.2", "tetA-NCBIFAM"),
    file.path(dir, "protein_AMRFinder.parquet")
  )
  arrow::write_parquet(
    tibble::tibble(proteinID = "fig|1.peg.1", proteinName = "Protein <one>"),
    file.path(dir, "protein_names.parquet")
  )
  arrow::write_parquet(
    tibble::tibble(Gene = "g1", Annotation = "Gene one"),
    file.path(dir, "gene_names.parquet")
  )
  arrow::write_parquet(
    tibble::tibble(
      genome_drug.antibiotic = c("ampicillin", "tetracycline"),
      drug_abbr = c("AMP", "TET"),
      drug_class = c("penicillins", "tetracyclines"),
      class_abbr = c("PEN", "TET")
    ),
    file.path(dir, "metadata.parquet")
  )
  file.path(dir, "metadata.parquet")
}

dyad_rows <- function(metadata_path) {
  tibble::tibble(
    species = "Sso", drug_or_class = "AMP",
    feature_type = c("proteins", "genes", "domains", "args", "genes", "genes"),
    Variable = c("fig|1.peg.1", "g1", "DDE_Tnp_Tn3", "tetA.NCBIFAM", "none", "f1~~~f2"),
    Importance = c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4),
    strat_label = NA_character_, cross_test = FALSE,
    metadata_path = metadata_path
  )
}

test_that("networkTables links features to their dyads", {
  rows <- dyad_rows(make_dyad_orb(withr::local_tempdir()))
  net <- networkTables(rows, "Sso", top_n = 10, include_dyads = TRUE)
  dyad_links <- net$links[grepl("^dyad:", net$links$to), ]
  expect_setequal(net$nodes$label[net$nodes$group == "dyad"], paste("Dyad", 1:3))
  expect_true(all(c(net$links$from, net$links$to) %in% net$nodes$id))
  # The protein is in two dyads, the gene in one, the domain in two proteins'
  # dyads, and the AMR hit (named as an R column) in one; "none" has no dyad.
  expect_equal(as.list(table(dyad_links$from)), list(
    `args:tetA` = 1L, `domains:DDE_Tnp_Tn3` = 3L, `genes:f1~~~f2` = 1L,
    `genes:g1` = 1L, `proteins:fig|1.peg.1` = 2L
  ))
  # Without the option there are no dyads.
  expect_false("dyad" %in% networkTables(rows, "Sso")$nodes$group)
})

test_that("networkTables gives each node hover text with its ID and name", {
  rows <- dyad_rows(make_dyad_orb(withr::local_tempdir()))
  net <- networkTables(rows, "Sso", top_n = 10, include_dyads = TRUE)
  info <- stats::setNames(net$nodes$info, net$nodes$id)
  expect_length(net$nodes$info, nrow(net$nodes))
  expect_match(info[["drug:AMP"]], "<b>ID:</b> AMP")
  expect_match(info[["drug:AMP"]], "<b>Type:</b> Drug<")
  expect_match(info[["drug:AMP"]], "<b>Name:</b> ampicillin")
  # Names come from the ORB, escaped for HTML.
  expect_match(info[["proteins:fig|1.peg.1"]], "<b>Name:</b> Protein &lt;one&gt;")
  expect_match(info[["genes:g1"]], "<b>Name:</b> Gene one")
  expect_match(info[["domains:DDE_Tnp_Tn3"]], "<b>Type:</b> Domain")
  # A feature without a name has no name line.
  expect_no_match(info[["genes:none"]], "Name")
  # A dyad shows its protein and gene, each with an ID and a name when known.
  expect_match(info[["dyad:1"]], "<b>Protein ID:</b> fig|1.peg.1", fixed = TRUE)
  expect_match(info[["dyad:1"]], "<b>Protein name:</b> Protein &lt;one&gt;")
  expect_match(info[["dyad:1"]], "<b>Gene ID:</b> g1")
  expect_match(info[["dyad:1"]], "<b>Gene name:</b> Gene one")
  expect_match(info[["dyad:2"]], "<b>Dyad ID:</b> 2")
  # A fused gene is "f1~~~f2" in the models and "f1...f2" in dyads.parquet.
  expect_match(info[["dyad:2"]], "<b>Gene ID:</b> f1~~~f2")
  expect_no_match(info[["dyad:3"]], "Gene")
})

test_that("a drug and a class sharing an abbreviation are separate nodes", {
  rows <- dyad_rows(make_dyad_orb(withr::local_tempdir()))
  rows <- rows[c(1, 2, 1, 2), ]
  rows$drug_or_class <- "TET"
  rows$drug_label <- c("drug", "drug", "drug_class", "drug_class")
  rows$Variable <- c("fig|1.peg.1", "g1", "fig|1.peg.1", "f1~~~f2")
  net <- networkTables(rows, "Sso")
  drugs <- net$nodes[net$nodes$group == "drug", ]
  expect_setequal(drugs$id, c("drug:TET", "drug_class:TET"))
  expect_setequal(drugs$label, c("TET (drug)", "TET (class)"))
  # Each has its own features.
  expect_setequal(net$links$to[net$links$from == "drug:TET"], c("proteins:fig|1.peg.1", "genes:g1"))
  expect_setequal(
    net$links$to[net$links$from == "drug_class:TET"],
    c("proteins:fig|1.peg.1", "genes:f1~~~f2")
  )
  info <- stats::setNames(drugs$info, drugs$id)
  expect_match(info[["drug:TET"]], "<b>Name:</b> tetracycline$")
  expect_match(info[["drug_class:TET"]], "<b>Type:</b> Drug class")
  expect_match(info[["drug_class:TET"]], "<b>Name:</b> tetracyclines")
  # Abbreviations used by only one unit keep their plain label.
  rows$drug_or_class[3:4] <- "TETC"
  expect_setequal(
    networkTables(rows, "Sso")$nodes$label[networkTables(rows, "Sso")$nodes$group == "drug"],
    c("TET", "TETC")
  )
})

test_that("makeDrugFeatureNetwork adds a hover box with each node's text", {
  skip_if_not_installed("networkD3")
  rows <- dyad_rows(make_dyad_orb(withr::local_tempdir()))
  net <- makeDrugFeatureNetwork(rows, "Sso", include_dyads = TRUE)
  hook <- net$jsHooks$render[[2]]
  expect_match(hook$code, "network-tip")
  expect_length(hook$data, length(net$x$nodes$name))
  expect_match(hook$data[[1]], "ID")
})

test_that("networkTables leaves dyads out when the ORB has none", {
  dir <- withr::local_tempdir()
  net <- networkTables(
    dyad_rows(file.path(dir, "metadata.parquet")), "Sso",
    include_dyads = TRUE
  )
  expect_false("dyad" %in% net$nodes$group)
})

# ── makeFeatureEgoNetwork ───────────────────────────────────────────────────

test_that("makeFeatureEgoNetwork returns NULL for NULL input", {
  expect_null(makeFeatureEgoNetwork(NULL, variable = "x"))
})

test_that("makeFeatureEgoNetwork returns NULL for zero-row data", {
  expect_null(makeFeatureEgoNetwork(tibble::tibble(), variable = "x"))
})

test_that("makeFeatureEgoNetwork returns NULL when variable is NULL or empty", {
  df <- tibble::tibble(Variable = "x", cluster = "c1", COG = "COG0001")
  expect_null(makeFeatureEgoNetwork(df, variable = NULL))
  expect_null(makeFeatureEgoNetwork(df, variable = ""))
})

test_that("makeFeatureEgoNetwork returns NULL when variable not in table", {
  df <- tibble::tibble(Variable = "x", cluster = "c1", COG = "COG0001")
  expect_null(makeFeatureEgoNetwork(df, variable = "not-here"))
})

test_that("makeFeatureEgoNetwork returns NULL when no clusters or COGs", {
  df <- tibble::tibble(
    Variable = "x",
    cluster = NA_character_,
    COG = NA_character_
  )
  expect_null(makeFeatureEgoNetwork(df, variable = "x"))
})

test_that("makeFeatureEgoNetwork returns forceNetwork for valid annotated row", {
  skip_if_not_installed("networkD3")

  df <- tibble::tibble(
    Variable = "x",
    cluster = "fig|42897.100.peg.487",
    COG = "COG0001, COG0002"
  )
  result <- makeFeatureEgoNetwork(df, variable = "x")
  expect_s3_class(result, "htmlwidget")
})

# ── networkUI ───────────────────────────────────────────────────────────────

test_that("networkUI returns a shiny tabPanel", {
  ui <- networkUI("network")
  expect_s3_class(ui, "shiny.tag")
})

# ── shared font + node-shape styling ─────────────────────────────────────────

test_that("network visualisations carry the shared font + node-shape hook", {
  skip_if_not_installed("networkD3")

  top_data <- suppressMessages(loadTopFeat(verbose = FALSE))
  skip_if(nrow(top_data) == 0, "No demo top-features data")

  net <- makeDrugFeatureNetwork(
    top_data,
    bug = unique(top_data$species)[1], top_n = 5
  )
  # Font matches the plotly/ggplot panels.
  expect_equal(net$x$options$fontFamily, "Arial, sans-serif")
  # Per-group node shapes are applied via an onRender hook.
  expect_false(is.null(net$jsHooks$render))
  expect_match(net$jsHooks$render[[1]]$code, "d3.symbol")
  # Shapes are mirrored in both the nodes and the legend.
  expect_match(net$jsHooks$render[[1]]$code, "node-shape")
  expect_match(net$jsHooks$render[[1]]$code, "legend-shape")
})
