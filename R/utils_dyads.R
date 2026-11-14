# Dyads: the protein/gene pairs an ORB records in `dyads.parquet`.
#
# Every feature type connects to dyads through the ORB's own files. Proteins
# and genes are listed in `dyads.parquet` itself; domains, COGs and AMRFinder
# hits belong to proteins (`protein_Pfam.parquet`, `protein_COG.parquet`,
# `protein_AMRFinder.parquet`), and so to the dyads those proteins are in.

# The protein-level annotation files, and the feature type each gives (the
# names the tabs use).
.DYAD_PROTEIN_FILES <- c(
  domains = "protein_Pfam.parquet",
  cogs = "protein_COG.parquet",
  args = "protein_AMRFinder.parquet"
)

.dyad_cache <- new.env(parent = emptyenv())


#' The dyads of one ORB
#'
#' @param orb_dir An ORB folder.
#' @return A tibble of `dyad_id`, `protein`, `gene` and `dyad_type`, or NULL
#'   when the ORB has no readable `dyads.parquet`. Cached per session.
#' @keywords internal
#' @noRd
.read_dyads <- function(orb_dir) {
  key <- paste0("dyads:", orb_dir)
  if (!is.null(.dyad_cache[[key]])) {
    return(.dyad_cache[[key]])
  }
  dyads <- .read_parquet_safe(file.path(orb_dir, "dyads.parquet"), verbose = FALSE)
  if (!all(c("dyad_id", "protein", "gene") %in% names(dyads))) {
    return(NULL)
  }
  if (!"dyad_type" %in% names(dyads)) {
    dyads$dyad_type <- NA_character_
  }
  dyads <- dyads[c("dyad_id", "protein", "gene", "dyad_type")]
  # A fused gene is "folP~~~chbB" in the importance and gene-name tables but
  # "folP...chbB" here; use the first so the three meet.
  dyads$gene <- gsub("...", "~~~", dyads$gene, fixed = TRUE)
  .dyad_cache[[key]] <- dyads
  dyads
}


#' Feature-to-dyad links for one ORB
#'
#' @param orb_dir An ORB folder.
#' @return A tibble of `feature_type`, `feature` (named as in the importance
#'   tables) and `dyad_id`, or NULL when the ORB has no readable
#'   `dyads.parquet`. Cached per session.
#' @keywords internal
#' @noRd
.read_dyad_links <- function(orb_dir) {
  key <- paste0("links:", orb_dir)
  if (!is.null(.dyad_cache[[key]])) {
    return(.dyad_cache[[key]])
  }
  dyads <- .read_dyads(orb_dir)
  if (is.null(dyads)) {
    return(NULL)
  }

  links <- list(
    proteins = tibble::tibble(
      feature_type = "proteins", feature = dyads$protein, dyad_id = dyads$dyad_id
    ),
    genes = tibble::tibble(
      feature_type = "genes", feature = dyads$gene, dyad_id = dyads$dyad_id
    )
  )
  for (type in names(.DYAD_PROTEIN_FILES)) {
    ann <- .read_parquet_safe(
      file.path(orb_dir, .DYAD_PROTEIN_FILES[[type]]),
      verbose = FALSE
    )
    if (!all(c("protein", "query_name") %in% names(ann))) next
    joined <- dplyr::inner_join(
      dplyr::distinct(ann, .data$protein, .data$query_name),
      dyads[c("protein", "dyad_id")],
      by = "protein", relationship = "many-to-many"
    )
    links[[type]] <- tibble::tibble(
      feature_type = type,
      # The importance tables name features as R column names, so hits such
      # as "blaTEM-NCBIFAM" appear there as "blaTEM.NCBIFAM".
      feature = make.names(joined$query_name),
      dyad_id = joined$dyad_id
    )
  }

  result <- dplyr::bind_rows(links) |>
    dplyr::filter(!is.na(.data$feature), !is.na(.data$dyad_id)) |>
    dplyr::distinct()
  .dyad_cache[[key]] <- result
  result
}


#' Feature-to-dyad links for several ORBs
#'
#' @param orb_dirs ORB folders.
#' @return A tibble as from `.read_dyad_links()`, or NULL when none of the
#'   ORBs has dyads.
#' @keywords internal
#' @noRd
load_dyad_links <- function(orb_dirs) {
  orb_dirs <- unique(orb_dirs[!is.na(orb_dirs) & nzchar(orb_dirs)])
  links <- dplyr::bind_rows(lapply(orb_dirs, .read_dyad_links))
  if (!nrow(links)) {
    return(NULL)
  }
  dplyr::distinct(links)
}


#' Descriptions of an ORB's domains, COGs and AMR hits
#'
#' @param orb_dir An ORB folder.
#' @param types Feature types to read, from `.DYAD_PROTEIN_FILES`.
#' @return A tibble of `feature_type`, `feature` (named as the network names
#'   features) and `name`; empty when the ORB has none.
#' @keywords internal
#' @noRd
.read_feature_descriptions <- function(orb_dir, types) {
  out <- lapply(intersect(types, names(.DYAD_PROTEIN_FILES)), function(type) {
    ann <- .read_parquet_safe(
      file.path(orb_dir, .DYAD_PROTEIN_FILES[[type]]),
      verbose = FALSE
    )
    if (!all(c("query_name", "description") %in% names(ann))) {
      return(NULL)
    }
    ann <- dplyr::distinct(ann, dplyr::across(dplyr::any_of(
      c("query_name", "query_accession", "description")
    )))
    keys <- ann$query_name
    # Older amRml named domains by their Pfam ID, without its version.
    if (type == "domains" && "query_accession" %in% names(ann)) {
      keys <- c(keys, sub("[.][0-9]+$", "", ann$query_accession))
      ann <- dplyr::bind_rows(ann, ann)
    }
    tibble::tibble(
      feature_type = type,
      feature = .network_feature_name(type, make.names(keys)),
      name = ann$description
    )
  })
  dplyr::bind_rows(out)
}


#' Full names of an ORB's drugs and drug classes, by abbreviation
#'
#' @param metadata_path An ORB's metadata parquet.
#' @return A tibble of `kind` ("drug" or "drug class"), `abbr` and `name`
#'   (underscores shown as spaces); empty when the metadata has no
#'   abbreviations. Cached per session.
#' @keywords internal
#' @noRd
.read_drug_names <- function(metadata_path) {
  key <- paste0("drugs:", metadata_path)
  if (!is.null(.dyad_cache[[key]])) {
    return(.dyad_cache[[key]])
  }
  cols <- c("genome_drug.antibiotic", "drug_abbr", "drug_class", "class_abbr")
  meta <- tryCatch(
    arrow::read_parquet(metadata_path, col_select = dplyr::any_of(cols)),
    error = function(e) tibble::tibble()
  )
  pairs <- function(kind, name, abbr) {
    if (!all(c(name, abbr) %in% names(meta))) {
      return(NULL)
    }
    tibble::tibble(
      kind = kind, abbr = meta[[abbr]], name = gsub("_", " ", meta[[name]], fixed = TRUE)
    )
  }
  result <- dplyr::bind_rows(
    tibble::tibble(kind = character(), abbr = character(), name = character()),
    pairs("drug", "genome_drug.antibiotic", "drug_abbr"),
    pairs("drug class", "drug_class", "class_abbr")
  ) |>
    dplyr::filter(!is.na(.data$abbr), !is.na(.data$name)) |>
    dplyr::distinct()
  .dyad_cache[[key]] <- result
  result
}


#' Names of some proteins and genes of an ORB
#'
#' @param orb_dir An ORB folder.
#' @param proteins Protein IDs to look up.
#' @param genes Gene IDs to look up.
#' @return A list of `proteins` and `genes`, each a named character vector of
#'   names (IDs without a name are left out).
#' @keywords internal
#' @noRd
.read_gene_protein_names <- function(orb_dir, proteins, genes) {
  proteins <- unique(proteins[!is.na(proteins)])
  path <- file.path(orb_dir, "protein_names.parquet")
  protein_names <- if (length(proteins) && file.exists(path)) {
    # The file holds millions of proteins; read only the ones asked for.
    tryCatch(
      arrow::open_dataset(path, format = "parquet") |>
        dplyr::filter(.data$proteinID %in% proteins) |>
        dplyr::collect(),
      error = function(e) NULL
    )
  }
  gene_names <- .read_parquet_safe(
    file.path(orb_dir, "gene_names.parquet"),
    verbose = FALSE
  )
  list(
    proteins = if (!is.null(protein_names) && nrow(protein_names)) {
      stats::setNames(protein_names$proteinName, protein_names$proteinID)
    } else {
      character(0)
    },
    genes = if (all(c("Gene", "Annotation") %in% names(gene_names))) {
      stats::setNames(gene_names$Annotation, gene_names$Gene)
    } else {
      character(0)
    }
  )
}


#' Hover text for the nodes of the network
#'
#' One block of HTML per node: its ID and name, and its type. A dyad also gets
#' its protein and gene, each with an ID and a name wherever the ORB has them.
#'
#' @param nodes Network nodes: `group` and `key` (the ID shown: a drug, a
#'   feature name, or a dyad ID).
#' @param orb_dirs The ORBs the models were fitted on.
#' @param metadata_paths Their metadata parquets, for the drug names.
#' @return A character vector with one HTML string per node.
#' @keywords internal
#' @noRd
network_node_info <- function(nodes, orb_dirs, metadata_paths = character(0)) {
  orb_dirs <- unique(orb_dirs[!is.na(orb_dirs) & nzchar(orb_dirs)])
  is_dyad <- nodes$group == "dyad"
  dyads <- dplyr::bind_rows(
    tibble::tibble(
      dyad_id = integer(), protein = character(), gene = character(),
      dyad_type = character()
    ),
    lapply(orb_dirs, .read_dyads)
  ) |>
    dplyr::distinct(.data$dyad_id, .keep_all = TRUE)
  dyad <- dyads[match(nodes$key[is_dyad], dyads$dyad_id), ]

  looked_up <- lapply(orb_dirs, .read_gene_protein_names,
    proteins = c(dyad$protein, nodes$key[nodes$group == "proteins"]),
    genes = c(dyad$gene, nodes$key[nodes$group == "genes"])
  )
  protein_names <- c(character(0), unlist(lapply(looked_up, `[[`, "proteins")))
  gene_names <- c(character(0), unlist(lapply(looked_up, `[[`, "genes")))
  described <- dplyr::bind_rows(
    tibble::tibble(
      feature_type = character(), feature = character(), name = character()
    ),
    lapply(orb_dirs, .read_feature_descriptions, types = unique(nodes$group))
  )

  # The name of each node: from its names table, by node type.
  name <- rep(NA_character_, nrow(nodes))
  name[nodes$group == "proteins"] <- protein_names[nodes$key[nodes$group == "proteins"]]
  name[nodes$group == "genes"] <- gene_names[nodes$key[nodes$group == "genes"]]
  other <- nodes$group %in% .DYAD_NAMED_TYPES
  name[other] <- described$name[match(
    paste(nodes$group[other], nodes$key[other]),
    paste(described$feature_type, described$feature)
  )]

  # A drug node is a drug or a drug class: its id starts with its unit.
  drug_names <- dplyr::bind_rows(
    tibble::tibble(kind = character(), abbr = character(), name = character()),
    lapply(metadata_paths[!is.na(metadata_paths)], .read_drug_names)
  ) |>
    dplyr::distinct()
  is_drug <- nodes$group == "drug"
  drug_kind <- ifelse(
    startsWith(nodes$id[is_drug], "drug_class:"), "drug class", "drug"
  )
  name[is_drug] <- vapply(seq_along(drug_kind), function(i) {
    hit <- drug_names$name[
      drug_names$abbr == nodes$key[is_drug][i] & drug_names$kind == drug_kind[i]
    ]
    if (length(hit)) paste(unique(hit), collapse = "; ") else NA_character_
  }, character(1))

  type <- unname(.scale_label_map()$labels[nodes$group])
  type[is_drug] <- ifelse(drug_kind == "drug class", "Drug class", "Drug")
  type[is_dyad] <- "Dyad"
  type[is.na(type)] <- nodes$group[is.na(type)]

  line <- function(label, value) {
    ok <- !is.na(value) & nzchar(value)
    ifelse(ok, paste0(
      "<b>", label, ":</b> ", htmltools::htmlEscape(ifelse(ok, value, "")), "<br>"
    ), "")
  }
  protein <- dyad$protein
  gene <- dyad$gene
  protein_name <- unname(protein_names[protein])
  gene_name <- unname(gene_names[gene])

  info <- paste0(
    line("Type", type),
    line("ID", nodes$key),
    line("Name", unname(name))
  )
  info[is_dyad] <- paste0(
    line("Type", paste("Dyad", dyad$dyad_type)),
    line("Dyad ID", nodes$key[is_dyad]),
    line("Protein ID", protein), line("Protein name", protein_name),
    line("Gene ID", gene), line("Gene name", gene_name)
  )
  sub("<br>$", "", info)
}

# Feature types whose names come from the ORB's protein annotation files.
.DYAD_NAMED_TYPES <- names(.DYAD_PROTEIN_FILES)
