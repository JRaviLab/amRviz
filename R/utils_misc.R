# Misc cross-cutting helpers: species-code handling, choice derivation.


# Species code regex (note: Esp. includes escaped period)
SPECIES_PATTERN <- "(Efa|Sau|Kpn|Aba|Pae|Esp\\.?)"



# normalize_species helper: make "Esp." and "Esp" equivalent by removing
# a single trailing dot for comparisons (preserves NA).
normalize_species <- function(x) {
  x_chr <- as.character(x)
  x_chr[is.na(x_chr)] <- NA_character_
  stringr::str_replace_all(x_chr, "\\.$", "")
}


# getHoldoutsDrugChoices: derive drug/class choices for the holdouts tab
# perf_data: combined performance tibble from loadMLResults()
# bug: optional 3-letter species code to filter to
getHoldoutsDrugChoices <- function(perf_data, bug = NULL) {
  if (is.null(perf_data) || !is.data.frame(perf_data) || !nrow(perf_data)) {
    return(character(0))
  }

  # Filter to stratified models only (country or year, not baseline)
  df <- perf_data |>
    dplyr::filter(!is.na(.data$strat_label) & nzchar(.data$strat_label))

  if ("species" %in% names(df)) {
    df <- df |> dplyr::mutate(species = normalize_species(.data$species))
  }

  if (!is.null(bug) && "species" %in% names(df)) {
    bug_norm <- normalize_species(bug)
    df <- dplyr::filter(df, .data$species %in% bug_norm)
  }

  if (!nrow(df)) {
    return(character(0))
  }

  x <- as.character(df$drug_or_class)
  x[!is.na(x) & nzchar(trimws(x))] |>
    unique() |>
    sort()
}

