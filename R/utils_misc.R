# Misc cross-cutting helpers: species-code handling, choice derivation.


# Species code regex (note: Esp. includes escaped period)
SPECIES_PATTERN <- "(Efa|Sau|Kpn|Aba|Pae|Esp\\.?)"


#' Normalise species codes for comparison
#'
#' Makes "Esp." and "Esp" equivalent by removing a single trailing dot, so
#' species filters match regardless of the trailing-dot convention. NA-safe.
#'
#' @param x Character vector (or coercible) of species codes.
#' @return A character vector with any single trailing dot removed.
#' @keywords internal
#' @noRd
normalize_species <- function(x) {
  x_chr <- as.character(x)
  x_chr[is.na(x_chr)] <- NA_character_
  stringr::str_replace_all(x_chr, "\\.$", "")
}


#' Derive drug/class choices for the holdouts tab
#'
#' Restricts to stratified (country/year, non-baseline) models, optionally
#' filters to one species, and returns the sorted unique drug/class labels.
#'
#' @param perf_data Combined performance tibble from loadMLResults().
#' @param bug Optional 3-letter species code to filter to.
#' @return A sorted character vector of unique drug/class choices (empty when
#'   there are no matching rows).
#' @keywords internal
#' @noRd
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
