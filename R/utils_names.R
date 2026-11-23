# Full names for the drug, drug class and country codes in the results.


#' Full names for the drug, drug class and country codes
#'
#' Read from the ORB metadata each row points to (`metadata_path`), which
#' pairs each code with its name: `drug_abbr` with `genome_drug.antibiotic`,
#' `class_abbr` with `drug_class`, and `country_abbr` with
#' `genome.isolation_country`. Underscores in names are shown as spaces.
#'
#' @param data Performance or top-features tibble with a `metadata_path`
#'   column.
#' @return A list of two named character vectors, code -> name: `drug` (drugs
#'   and drug classes) and `country`. Empty when there is no metadata.
#' @keywords internal
#' @noRd
amrNames <- function(data) {
  out <- list(drug = character(0), country = character(0))
  if (!is.data.frame(data) || !"metadata_path" %in% names(data)) {
    return(out)
  }
  paths <- unique(stats::na.omit(as.character(data$metadata_path)))
  paths <- paths[file.exists(paths)]
  if (!length(paths)) {
    return(out)
  }
  cols <- c(
    "genome_drug.antibiotic", "drug_abbr", "drug_class", "class_abbr",
    "genome.isolation_country", "country_abbr"
  )
  meta <- dplyr::bind_rows(lapply(paths, function(p) {
    tryCatch(
      arrow::read_parquet(p, col_select = dplyr::any_of(cols)),
      error = function(e) tibble::tibble()
    )
  }))
  pairs <- function(name, code) {
    if (!all(c(name, code) %in% names(meta))) {
      return(character(0))
    }
    d <- meta[!is.na(meta[[name]]) & !is.na(meta[[code]]), c(name, code)]
    d <- d[!duplicated(d[[code]]), ]
    stats::setNames(gsub("_", " ", d[[name]]), d[[code]])
  }
  list(
    drug = c(
      pairs("genome_drug.antibiotic", "drug_abbr"),
      pairs("drug_class", "class_abbr")
    ),
    country = pairs("genome.isolation_country", "country_abbr")
  )
}


#' Label codes with their full names, as "gentamicin (GEN)"
#'
#' @param code Codes to label.
#' @param lookup Named character vector, code -> name, from amrNames().
#' @return Labels; a code with no known name is returned as is.
#' @keywords internal
#' @noRd
.amr_label <- function(code, lookup) {
  code <- as.character(code)
  if (is.null(lookup)) {
    lookup <- character(0)
  }
  name <- unname(lookup[code])
  ifelse(is.na(name), code, paste0(name, " (", code, ")"))
}


#' Menu choices for codes, labeled with their full names
#'
#' @inheritParams .amr_label
#' @return `code`, named by its label, for `selectInput(choices = )`.
#' @keywords internal
#' @noRd
.amr_choices <- function(code, lookup) {
  stats::setNames(code, .amr_label(code, lookup))
}
