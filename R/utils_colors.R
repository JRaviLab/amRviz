# Colour palettes shared across the dashboard plots.

# Unified palette for molecular feature scales (gene, protein, domain, struct).
# Keyed to match the values in the `feature_type` column.
SCALE_COLORS <- c(
  domains  = "#66a61e",
  genes    = "#e6ab80",
  proteins = "#87ceeb",
  struct   = "#a52a2a"
)

# AMR phenotype palette (R/S/I plus full-word + lowercase variants so it works
# regardless of how the column is encoded). Susceptible is intentionally
# neutral grey so Resistant amber stands out as the signal of interest.
PHENOTYPE_COLORS <- c(
  R = "#d4872a", Resistant = "#d4872a", resistant = "#d4872a",
  S = "#8a8a8a", Susceptible = "#8a8a8a", susceptible = "#8a8a8a",
  I = "#e6ab80", Intermediate = "#e6ab80", intermediate = "#e6ab80"
)

# Categorical palette for metadata plots (hosts, isolation sources, etc.).
# Use `meta_palette(n)` at call sites rather than referencing META_COLORS
# directly, so plots handle species with more than length(META_COLORS) unique
# values without falling back to NA / gray70.
META_COLORS <- c(
  "#4a6b8a", "#87ceeb", "#e6ab80", "#8b6b7a",
  "#4e9a9a", "#c4a35a", "#9b7fba", "#5b8db8",
  "#d4735e", "#d4872a", "#b5b5b5"
)


# Return a vector of `n` muted categorical colors. For n <= length(META_COLORS)
# returns the first n curated colors verbatim; for n above that, interpolates
# via colorRampPalette() so callers get a full palette no matter how many
# unique values their data has.
meta_palette <- function(n = length(META_COLORS)) {
  if (n <= 0) {
    return(character(0))
  }
  if (n <= length(META_COLORS)) {
    return(META_COLORS[seq_len(n)])
  }
  grDevices::colorRampPalette(META_COLORS)(n)
}
