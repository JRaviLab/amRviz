# Build the packaged demo data (inst/extdata/Shigella_flexneri/) from an
# amRdata ORB with one amRml run inside it.
#
# The demo is a small copy of the ORB: its manifest and metadata, the run's
# tables, every performance row, the top features of each real-label model,
# and the gene and protein names those features use. The ORB's files sit
# directly in the species folder, so the run is at
# inst/extdata/Shigella_flexneri/amRml/<run>/.
#
# Usage, from the package root:
#   Rscript data-raw/demo_orb.R path/to/orb [run_name]
# The source used for the packaged demo is Evan Brenner's S. flexneri ORB
# (sfl_models.zip, amRml run 2026-10-18_130824).

args <- commandArgs(trailingOnly = TRUE)
orb <- args[1]
run_name <- if (length(args) > 1) args[2] else basename(
  list.dirs(file.path(orb, "amRml"), recursive = FALSE)[1]
)
top_n <- 25
out <- file.path("inst", "extdata", "Shigella_flexneri")

src_run <- file.path(orb, "amRml", run_name)
out_run <- file.path(out, "amRml", run_name)
unlink(out, recursive = TRUE)
dir.create(file.path(out_run, "results"), recursive = TRUE)

# Manifest and metadata, as amRdata wrote them.
manifest <- list.files(orb, pattern = "^manifest_.*[.]json$", full.names = TRUE)
file.copy(manifest, out)
file.copy(file.path(orb, "metadata.parquet"), out)

# The run's tables, and every performance row.
for (f in c("tasks.parquet", "matrices.parquet", "results/performance.parquet")) {
  file.copy(file.path(src_run, f), file.path(out_run, f))
}

# Top features of each real-label model.
tasks <- arrow::read_parquet(file.path(src_run, "tasks.parquet"))
real <- tasks$task_id[tasks$labels == "real"]
importance <- arrow::read_parquet(file.path(src_run, "results", "importance.parquet")) |>
  dplyr::filter(.data$task_id %in% real) |>
  dplyr::group_by(.data$task_id) |>
  dplyr::slice_max(abs(.data$Importance), n = top_n, with_ties = FALSE) |>
  dplyr::ungroup()
arrow::write_parquet(importance, file.path(out_run, "results", "importance.parquet"))

# Names for the features kept. Older amRml wrote protein IDs with "." in
# place of "|", so either form is matched.
features <- unique(importance$Variable)
genes <- arrow::read_parquet(file.path(orb, "gene_names.parquet")) |>
  dplyr::filter(.data$Gene %in% features)
proteins <- arrow::read_parquet(file.path(orb, "protein_names.parquet")) |>
  dplyr::filter(
    .data$proteinID %in% features |
      sub("|", ".", .data$proteinID, fixed = TRUE) %in% features
  )
arrow::write_parquet(genes, file.path(out, "gene_names.parquet"))
arrow::write_parquet(proteins, file.path(out, "protein_names.parquet"))

message("Wrote ", out, ": ", nrow(importance), " importance rows, ",
        nrow(genes), " gene names, ", nrow(proteins), " protein names")
