# Build the packaged demo data (inst/extdata/<species>/) from an amRdata ORB
# and one amRml run.
#
# The demo is a small ORB laid out as amRviz reads it: the ORB's manifest and
# metadata at the top, and the run at amRml/<run>/ with its manifest.json.
# It keeps one seed; the baseline, stratified and cross-tested country and
# year models (real and shuffled labels); the top features of the real
# baseline and stratified models; and the gene and protein names those
# features use. Leave-one-out, cross-drug and MDR models, which no tab shows
# yet, are left out.
#
# Usage, from the package root:
#   Rscript data-raw/demo_orb.R path/to/orb path/to/run
# The source used for the packaged demo is Evan Brenner's S. sonnei ORB and
# its full amRml run (test-full), shared 2026-11-10.

args <- commandArgs(trailingOnly = TRUE)
orb <- args[1]
src_run <- args[2]
top_n <- 25
keep_modes <- c(
  "vanilla", "stratified_country", "stratified_year",
  "cross_country", "cross_year"
)
feature_modes <- c("vanilla", "stratified_country", "stratified_year")

orb_manifest_path <- list.files(orb, pattern = "^manifest_.*[.]json$", full.names = TRUE)[1]
orb_manifest <- jsonlite::read_json(orb_manifest_path)
run_manifest <- jsonlite::read_json(file.path(src_run, "manifest.json"))
input <- run_manifest$artifacts$amRviz_input
src_file <- function(section, entry) {
  file.path(src_run, input$directory %||% ".", input[[section]][[entry]])
}

out <- file.path(
  "inst", "extdata", gsub(" ", "_", run_manifest$dataset$label)
)
out_run <- file.path(out, "amRml", basename(src_run))
unlink(out, recursive = TRUE)
dir.create(file.path(out_run, "metadata"), recursive = TRUE)
dir.create(file.path(out_run, "results"), recursive = TRUE)

# The models kept: one seed, the modes and units the tabs show.
matrices <- arrow::read_parquet(src_file("metadata", "matrices")) |>
  dplyr::filter(
    .data$mode_id %in% keep_modes,
    .data$unit %in% c("drug", "drug_class")
  )
tasks <- arrow::read_parquet(src_file("results", "tasks"))
keep_seed <- min(tasks$seed)
tasks <- tasks |>
  dplyr::filter(.data$matrix_id %in% matrices$matrix_id, .data$seed == keep_seed)
performance <- arrow::read_parquet(src_file("results", "performance")) |>
  dplyr::filter(.data$task_id %in% tasks$task_id)

# Top features of the real baseline and stratified models.
feature_tasks <- tasks$task_id[
  tasks$labels == "real" &
    tasks$matrix_id %in% matrices$matrix_id[matrices$mode_id %in% feature_modes]
]
importance <- arrow::open_dataset(src_file("results", "importance")) |>
  dplyr::filter(.data$task_id %in% feature_tasks) |>
  dplyr::select("task_id", "Variable", "Importance", "Sign") |>
  dplyr::collect() |>
  dplyr::group_by(.data$task_id) |>
  dplyr::slice_max(abs(.data$Importance), n = top_n, with_ties = FALSE) |>
  dplyr::ungroup()

arrow::write_parquet(matrices, file.path(out_run, "metadata", "matrices.parquet"))
arrow::write_parquet(tasks, file.path(out_run, "results", "tasks.parquet"))
arrow::write_parquet(performance, file.path(out_run, "results", "performance.parquet"))
arrow::write_parquet(importance, file.path(out_run, "results", "importance.parquet"))

# The run's manifest, listing only the tables written here.
run_manifest$artifacts <- list(amRviz_input = c(
  input[c("status", "producer", "producer_run_id")],
  list(
    directory = ".",
    path_mode = "relative_to_manifest",
    metadata = list(matrices = "metadata/matrices.parquet"),
    results = list(
      tasks = "results/tasks.parquet",
      performance = "results/performance.parquet",
      importance = "results/importance.parquet"
    ),
    modes = as.list(keep_modes),
    seeds = list(keep_seed)
  )
))
run_manifest$dataset$input_orb <- list(directory = "../..")
jsonlite::write_json(run_manifest, file.path(out_run, "manifest.json"),
  auto_unbox = TRUE, null = "null", pretty = TRUE, digits = NA
)

# The ORB's manifest and metadata, as amRdata wrote them.
file.copy(orb_manifest_path, out)
file.copy(
  file.path(orb, orb_manifest$artifacts$amRml_input$metadata_parquet), out
)

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

message(
  "Wrote ", out, ": ", nrow(tasks), " models, ", nrow(performance),
  " performance rows, ", nrow(importance), " importance rows, ",
  nrow(genes), " gene names, ", nrow(proteins), " protein names"
)
