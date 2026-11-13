# Small ORBs and amRml runs on disk, for the ORB and results tests.

# Build {root}/{folder}/orb/ with a manifest and (optionally) its metadata file.
make_orb <- function(root, folder = "Haemophilus_influenzae", user_bacs = NULL,
                     metadata = TRUE, manifest_type = "amR_dataset",
                     schema_version = 1L, dataset_id = "Hin",
                     metadata_file = "metadata.parquet") {
  orb_dir <- file.path(root, folder, "orb")
  dir.create(orb_dir, recursive = TRUE)
  if (metadata) {
    arrow::write_parquet(
      tibble::tibble(genome.genome_id = "1.1", drug_class = "x"),
      file.path(orb_dir, metadata_file)
    )
  }
  manifest <- list(
    schema_version = schema_version,
    manifest_type = manifest_type,
    dataset_id = dataset_id,
    dataset = list(selection = list(user_bacs = as.list(user_bacs))),
    artifacts = list(amRml_input = list(metadata_parquet = metadata_file))
  )
  mp <- file.path(orb_dir, "manifest_run_20260101T000000000_pid1.json")
  jsonlite::write_json(manifest, mp, auto_unbox = TRUE, null = "null")
  mp
}

# Write an amRml run into `run_dir` (normally `{orb}/amRml/{run}/`): two
# matrices (a vanilla drug model and a cross-country model), each fitted with
# real and shuffled labels on one seed, in the layout amRml writes, with a
# manifest.json listing the tables. `manifest` entries are added to the
# manifest; `importance` replaces the importance rows.
make_amrml_run <- function(run_dir, manifest = list(), importance = NULL) {
  dir.create(file.path(run_dir, "results"), recursive = TRUE)
  dir.create(file.path(run_dir, "metadata"), recursive = TRUE)
  matrices <- tibble::tibble(
    matrix_id = c("m1", "m2"),
    mode_id = c("vanilla", "cross_country"),
    unit = c("drug", "drug_class"),
    target = c("AMP", "PEN"),
    train_group = c(NA, "USA"),
    test_group = c(NA, "GBR"),
    feature_type = c("Pfam", "gene"),
    encoding = c("binary", "counts")
  )
  tasks <- tibble::tibble(
    task_id = c("t1", "t2", "t3", "t4"),
    matrix_id = c("m1", "m1", "m2", "m2"),
    labels = c("real", "shuffled", "real", "shuffled"),
    seed = 31L
  )
  performance <- tibble::tibble(
    task_id = tasks$task_id, model = "LR", mcc = c(0.9, 0.1, 0.7, 0),
    auprc = 0.8, seed = 31L
  )
  importance <- importance %||% tibble::tibble(
    task_id = c("t1", "t3"), Variable = c("PF00001", "blaTEM"),
    Importance = c(0.5, 0.4), Sign = "POS"
  )
  arrow::write_parquet(matrices, file.path(run_dir, "metadata", "matrices.parquet"))
  arrow::write_parquet(tasks, file.path(run_dir, "results", "tasks.parquet"))
  arrow::write_parquet(performance, file.path(run_dir, "results", "performance.parquet"))
  arrow::write_parquet(importance, file.path(run_dir, "results", "importance.parquet"))
  run_manifest <- utils::modifyList(list(
    schema_version = 1L,
    manifest_type = "amR_ml",
    dataset_id = "Hin",
    dataset = list(label = "Haemophilus influenzae"),
    artifacts = list(amRviz_input = list(
      directory = ".",
      path_mode = "relative_to_manifest",
      metadata = list(matrices = "metadata/matrices.parquet"),
      results = list(
        tasks = "results/tasks.parquet",
        performance = "results/performance.parquet",
        importance = "results/importance.parquet"
      )
    )),
    runs = list(list(started_at = "2026-11-09 12:48:10.52974"))
  ), manifest)
  jsonlite::write_json(run_manifest, file.path(run_dir, "manifest.json"),
    auto_unbox = TRUE, null = "null"
  )
  run_dir
}
