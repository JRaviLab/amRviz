# Global variables for NSE (Non-Standard Evaluation)
# These are column names used in dplyr/data.table operations
# Declaring them here silences R CMD check notes about "no visible binding"

utils::globalVariables(c(
  # Column names used across the package
  "ARG_name",
  "COG_name",
  "Gene",
  "Importance",
  "Variable",
  "accession",
  "bal_acc",
  "configurations",
  "count",
  "data_type",
  "drug.antibiotic_name",
  "drug_class",
  "drug_classes",
  "drug_count",
  "drug_level",
  "drug_or_class_name",
  "drug_or_drug_class",
  "feature_type",
  "genome.collection_year",
  "genome.host_common_name",
  "genome.isolation_country",
  "genome.isolation_source",
  "genome.isolation_source_",
  "genome_drug.antibiotic",
  "genome_drug.resistant_phenotype",
  "group",
  "max_imp",
  "mcc",
  "predicted_drug",
  "proteinID",
  "source_file",
  "species",
  "species_count",
  "tested_country",
  "tested_year",
  "trained_country",
  "trained_year",

  # data.table operator
  ":="
))
