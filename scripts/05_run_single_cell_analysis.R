# ==============================================================================
# 05_run_single_cell_analysis.R
#
# Map bulk STaN-associated genes onto the Chen single-cell reference.
# ==============================================================================

source("config/config.R")

# --- packages ----------------------------------------------------------------

library(Seurat)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

# --- functions ---------------------------------------------------------------

source("R/data_io.R")
source("R/singlecell_preprocessing.R")
source("R/singlecell_heatmaps.R")

# --- optional preprocessing ---------------------------------------------------

if ("--preprocess" %in% commandArgs(trailingOnly = TRUE)) {
  # The original multiplet-rate table has unequal column lengths. Selecting a
  # replacement affects doublet filtering and must match the original analysis.
  stop(
    "Chen preprocessing requires resolution of the multiplet-rate table ",
    "and reference environment; see docs/REPRODUCIBILITY.md. ",
    "Use the published annotated Seurat object for cell-type mapping."
  )
}

# --- data and figures ---------------------------------------------------------

deg_file <- file.path(
  cfg$path, cfg$paths$results$bulk_rnaseq, "tables",
  paste0("DESeq2_Results_Cohorts_AB_adj_", paste(cfg$deseq2_covariates, collapse = "_"), ".csv")
)
require_input_files(c(cfg$paths$single_cell, deg_file, cfg$paths$single_cell_annotations))
single_cell_results <- run_single_cell_heatmaps(
  cfg$paths$single_cell, deg_file, cfg$paths$single_cell_annotations,
  file.path(cfg$path, cfg$paths$results$single_cell, "figures")
)
