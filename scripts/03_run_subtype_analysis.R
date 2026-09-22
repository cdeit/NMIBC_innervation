# ==============================================================================
# 03_run_subtype_analysis.R
#
# Molecular subtype associations and survival using external gotNeRve predictions.
# ==============================================================================

source("config/config.R")

# --- packages ----------------------------------------------------------------

library(ComplexHeatmap)
library(cowplot)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(survival)

# --- functions ---------------------------------------------------------------

source("R/data_io.R")
source("R/survival.R")
source("R/clinical_plots.R")

# --- data and analyses --------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if ("--predictions" %in% args) {
  # Prediction implementation and all classifier development remain external.
  library(BRSpred)
  library(gotNeRve)
  unseen_pred <- gotNeRve::gotNeRve(gotNeRve::unseen_log2cpm, input = "log2cpm")
  unseen_pred_clin <- merge(unseen_pred, BRSpred::erasmus_clinical, by.x = "sample_id", by.y = 0)
  plot_km_subtype_overlay(unseen_pred_clin, cfg,
    filename = "KM_overlay_nerve_preds_x_molecular_subtype.pdf"
  )
  plot_km_subtype_faceted(unseen_pred_clin, cfg,
    filename = "KM_faceted_nerve_preds_x_molecular_subtype.pdf"
  )
} else {
  require_input_files(cfg$paths$patient_data)
  df <- readRDS(cfg$paths$patient_data)
  # Subtype associations use the same adequate-BCG cohort as clinical outcomes.
  df <- dplyr::filter(
    df, Adequate_BCG == "Yes",
    sum_stroma_area_um2 >= if (is.null(cfg$min_stroma_area_um2)) 250000 else cfg$min_stroma_area_um2
  )
  plot_stan_heatmaps(df, cfg, plots = "subtypes", complete_cases = TRUE)
  plot_nerves_by_subtype(df, cfg)
  if ("--imaging-survival" %in% args) {
    stop(
      "The active subtype KM configuration describes Pred/EMC.subtype/BCG. ",
      "Imaging-cohort endpoints and strata are not configured; ",
      "see docs/REPRODUCIBILITY.md. No alternative endpoints have been inferred."
    )
  }
}
