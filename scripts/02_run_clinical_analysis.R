# ==============================================================================
# 02_run_clinical_analysis.R
#
# STaN distributions, clinical associations, and clinical heatmaps.
# ==============================================================================

source("config/config.R")

# --- packages ----------------------------------------------------------------

library(ComplexHeatmap)
library(cowplot)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(patchwork)

# --- functions ---------------------------------------------------------------

source("R/data_io.R")
source("R/clinical_plots.R")

# --- data --------------------------------------------------------------------

require_input_files(cfg$paths$patient_data)
df <- readRDS(cfg$paths$patient_data)
# Clinical summaries use the adequate-BCG cohort defined for treatment-outcome analyses.
# Preserve the run-script >= boundary separately from function-level > filters.
df <- dplyr::filter(
  df, Adequate_BCG == "Yes",
  sum_stroma_area_um2 >= if (is.null(cfg$min_stroma_area_um2)) 250000 else cfg$min_stroma_area_um2
)

# --- figures -----------------------------------------------------------------

plot_tissue_histogram(df, xvar = "nerve_obj_area_percent", ylab = "Number of patients", cfg = cfg)
plot_stan_heatmaps(df, cfg, plots = "clinical", complete_cases = TRUE)
df_fig <- prep_patient_figure_data(df)
# Standardize Relapsing? to Relapsing; retain all other labels.
plot_stan_clinical_associations(df_fig, cfg)
plot_stan_heatmaps(df, cfg, plots = "full", complete_cases = FALSE)
