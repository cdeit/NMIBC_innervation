# ==============================================================================
# 08_run_spatial_analysis.R
#
# Analyze archived spatial measurements without image reprocessing.
# ==============================================================================

source("config/config.R")
source("config/spatial_config.R")
# Shared output roots are defined once in the main config.
spatial_cfg$paths <- cfg$paths

# --- packages ----------------------------------------------------------------

library(dplyr)
library(glmmTMB)
library(emmeans)
library(ggplot2)

# --- functions ---------------------------------------------------------------

source("R/data_io.R")
source("R/read_table_s1.R")
source("R/spatial_analysis.R")

# --- data --------------------------------------------------------------------

require_input_files(c(spatial_cfg$patient_data_file, spatial_cfg$precomputed_file, cfg$paths$table_s4))
patient_df <- readRDS(spatial_cfg$patient_data_file)
spatial <- readRDS(spatial_cfg$precomputed_file)
# Publication image-to-patient mapping replaces image lists embedded in the RDS.
# Existing patient-level strata remain in the RDS until supplied in Table S1.
image_annotations <- read_table_s4_image_annotations(cfg$paths$table_s4)
image_df <- dplyr::transmute(
  image_annotations,
  im3_filename = `Sample Name`,
  include_image_final_spatial_cohort = include_final_spatial_cohort
)
patient_images <- image_annotations |>
  dplyr::filter(!is.na(StudyID)) |>
  dplyr::group_by(StudyID) |>
  dplyr::summarise(consolidated_im3s = paste(`Sample Name`, collapse = ", "), .groups = "drop")
patient_df <- patient_df |>
  dplyr::select(-consolidated_im3s) |>
  dplyr::left_join(patient_images, by = "StudyID")

# --- analysis ----------------------------------------------------------------

res_spatial_bcg <- run_spatial_tweedie(
  patient_df, image_df, spatial_cfg,
  outcome = "BCG_failure",
  band_data = spatial$band_enrichment
)
plot_spatial_tweedie(res_spatial_bcg, spatial_cfg, outcome = "BCG_failure")
res_spatial_prog <- run_spatial_tweedie(
  patient_df, image_df, spatial_cfg,
  outcome = "Progression",
  band_data = spatial$band_enrichment
)
plot_spatial_tweedie(res_spatial_prog, spatial_cfg, outcome = "Progression")

# --- object-distance distributions --------------------------------------------

master_nerve_vessel <- dplyr::bind_rows(
  spatial$nerve_dist_to_tumor, spatial$vessel_dist_to_tumor
)
# The distribution plots require the original nerve_color_map. Set its
# exact values in config/spatial_config.R; do not substitute another cluster palette.
if (is.null(spatial_cfg$distribution_palette)) {
  stop(
    "Spatial models completed, but distribution figures require the original ",
    "nerve_color_map at spatial_cfg$distribution_palette; see docs/REPRODUCIBILITY.md."
  )
}
plot_spatial_distributions(master_nerve_vessel, res_spatial_prog, spatial_cfg,
  nerve_color_map = spatial_cfg$distribution_palette
)
