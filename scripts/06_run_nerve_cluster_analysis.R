# ==============================================================================
# 06_run_nerve_cluster_analysis.R
#
# Cluster abundance, composition, and Cox proportional hazards analyses.
# ==============================================================================

source("config/config.R")

# --- packages ----------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(survival)

# --- functions ---------------------------------------------------------------

source("R/data_io.R")
source("R/nerve_cluster_analysis.R")
source("R/survival.R")

# --- data --------------------------------------------------------------------

require_input_files(cfg$paths$patient_data)
df <- readRDS(cfg$paths$patient_data)
# Cluster associations and survival describe outcomes in the adequate-BCG cohort.
# Apply eligibility here, after upstream nerve phenotyping, to preserve its population.
df <- dplyr::filter(df, Adequate_BCG == "Yes")

# --- analysis and figures -----------------------------------------------------

# Each function applies its additional analysis-specific eligibility criteria.
for (outcome in c("BCG_failure", "Progression")) {
  plot_cluster_abundance(df, cfg, outcome_col = outcome)
  plot_cluster_composition(df, cfg, outcome_col = outcome, relative = TRUE)
}
run_cluster_cox(df, cfg)
