# ==============================================================================
# 04_run_DESeq2_GSEA.R
#
# Bulk DESeq2 analysis, volcano plot, and configured GMT enrichment analyses.
# ==============================================================================

source("config/config.R")

# --- packages ----------------------------------------------------------------

library(DESeq2)
library(dplyr)
library(edgeR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOSemSim)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(ggpubr)
library(stringr)
library(ggbreak)
library(ggplotify)
library(patchwork)

# --- functions ---------------------------------------------------------------

source("R/data_io.R")
source("R/read_table_s1.R")
source("R/DESeq2_and_GSEA.R")

# --- data --------------------------------------------------------------------

require_input_files(c(cfg$paths$table_s1, cfg$paths$a_counts, cfg$paths$b_counts, cfg$paths$gene_annotations))
df <- read_table_s1(cfg$paths$table_s1)
# The existing DESeq2 helper assigns groups from within-cohort abundance quantiles.
# S1 contains continuous abundance, not stored STaN strata.
df$nerve_group <- factor(NA_character_, levels = c("Low", "High"))
# This comparison asks how expression differs by STaN abundance, not BCG outcome.
# Require sufficient stroma for imaging quantification, but retain patients regardless
# of BCG adequacy; RNA-seq matching and STaN stratification occur in the helper.
df <- dplyr::filter(
  df,
  sum_stroma_area_um2 >= if (is.null(cfg$min_stroma_area_um2)) 250000 else cfg$min_stroma_area_um2
)
a_counts <- readRDS(cfg$paths$a_counts)
b_counts <- readRDS(cfg$paths$b_counts)
gene_anns <- read_workspace_object(cfg$paths$gene_annotations, "gene_anns")

# --- analysis and figures -----------------------------------------------------

res_deseq2 <- run_combined_deseq2(df, a_counts, b_counts, gene_anns, cfg)
plot_deg_volcano(res_deseq2$res, cfg, save = TRUE, print = TRUE)
# Preserve the original invocation: gene annotations were not passed to GSEA.
gsea_results <- run_all_gsea(res_deseq2$res, cfg, print = FALSE)
