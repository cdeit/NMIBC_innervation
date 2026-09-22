# ==============================================================================
# 01_run_nerve_clustering.R
#
# Identify nerve candidates, apply published phenotypes, and aggregate patients.
# ==============================================================================

source("config/config.R")

# --- packages ----------------------------------------------------------------

library(dplyr)
library(tidyr)
library(rlang)
library(stringr)
library(readxl)
library(ggplot2)
library(patchwork)

# --- functions ---------------------------------------------------------------

source("R/data_io.R")
source("R/read_table_s4.R")
source("R/nerve_clustering.R")
source("R/nerve_aggregation.R")

# --- data --------------------------------------------------------------------

require_input_files(c(cfg$paths$table_s4, cfg$paths$clin_matched, cfg$paths$clin_compact))
s4 <- load_table_s4(cfg$paths$table_s4)
clin_matched <- read.delim(cfg$paths$clin_matched, check.names = FALSE, na.strings = c("", "NA"))
clin_compact <- read.delim(cfg$paths$clin_compact, check.names = FALSE, na.strings = c("", "NA"))

# --- analysis ----------------------------------------------------------------

# Stage 1 uses the broader reviewed object population. Final image eligibility
# is applied during aggregation, not before marker scaling or phenotyping.
clustering <- run_nerve_object_kmeans_pipeline(
  object_count_df = s4$object_count_df,
  tissue_seg_df = s4$tissue_seg_df,
  publication_reference = s4$publication_reference,
  nerve_subclusters_final = cfg$clustering$nerve_subclusters_final
)
patient_df <- run_patient_nerve_quant_pipeline(
  ocd = s4$nerve_object_df,
  tsd = s4$tissue_seg_df,
  tsds = s4$tissue_seg_summary,
  nerve_obj_table = clustering$nerve_obj_table,
  clin_df_matched = clin_matched,
  clin_df_compact = clin_compact,
  image_annotations = s4$image_annotations,
  regions_to_omit = clustering$regions_to_omit
)

# --- save --------------------------------------------------------------------

dir.create(cfg$paths$processed_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(patient_df, cfg$paths$patient_data)
saveRDS(clustering$nerve_obj_table, cfg$paths$nerve_objects)
saveRDS(clustering, cfg$paths$clustering_results)
object_cols <- intersect(
  c("object_tag", "Sample Name", "TMA", "subcluster_num", "phenotype", "cell_type"),
  names(clustering$nerve_obj_table)
)
write.table(clustering$nerve_obj_table[, object_cols],
  file.path(cfg$paths$processed_dir, "nerve_object_assignments.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE, na = ""
)
write.table(patient_df, file.path(cfg$paths$processed_dir, "patient_nerve_quantification.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE, na = ""
)

# --- figures -----------------------------------------------------------------

# Diagnostic embeddings and manuscript embeddings use distinct populations and
# parameters. The published UMAP requires the original image list explicitly.
make_figures <- !"--no-figures" %in% commandArgs(trailingOnly = TRUE)
if (make_figures) {
  nerve_objects <- dplyr::mutate(clustering$nerve_obj_table, subcluster = factor(subcluster_num))
  plots <- list(
    stage1_wcss_silhouette.png = plot_wcss_silhouette_diagnostics(
      clustering$objs4pca_og_clustered,
      subtitle = "All segmented objects (Stroma)"
    ),
    stage2_wcss_silhouette.png = plot_wcss_silhouette_diagnostics(
      clustering$nerve_obj_table,
      subtitle = "Nerve-candidate objects"
    )
  )
  stage1 <- plot_og_clustering_diagnostics(
    clustering$objs4pca_og_clustered, clustering$scaled_objs4pca_og_clustered,
    clustering$og_cluster_marker_scores
  )
  plots$stage1_pca_umap_enrichment.png <- stage1$pca_umap_enrichment
  plots$stage1_marker_umaps.png <- stage1$marker_umaps
  plots$stage2_marker_enrichment.png <- plot_subcluster_enrichment(nerve_objects, "subcluster")
  # Retain the general diagnostic UMAP and its historical output filename.
  scaled_nerve <- dplyr::filter(
    clustering$scaled_objs4pca_og_clustered,
    object_tag %in% nerve_objects$object_tag
  )
  plots$stage2_umap.png <- plot_subcluster_umaps(nerve_objects, scaled_nerve, "subcluster")$combined_plot
  dimensions <- list(c(10, 4.5), c(10, 4.5), c(14, 8), c(10, 10), c(5, 10), c(13, 7))
  dir.create(file.path(cfg$path, cfg$paths$results$nerve_clusters, "figures"), recursive = TRUE, showWarnings = FALSE)
  for (i in seq_along(plots)) {
    ggplot2::ggsave(file.path(cfg$path, cfg$paths$results$nerve_clusters, "figures", names(plots)[i]), plots[[i]],
      width = dimensions[[i]][1], height = dimensions[[i]][2], dpi = 300
    )
  }
  if ("--manuscript-umap" %in% commandArgs(trailingOnly = TRUE)) {
    require_input_files(cfg$paths$original_umap_images)
    original_images <- readRDS(cfg$paths$original_umap_images)
    final_images <- unique(trimws(unlist(strsplit(patient_df$consolidated_im3s, ",\\s*"))))
    manuscript_umap <- plot_subcluster_umap_figure(
      clustering$scaled_objs4pca_og_clustered, nerve_objects, "subcluster",
      final_im3s = final_images, og_im3s = original_images
    )
    # Retain the manuscript preparation separately; panel/export selection needs
    # confirmation against the reference figure (docs/REPRODUCIBILITY.md).
    saveRDS(manuscript_umap, file.path(cfg$path, cfg$paths$results$nerve_clusters, "figures", "manuscript_umap.rds"))
  }
}
