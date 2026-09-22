# ==============================================================================
# Non-spatial manuscript analysis configuration
# Run scripts from the repository root. Paths may be edited for external inputs.
# Scientific settings are stage-specific: matching values do not imply a shared cohort.
# Figure-specific cosmetics are defined beside their plotting functions.
# ==============================================================================

cfg <- list()

# Recommended for reviewers: use the supplied precomputed spatial measurements
# to run downstream analyses without regenerating large image-derived datasets.
# Regenerating TIFFs and spatial measurements is computationally intensive,
# time-consuming, and requires substantial temporary/output storage.
# TRUE skips all image processing in workflow 07; workflow 08 reads the dataset
# configured in spatial_config.R. The precomputed data must be available locally.
# Set FALSE only to regenerate measurements from segmentation TIFFs, using the
# documented Alpine environment and required inputs (see docs/spatial/README.md).
cfg$use_precomputed_spatial <- TRUE

# --- path ---
cfg$path <- "."

# --- paths ---
cfg$paths <- list(
  table_s4 = "data/supp_tables/Table_S4_segmented_PGP95_SYP_objects.xlsx",
  table_s1 = "data/supp_tables/Table_S1_cohort_demographics.xlsx",
  clin_matched = "data/processed/clin_msi_matched.tsv", clin_compact = "data/processed/clin_compact.tsv",
  patient_data = "data/processed/patient_nerve_quant.rds",
  nerve_objects = "data/processed/nerve_objects_clustered.rds",
  clustering_results = "data/processed/nerve_clustering_results.rds",
  original_umap_images = "data/processed/original_umap_images.rds",
  publication_stage2_reference = "data/processed/publication_stage2_reference.rds",
  processed_dir = "data/processed",
  # Spatial measurement inputs/intermediates; statistical outputs use results$spatial.
  spatial = list(
    tissue_nerve_masks = "data/raw/tumor_stroma_nerve_masks",
    vessel_masks = "data/raw/vessel_masks",
    generated_masks = "data/processed/spatial_masks",
    measurement_runs = "data/processed/spatial_runs",
    measurements = "data/processed/spatial",
    image_manifest = "data/metadata/spatial/goal173.txt"
  ),
  results = list(
    clinical = "results/clinical",
    subtypes = "results/subtypes",
    bulk_rnaseq = "results/bulk_RNAseq",
    single_cell = "results/single_cell",
    nerve_clusters = "results/nerve_clusters",
    spatial = "results/spatial"
  ),
  a_counts = "data/external/BRS_cohort_A_counts.rds",
  b_counts = "data/external/BRS_cohort_B_counts.rds", gene_annotations = "data/external/BioMart_gene_anns.Rds",
  single_cell = "data/processed/scRNAseq_chen_annotated.rds",
  single_cell_10x = "data/external/PRJNA662018",
  single_cell_annotations = "data/metadata/nerve_communication_genes_ann.rds"
)

# --- clustering ---
cfg$clustering <- list(
  quant_cols = c(
    "Object PGP9.5 (Opal 480) Mean (Normalized Counts, Total Weighting).obj_count_data", "Object CD31 (Opal 520) Mean (Normalized Counts, Total Weighting).obj_count_data",
    "Object VACHT (Opal 540) Mean (Normalized Counts, Total Weighting).obj_count_data", "Object VGLUT1 (Opal 570) Mean (Normalized Counts, Total Weighting).obj_count_data",
    "Object TH (Opal 620) Mean (Normalized Counts, Total Weighting).obj_count_data", "Object SYP (Opal 650) Mean (Normalized Counts, Total Weighting).obj_count_data",
    "Object SP (Opal 690) Mean (Normalized Counts, Total Weighting).obj_count_data", "Object panCK (Opal 780) Mean (Normalized Counts, Total Weighting).obj_count_data",
    "Object DAPI Mean (Normalized Counts, Total Weighting).obj_count_data", "Object Autofluorescence Mean (Normalized Counts, Total Weighting).obj_count_data"
  ),
  seed = 222,
  k_og = 3,
  which_k_og = 3,
  default_region_filter_expr = quote(`Region Area (square microns).tissue_seg_data` > 2000),
  k_sub_min = 3,
  k_sub_max = 7,
  use_publication_cluster_labels = TRUE,
  omit_nonneuronal_cluster_3_from_publication = TRUE,
  nstart_denovo = 50,
  iter_max_denovo = 100,
  nerve_subclusters_final = c(1, 2, 4, 5, 6, 7)
)

# --- aggregation ---
cfg$aggregation <- list(
  min_obj_area_um2 = 0, max_obj_area_um2 = 3000, min_region_area_um2 = 2500,
  max_region_obj_density = 3000, min_stroma_area_percent = 20,
  min_tumor_area_percent = 5, min_combined_tissue_area_um2 = 1e+05,
  max_n_objs = Inf, max_obj_density_per_mm2 = 750, nerve_clusters_final = c(
    1,
    2, 4, 5, 6, 7
  ), patient_filters = list(
    min_sum_stroma_area_um2 = 250000,
    min_sum_tumor_area_um2 = NULL, min_sum_tissue_area_um2 = NULL,
    adequate_bcg = "Yes"
  )
)

# --- min stroma area um2 ---
cfg$min_stroma_area_um2 <- 250000

# --- subtype plots ---
cfg$subtype_plots <- list(
  nerve_col = "nerve_obj_area_percent", filter_col = "Adequate_BCG",
  filter_value = "Yes", followup_col = "Time_to_prog_FUend",
  min_per_level = 10, min_per_outcome = 3, subtype_cols = c(
    EAU_NMIBC = "EAU",
    EORTC_NMIBC = "EORTC", UROMOL_NMIBC = "Uromolclass", BRS_NMIBC = "BRS.merge",
    Chicago_T1_BCa = "TLumGU_sub", Lund_BCa = "Lund.subtype"
  ),
  outcomes = list(
    BCG_failure = c(No = "No", Yes = "Yes"),
    Progression = c(`No Progression` = "No", Progression = "Yes")
  )
)

# --- subtype km plots ---
cfg$subtype_km_plots <- list(
  nerve_col = "Pred", strata = "pred", min_per_level = 10,
  filters = NULL, subtype_cols = c(
    EORTC_NMIBC = "EORTC", UROMOL_NMIBC = "Uromolclass",
    BRS_NMIBC = "EMC.subtype", Chicago_T1_BCa = "TLumGU_sub",
    Lund_BCa = "Lund.subtype"
  ), outcomes = list(`HG-RFS probability` = list(
    time = "Time_to_BCG_fail_FUend", event = "BCG", event_map = c(
      Responder = 0,
      Failure = 1
    )
  ), `PFS probability` = list(
    time = "Time_to_prog_FUend",
    event = "Progression", event_map = c(
      `No Progression` = 0,
      Progression = 1
    )
  ))
)

# --- deseq2 condition ---
cfg$deseq2_condition <- "nerve_group"

# --- deseq2 ref ---
cfg$deseq2_ref <- "Low"

# --- deseq2 covariates ---
cfg$deseq2_covariates <- c("Sex", "Cohort")

# --- deseq2 cohorts ---
cfg$deseq2_cohorts <- "AB"

# --- quantile threshold ---
cfg$quantile_threshold <- 0.33

# --- min count ---
cfg$min_count <- 8

# --- min samples ---
cfg$min_samples <- 5

# --- gene type filter ---
cfg$gene_type_filter <- c("protein_coding", "lncRNA", "miRNA")

# --- deg padj cutoff ---
cfg$deg_padj_cutoff <- 0.05

# --- deg log2fc cutoff ---
cfg$deg_log2fc_cutoff <- 1

# --- gmt files ---
cfg$gmt_files <- list(
  Hallmarks = "h.all.v2025.1.Hs.symbols.gmt.txt", KEGG_2026 = "KEGG_2026.gmt.txt",
  PanglaoDB = "PanglaoDB_Augmented_2021.gmt.txt"
)

# --- gsea min size ---
cfg$gsea_min_size <- 15

# --- gsea max size ---
cfg$gsea_max_size <- 500

# --- gsea eps ---
cfg$gsea_eps <- 1e-10

# --- jaccard cutoff ---
cfg$jaccard_cutoff <- 0.25

# --- semantic cutoff ---
cfg$semantic_cutoff <- 0.25

# --- top gsets plot n ---
cfg$top_gsets_plot_n <- 30

# --- sig gsets only ---
cfg$sig_gsets_only <- TRUE

# --- gsea nes cutoff ---
cfg$gsea_nes_cutoff <- 1

# --- gsea padj cutoff ---
cfg$gsea_padj_cutoff <- 0.05

# --- nerve subclusters ---
cfg$nerve_subclusters <- list(
  marker_result = "n_clusters_7",
  marker_col = "marker_clean", delta_col = "delta", cluster_col = "cluster_when_n_clusters_7_and_OGnclust3",
  markers_keep = c("PGP9.5", "SYP", "SP", "TH", "VACHT", "VGLUT1"), sample_col = "Sample Name", object_col = "object_tag",
  im3_col = "consolidated_im3s", umap_seed = 444, umap_n_neighbors = 15,
  umap_min_dist = 0.0125, umap_repulsion = 2, umap_input_markers = c(
    "PGP9.5",
    "SYP", "TH", "VACHT", "SP", "VGLUT1", "CD31", "DAPI", "panCK",
    "AF"
  ), umap_plot_markers = c(
    "PGP9.5", "SYP", "TH", "VACHT",
    "SP", "VGLUT1"
  ), marker_cap_quantile = 0.98, stroma_area_col = "sum_stroma_area_um2",
  min_stroma_area = 250000, adequate_bcg_col = "Adequate_BCG",
  adequate_bcg_value = "Yes", total_nerve_col = "nerve_obj_area_percent",
  composition_min_total = 0.05, subcluster_cols = c(
    subcluster_1 = "nerve_obj_area_percent_subcluster_1",
    subcluster_2 = "nerve_obj_area_percent_subcluster_2", subcluster_4 = "nerve_obj_area_percent_subcluster_4",
    subcluster_5 = "nerve_obj_area_percent_subcluster_5", subcluster_6 = "nerve_obj_area_percent_subcluster_6",
    subcluster_7 = "nerve_obj_area_percent_subcluster_7"
  ), subcluster_levels = c(
    "subcluster_1",
    "subcluster_2", "subcluster_4", "subcluster_5", "subcluster_6",
    "subcluster_7"
  ), abundance_test = "wilcox", p_adjust = "BH",
  cox_predictors = c(
    `Cluster 2 (TH-mid)` = "nerve_obj_area_percent_subcluster_2",
    `Cluster 4 (VACHT+)` = "nerve_obj_area_percent_subcluster_4",
    `Cluster 5 (TH-high)` = "nerve_obj_area_percent_subcluster_5",
    `Cluster 6 (SYP-high)` = "nerve_obj_area_percent_subcluster_6",
    `Cluster 7 (Indeterminate)` = "nerve_obj_area_percent_subcluster_7",
    `Total STaN` = "nerve_obj_area_percent"
  ), cox_zero_threshold = 0.01,
  cox_adjust_covariates = "Age", cox_endpoints = list(
    `BCG Failure` = list(
      time = "Time_to_BCG_fail_FUend", event = "BCGfailR"
    ),
    Progression = list(time = "Time_to_prog_FUend", event = "ProgR")
  ),
  missing_as_level = c("Smoking", "CIS"), factor_refs = c(
    CIS = "NoCIS",
    LVI = "No", FGFR3_mut = "Yes", Tumor_focality = "Unifocal",
    Tumor_size = "Small"
  )
)

# --- heatmaps ---
cfg$heatmaps <- list(
  filter_col = "Adequate_BCG", filter_value = "Yes", stroma_area_col = "sum_stroma_area_um2",
  min_stroma_area = 250000, patient_id_col = "StudyID", nerve_col = "nerve_obj_area_percent",
  nerve_group_col = "nerve_group", composition_prefix = "proportion_obj_in_",
  composition_zero_threshold = 0.01, subtype_cols = c(
    `EAU NMIBC` = "EAU",
    `EORTC NMIBC` = "EORTC", `UROMOL NMIBC` = "Uromolclass",
    `BRS NMIBC` = "BRS.merge", `Chicago T1 BCa` = "TLumGU_sub",
    `Lund BCa` = "Lund.subtype"
  ), clin_cols = c(
    `BRS cohort` = "Cohort",
    Age = "Age", Sex = "Sex", Smoking = "Smoking", `reTURBT prior to BCG` = "Re_TUR",
    `Stage/grade` = "Stage_Grade", Substage = "Substage", CIS = "CIS",
    `Variant histology` = "Variant", `Tumor focality` = "Tumor_focality",
    `Tumor size` = "Size", LVI = "LVI", TILs = "TILs", `FGFR3 mutation` = "FGFR3_mut"
  ), surv_cols = c(
    BCG_failure = "BCG_failure", Criteria_failure = "Criteria_failure",
    Progression = "Progression"
  )
)

# --- color pals ---
cfg$color_pals <- list(nerve_subclusters = c(
  `1` = "#C4B8F2", `2` = "#F7C9D5",
  `3` = "#5D46F2", `4` = "#96E082", `5` = "#F2178E", `6` = "#FEC83E",
  `7` = "#F28700", `Cluster 1` = "#CABFFA", `Cluster 2` = "#FFD1DC",
  `Cluster 3` = "#654CFF", `Cluster 4` = "#B3FB8D", `Cluster 5` = "deeppink",
  `Cluster 6` = "goldenrod1", `Cluster 7` = "darkorange", subcluster_1 = "#CABFFA",
  `Cluster 1 (SP+ VGLUT1+)` = "#CABFFA", subcluster_2 = "#FFD1DC",
  `Cluster 2 (TH-mid)` = "#FFD1DC", subcluster_3 = "#654CFF", `Cluster 3` = "#654CFF",
  subcluster_4 = "#96E082", `Cluster 4 (VACHT+)` = "#96E082", subcluster_5 = "deeppink",
  `Cluster 5 (TH-high)` = "deeppink", subcluster_6 = "goldenrod1",
  `Cluster 6 (SYP-high)` = "goldenrod1", subcluster_7 = "darkorange",
  `Cluster 7 (Indeterminate)` = "darkorange", `Total STaN` = "black"
), nerve_strata = list(
  nerve_lo_hi = c("orange", "#E85D04"),
  stan_lo_hi = c(
    `STaN-low` = "orange", `STaN-high` = "#E85D04",
    Low = "orange", High = "#E85D04"
  ), median = c(
    `Below median` = "#FFB347",
    `Above median` = "#D94801"
  ), tertile = c(
    T1 = "#FFB347",
    T2 = "#D94801", T3 = "#8B0000", `1st (low)` = "#FFB347",
    `2nd` = "#D94801", `3rd (high)` = "#8B0000"
  ), pred = c(
    Low = "orange",
    High = "#E85D04", `Pred STaN-low` = "orange", `Pred STaN-high` = "#E85D04",
    `Predicted STaN-low` = "orange", `Predicted STaN-high` = "#E85D04"
  )
), subtypes = list(`BRS NMIBC` = c(
  BRS1 = "#009E73", BRS2 = "#0047AB",
  BRS3 = "#A4133C"
), `EAU NMIBC` = c(`high risk` = "#0077B6", `highest risk` = "#C23A2B"), `Chicago T1 BCa` = c(
  Early = "#3B6EA8", Inflam = "#E17A5F",
  LumGU = "#6A994E", MYC = "#A4133C", Tlum = "#7A5C9E"
), UROMOL.hc = c(
  `Class 1` = "#4E79A7",
  `Class 2` = "#C58A2E", `Class 3` = "#2A9D8F"
), `UROMOL NMIBC` = c(
  Class_1 = "#2F5597",
  Class_2a = "#5FA85F", Class_2b = "#F5B041", Class_3 = "#B23A48"
), `EORTC NMIBC` = c(`High risk` = "#00767b", `Highest risk` = "#F05A50"), `Lund BCa` = c(
  `UroA-Prog` = "#D55E00", `GU-Inf` = "#00A9A5",
  UroB = "#0072B2", GU = "#66A61E", `Mes-like` = "#ffaabb", `Ba/Sq-Inf` = "#332288",
  UroC = "#0B525B", `Uro-Inf` = "#56B4E9", `Sc/NE-like` = "#A4133C",
  `Ba/Sq` = "#882255"
), BRS_NMIBC = c(
  BRS1 = "#009E73", BRS2 = "#0047AB",
  BRS3 = "#A4133C"
), EAU_NMIBC = c(`high risk` = "#0077B6", `highest risk` = "#C23A2B"), Chicago_T1_BCa = c(
  Early = "#3B6EA8", Inflam = "#A06A2C",
  LumGU = "#6A994E", MYC = "#A4133C", Tlum = "#7A5C9E"
), UROMOL.hc = c(
  `Class 1` = "#4E79A7",
  `Class 2` = "#C58A2E", `Class 3` = "#2A9D8F"
), UROMOL_NMIBC = c(
  Class_1 = "#2F5597",
  Class_2a = "#5FA85F", Class_2b = "#F5B041", Class_3 = "#B23A48"
), EORTC_NMIBC = c(`High risk` = "#00767b", `Highest risk` = "#F05A50"), Lund_BCa = c(
  `Ba/Sq` = "#882255", `Ba/Sq-Inf` = "#332288",
  GU = "#66A61E", `GU-Inf` = "#00A9A5", `Mes-like` = "#ffaabb",
  `Sc/NE-like` = "#A4133C", `Uro-Inf` = "#56B4E9", `UroA-Prog` = "#D55E00",
  UroB = "#0072B2", UroC = "#0B525B"
)))
