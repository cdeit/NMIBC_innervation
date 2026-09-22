# ==============================================================================
# Spatial preprocessing and analysis configuration (independent of cfg)
# Numerical settings and manuscript palettes are preserved.
# Paths are relative to the repository root; see docs/REVIEW_ITEMS.md before
# regenerating spatial data. Published image-QC flags are applied per row.
# ==============================================================================

spatial_cfg <- list(
  raw_cd31_dir = "data/external/02_segmentation_masks/vessel_masks",
  raw_nerve_dir = "data/external/02_segmentation_masks/tissue_and_nerve_masks",
  nerve_rdata_file = "data/processed/nerve_obj.rda",
  nerve_obj_name = "nerve_objs",
  images_include_name = "images_include",
  max_match_dist_px = 10,
  nerve_tumor_class = 0L,
  nerve_stroma_class = 1L,
  nerve_glass_class = 2L,
  min_area_px = 400,
  max_area_px = 1e+05,
  min_circularity = 0.2,
  min_lumen_frac = 0.035,
  min_lumen_frac_confident = 0.1,
  min_solidity_nolumen = 0.93,
  min_area_px_nolumen_round = 400,
  max_round_nolumen_branches = 100,
  min_elongated_area_px = 400,
  max_radius_cv_elongated = 0.6,
  min_solidity_elongated = 0.85,
  min_circularity_elongated = 0.1,
  min_eccentricity_elongated = 0.85,
  max_pa_ecc_ratio = 0.2,
  max_roundness_cshape = 0.4,
  max_radius_cv_cshape = 0.5,
  min_solidity_cshape = 0.87,
  min_area_cshape = 500,
  min_solidity_auto = 0.97,
  min_area_linear = 100,
  max_roundness_linear = 0.6,
  max_radius_cv_linear = 0.75,
  min_solidity_linear = 0.7,
  min_eccentricity_linear = 0.85,
  nerve_phenotype_lookup = structure(list(phenotype_code = c(1L, 2L, 4L, 5L, 6L, 7L), phenotype = c(
    "SP_pos_VGLUT1_pos_sc1",
    "TH_mid_sc2", "VACHT_pos_sc4", "TH_high_sc5", "SYP_high_sc6",
    "TH_neg_VACHT_neg_SP_neg_VGLUT1_neg_sc7"
  )), class = "data.frame", row.names = c(
    NA,
    -6L
  )),
  masks_dir = "data/external/03_processed_spatial_masks",
  clinical_rdata = "data/processed/spatial_clinical.RData",
  clinical_obj_name = "patient_nerve_quant_07302025",
  manual_eval_path = "data/processed/spatial_manual_eval.xlsx",
  spatial_manual_exclude_images = c(
    "120924_p9huP98_8_slide03_TMA_3_Cdj_Core1_13_G_16255_48862",
    "120924_p9huP98_8_slide04_TMA_5_Cdj_Core1_13_C_8007_47986", "120924_p9huP98_8_slide04_TMA_5_Cdj_Core1_6_I_18214_37024",
    "120924_p9huP98_8_slide04_TMA_5_Cdj_Core1_17_B_6210_54199"
  ),
  file_tumor = "tumor.tif",
  file_stroma = "stroma.tif",
  file_glass = "glass.tif",
  file_vessel = "vessel_filtered.tif",
  file_vessel_unfiltered = "vessel_unfiltered.tif",
  file_nerve_phenotyped = "nerve_phenotyped.tif",
  nerve_codes = c(1L, 2L, 4L, 5L, 6L, 7L),
  nerve_subtype_files = c(
    `1` = "nerve_SP_pos_VGLUT1_pos_sc1.tif", `2` = "nerve_TH_mid_sc2.tif",
    `4` = "nerve_VACHT_pos_sc4.tif", `5` = "nerve_TH_high_sc5.tif",
    `6` = "nerve_SYP_high_sc6.tif", `7` = "nerve_TH_neg_VACHT_neg_SP_neg_VGLUT1_neg_sc7.tif"
  ),
  nerve_labels = c(
    `1` = "SP+/VGLUT1+ sc1", `2` = "TH.mid sc2", `4` = "VACHT+ sc4",
    `5` = "TH.high sc5", `6` = "SYP.high sc6", `7` = "TH-/VACHT-/SP-/VGLUT1- sc7"
  ),
  um_per_px = 0.5,
  max_glass_contact_frac = 0.9,
  min_tissue_contact_frac = 0.1,
  min_tumor_area_px = 25000L,
  min_stroma_area_px = 25000L,
  min_nerve_object_area_px = 5L,
  min_vessel_object_area_px = 30L,
  border_dilate_px = 3L,
  distmap_downsample_factor = 2L,
  border_breaks = c(0, 25, 50, 100, 200, 300, Inf),
  border_labels = c(
    "0-25 um", "26-50 um", "51-100 um", "101-200 um", "200-300 um",
    ">300 um"
  ),
  col_stroma = "#01005C",
  col_tumor = "#1874CD",
  col_vessel = "#6A0DAD",
  nerve_cols = c(
    `1` = "#C4B8F2", `2` = "#FFD1DC", `4` = "#96E082", `5` = "#F2178E",
    `6` = "#FFC125", `7` = "#FE6100"
  ),
  qc_border_width_px = 9L,
  qc_band_ring_width_px = 7L,
  qc_band_ring_dash_period_px = 16L,
  qc_nerve_fill_dilate_px = 5L,
  processed_dir = "data/processed/spatial",
  path = ".",
  spatial_tweedie = list(
    max_radius_um = 300, band_levels = c(
      "0-25 um", "26-50 um",
      "51-100 um", "101-200 um", "201-300 um"
    ), min_outer_band_area_um2 = 10000,
    min_band_area_um2 = 5000, nerve_group_keep = "High", exclude_subclusters = "subcluster_1",
    min_radius_nerve_area_um2 = 0, p_adjust_scope = "subcluster",
    weight_mode = "both", weight_transform = "sqrt", nerve_weight_cap_quantile = 0.99,
    outcome_levels = list(BCG_failure = c("No", "Yes"), Progression = c(
      "No Progression",
      "Progression"
    )), dharma_n_sim = 500, dharma_seed = 222
  ),
  precomputed_file = "data/processed/spatial_measurements.rds",
  patient_data_file = "data/processed/patient_nerve_quant.rds"
)
