#' nerve_aggregation.R
#'
#' Object -> region -> tissue -> image -> patient filtering and aggregation
#' pipeline that turns a clustered nerve-object table (the output of
#' `nerve_clustering.R`) into per-patient nerve quantification
#' summaries.
#'
#' @section Stage 1:
#' Filters segmented objects down to the retained nerve subclusters,
#' applies object-area and region-area/density thresholds, and aggregates
#' object counts and areas up to the tissue-category level for each image.
#'
#' @section Stage 2:
#' Rolls per-image tissue summaries up to per-patient nerve quantification.
#' Images are restricted to those annotated `include_final_imaging_cohort`
#' and screened against minimum stroma/tumor area and maximum
#' object-density thresholds; images from the same tissue punch (multiple
#' scans of one physical core) are consolidated before aggregating to the
#' patient level.
#'
#' @section Cohort definition:
#' Which images and patients are included, and why others are excluded, is
#' recorded as data in the supplementary tables
#' (`include_final_imaging_cohort` and `reason_exclude_final_imaging_cohort`
#' in the image and patient tables). This script applies the image
#' annotations it is given, and then the patient-level cutoffs set in
#' `cfg$aggregation$patient_filters`, which remove patients that do not meet them. The
#' cutoffs only filter; no inclusion flags or reasons are generated.
#'
#' This script does not produce any figures; visualization of the
#' resulting per-patient quantification against clinical outcomes is
#' performed separately.
NULL

library(dplyr)
library(tidyr)
library(stringr)

## Parameters (object, region, image and patient-level cutoffs, including
## cfg$aggregation$patient_filters) are defined in config.R, which must be sourced first.


#' Stage 1: filter and aggregate object counts to the image level
#'
#' Filters objects to the retained nerve subclusters plus object/region
#' size and density thresholds, then aggregates object counts and areas
#' up to the tissue-category level for each image.
#'
#' @param ocd Unfiltered per-object marker/count table (object-count data).
#' @param tsd Per-region tissue segmentation table.
#' @param tsds Per-image tissue segmentation summary table.
#' @param nerve_obj_table Clustered per-object table (output of
#'   `nerve_clustering.R`'s `run_nerve_object_kmeans_pipeline()`),
#'   containing either a `subcluster_num` column or a `phenotype` column
#'   formatted as `"subcluster_N"`, and an `object_tag` column built with
#'   the same formula (from the same `ocd`, same row order).
#' @param tissue_category Tissue category to retain (default `"Stroma"`).
#' @param nerve_subclusters Numeric vector of subcluster indices to retain.
#'   Only used as a fallback when `nerve_obj_table` has no `cell_type`
#'   column -- when it does (the normal case, from
#'   `nerve_clustering.R`'s `run_nerve_object_kmeans_pipeline()`),
#'   `cell_type == "Nerve"` is used directly instead, since that's already
#'   the authoritative call. Defaults to [cfg$aggregation$nerve_clusters_final] (valid only
#'   for k = 7 tables lacking `cell_type`). An error is raised if any value
#'   here doesn't appear in `nerve_obj_table`'s subcluster assignments.
#' @param min_region_area_um2,max_region_obj_density Region-level area
#'   and density thresholds. Default [cfg$aggregation$min_region_area_um2]/[cfg$aggregation$max_region_obj_density].
#' @param min_obj_area_um2,max_obj_area_um2 Object-level area window.
#'   Default [cfg$aggregation$min_obj_area_um2]/[cfg$aggregation$max_obj_area_um2].
#' @param regions_to_omit Optional data frame with a `region_tag` column
#'   of regions to exclude (e.g. from `nerve_clustering.R`'s
#'   `compute_regions_to_omit()`). `NULL` to skip this exclusion.
#' @param include_obj_from_NAregions Whether to fold objects in
#'   unassigned regions back into `filtered_region_level_summary` (default
#'   `FALSE`).
#'
#' @return A list with:
#' \describe{
#'   \item{filtered_obj}{Object-level table after tagging and filtering.}
#'   \item{filtered_obj_summary_per_region}{Object counts/areas per region,
#'     overall and per subcluster.}
#'   \item{all_regions}{Full region list for `tissue_category`, from `tsd`.}
#'   \item{filtered_region_level_summary}{Region-level summary after the
#'     region-area/density/omission filters.}
#'   \item{unassigned_region_summary}{Summary of objects in unassigned
#'     regions (not folded in unless `include_obj_from_NAregions = TRUE`).}
#'   \item{filtered_tissue_summary}{Image-level (`Sample Name` x
#'     `Tissue Category`) summary -- the primary input to
#'     [summarize_summarized_msi_data_per_patient_v2()].}
#' }
filter_obj_region_tissue_level_data <- function(ocd,
                                                tsd,
                                                tsds,
                                                nerve_obj_table,
                                                tissue_category = "Stroma",
                                                nerve_subclusters = cfg$aggregation$nerve_clusters_final,
                                                min_region_area_um2 = cfg$aggregation$min_region_area_um2,
                                                max_region_obj_density = cfg$aggregation$max_region_obj_density,
                                                min_obj_area_um2 = cfg$aggregation$min_obj_area_um2,
                                                max_obj_area_um2 = cfg$aggregation$max_obj_area_um2,
                                                regions_to_omit = NULL,
                                                include_obj_from_NAregions = FALSE) {

  # Resolve a numeric subcluster identifier from whichever column shape
  # nerve_obj_table provides (nerve_clustering.R's `subcluster_num`, or
  # a `phenotype` column formatted as "subcluster_N"/"cluster_N"), used
  # below for the per-subcluster breakdown columns regardless of how the
  # Nerve/Stroma inclusion decision itself is made.
  if ("subcluster_num" %in% names(nerve_obj_table)) {
    nerve_obj_table$.subcluster_num <- nerve_obj_table$subcluster_num
  } else if ("phenotype" %in% names(nerve_obj_table)) {
    nerve_obj_table$.subcluster_num <- as.integer(gsub("^(subcluster|cluster)_", "", nerve_obj_table$phenotype))
  } else {
    stop("nerve_obj_table must contain either a `subcluster_num` column ",
         "or a `phenotype` column formatted as 'subcluster_N'.")
  }
  if (!"object_tag" %in% names(nerve_obj_table)) {
    stop("nerve_obj_table must contain an `object_tag` column built with the ",
         "same formula, from the same raw object-count data frame (same row ",
         "order), as `ocd` below.")
  }

  # Keep only objects retained as Nerve. Prefer nerve_obj_table's own
  # `cell_type` column when present -- it's already the authoritative
  # Nerve/Stroma call from nerve_clustering.R (by default anchored to
  # the published per-object assignment via object_tag, not to a
  # run-dependent numeric cluster index; see
  # run_nerve_object_kmeans_pipeline()'s `omit_nonneuronal_cluster_3_from_publication`).
  # Falls back to filtering by `nerve_subclusters` against `.subcluster_num`
  # only for tables that don't carry a `cell_type` column.
  if ("cell_type" %in% names(nerve_obj_table)) {
    nerve_objs_filtered <- nerve_obj_table %>%
      dplyr::filter(cell_type == "Nerve")
  } else {
    # Guard against a stale nerve_subclusters index set: subcluster identity
    # is specific to the k used in nerve_clustering.R, so a mismatch
    # here would otherwise silently filter to the wrong objects.
    observed_subclusters <- sort(unique(stats::na.omit(nerve_obj_table$.subcluster_num)))
    if (!all(nerve_subclusters %in% observed_subclusters)) {
      stop("nerve_subclusters (", paste(nerve_subclusters, collapse = ", "),
           ") includes values not present in nerve_obj_table's subcluster ",
           "assignments (observed: ", paste(observed_subclusters, collapse = ", "),
           "). Subcluster identity is specific to the k used in ",
           "nerve_clustering.R -- if you changed k there, update ",
           "nerve_subclusters/cfg$aggregation$nerve_clusters_final to match.")
    }
    nerve_objs_filtered <- nerve_obj_table %>%
      dplyr::filter(.subcluster_num %in% nerve_subclusters)
  }

  # Tag objects/regions on the raw object-count table, then restrict to the
  # target tissue category, retained nerve objects, and the object-area
  # window.
  # An `object_tag` column supplied in `ocd` (e.g. from Table S4) is used
  # as-is; otherwise it is built from the source table's row names.
  has_tag <- "object_tag" %in% names(ocd)

  filtered_obj <- ocd %>%
    dplyr::mutate(
      object_tag = if (has_tag) .data[["object_tag"]] else
        paste0(`Sample Name`, "_", `Tissue Category`, "_region",
               `Category Region ID.obj_count_data`, "_object",
               `Object ID.obj_count_data`, "_", row.names(.)),
      category_region_id_clean = dplyr::if_else(
        is.na(`Category Region ID.obj_count_data`), "unassigned",
        as.character(`Category Region ID.obj_count_data`)),
      region_tag = paste0(`Sample Name`, "_", `Tissue Category`, "_", category_region_id_clean)
    ) %>%
    dplyr::filter(
      `Tissue Category` == tissue_category,
      `object_tag` %in% nerve_objs_filtered$object_tag,
      `Object Area (square microns).obj_count_data` > min_obj_area_um2,
      `Object Area (square microns).obj_count_data` < max_obj_area_um2
    ) %>%
    dplyr::right_join(
      nerve_objs_filtered %>% dplyr::select(object_tag, .subcluster_num),
      by = "object_tag")

  # Total objects and area per region.
  filtered_obj_summary_per_region <- filtered_obj %>%
    dplyr::group_by(`Sample Name`, `Tissue Category`, category_region_id_clean, region_tag) %>%
    dplyr::summarize(
      region_filtered_n_total_obj = n(),
      region_filtered_obj_area_um2 = sum(`Object Area (square microns).obj_count_data`, na.rm = TRUE),
      .groups = "drop") %>%
    dplyr::ungroup()

  # Objects and area per subcluster per region, in wide format.
  subcluster_summary_per_region <- filtered_obj %>%
    dplyr::mutate(subcluster = .subcluster_num) %>%
    dplyr::group_by(`Sample Name`, `Tissue Category`, category_region_id_clean, region_tag, subcluster) %>%
    dplyr::summarize(
      n = dplyr::n(),
      obj_area_um2 = sum(`Object Area (square microns).obj_count_data`, na.rm = TRUE),
      .groups = "drop") %>%
    tidyr::pivot_wider(
      names_from  = subcluster,
      values_from = c(n, obj_area_um2),
      names_glue  = "region_filtered_{.value}_obj_in_subcluster_{subcluster}",
      values_fill = 0) %>%
    dplyr::ungroup()

  filtered_obj_summary_per_region <- filtered_obj_summary_per_region %>%
    dplyr::left_join(subcluster_summary_per_region,
                     by = c("Sample Name", "Tissue Category", "category_region_id_clean", "region_tag"))

  # Full list of valid regions, so regions with zero retained objects are
  # still represented (as zero counts) rather than dropped.
  all_regions <- tsd %>%
    dplyr::mutate(region_tag = paste0(`Sample Name`, "_", `Tissue Category`, "_",
                                       `Region ID.tissue_seg_data`),
                  category_region_id_clean = as.character(`Region ID.tissue_seg_data`)) %>%
    dplyr::select(`Sample Name`, `Tissue Category`, category_region_id_clean,
                  `Region Area (square microns).tissue_seg_data`, region_tag) %>%
    dplyr::filter(`Tissue Category` == tissue_category)

  # Join filtered counts onto the full region list, then apply region-level
  # area and density thresholds.
  filtered_region_level_summary <- all_regions %>%
    dplyr::left_join(filtered_obj_summary_per_region,
                     by = c("Sample Name", "Tissue Category", "category_region_id_clean", "region_tag")) %>%
    dplyr::mutate(
      dplyr::across(dplyr::starts_with("region_filtered_") & dplyr::where(is.numeric), ~ dplyr::coalesce(.x, 0)),
      region_filtered_obj_density_per_mm2 = dplyr::coalesce(
        ((region_filtered_n_total_obj / `Region Area (square microns).tissue_seg_data`) * 1e6), 0),
      region_filtered_obj_area_percent = dplyr::coalesce(
        (region_filtered_obj_area_um2 / `Region Area (square microns).tissue_seg_data`), 0)
    ) %>%
    dplyr::filter(
      `Region Area (square microns).tissue_seg_data` > min_region_area_um2,
      is.null(regions_to_omit) | !region_tag %in% regions_to_omit$region_tag,
      region_filtered_obj_density_per_mm2 < max_region_obj_density
    ) %>%
    dplyr::ungroup()

  # Objects in unassigned regions, summarized separately.
  unassigned_region_summary <- filtered_obj %>%
    dplyr::filter(category_region_id_clean == "unassigned") %>%
    dplyr::group_by(`Sample Name`, `Tissue Category`, category_region_id_clean, region_tag) %>%
    dplyr::summarize(
      region_filtered_n_total_obj = n(),
      region_filtered_obj_area_um2 = sum(`Object Area (square microns).obj_count_data`, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(`Region Area (square microns).tissue_seg_data` = NA,
                  region_filtered_obj_density_per_mm2 = NA,
                  region_filtered_obj_area_percent = NA) %>%
    dplyr::relocate(`Region Area (square microns).tissue_seg_data`, .before = region_tag) %>%
    as.data.frame()

  if (include_obj_from_NAregions) {
    filtered_region_level_summary <- dplyr::bind_rows(filtered_region_level_summary, unassigned_region_summary) %>%
      as.data.frame() %>%
      unique.data.frame()
  }

  # Aggregate region-level summaries up to tissue-category x image.
  filtered_tissue_summary <- filtered_region_level_summary %>%
    dplyr::group_by(`Sample Name`, `Tissue Category`) %>%
    dplyr::summarize(
      tissue_area_um2 = sum(`Region Area (square microns).tissue_seg_data`, na.rm = TRUE),
      dplyr::across(dplyr::starts_with("region_filtered_n_"), ~ sum(.x, na.rm = TRUE), .names = "{.col}"),
      dplyr::across(dplyr::starts_with("region_filtered_obj_area_um2"), ~ sum(.x, na.rm = TRUE), .names = "{.col}"),
      .groups = "drop") %>%
    dplyr::ungroup() %>%
    dplyr::rename_with(~ sub("^region_", "tissue_", .x),
                       c(dplyr::starts_with("region_filtered_n_"),
                         dplyr::starts_with("region_filtered_obj_area_um2"))) %>%
    dplyr::mutate(
      tissue_filtered_obj_density_per_mm2 = dplyr::if_else(
        tissue_area_um2 > 0, (tissue_filtered_n_total_obj / tissue_area_um2) * 1e6, 0))

  filtered_tissue_summary <- merge(
    filtered_tissue_summary,
    tsds[, c("Sample Name", "Tissue Category", "Tumor Area (um^2).tissue_seg_data_summary")],
    by = c("Sample Name", "Tissue Category"))

  list(
    filtered_obj = filtered_obj,
    filtered_obj_summary_per_region = filtered_obj_summary_per_region,
    all_regions = all_regions,
    filtered_region_level_summary = filtered_region_level_summary,
    unassigned_region_summary = unassigned_region_summary,
    filtered_tissue_summary = filtered_tissue_summary
  )
}


#' Stage 2: roll per-image summaries up to per-patient nerve quantification
#'
#' Restricts to images annotated for inclusion in the final imaging cohort,
#' applies image-level area/density/count thresholds, consolidates
#' duplicate scans of the same tissue punch (averaged), then aggregates to
#' the patient (`PA_number`) level.
#'
#' @param obj_df Image-level object summary, i.e.
#'   [filter_obj_region_tissue_level_data()]`$filtered_tissue_summary`.
#' @param tissue_seg_df Image-level tissue summary (e.g. tissue
#'   segmentation summary table).
#' @param clin_df_matched One row per image, with
#'   `im3_filename`/`PA_number`/`StudyID`/`MSI_tag`/`Cdj_tag`.
#' @param clin_df_compact One row per patient.
#' @param image_annotations Optional data frame with one row per image and
#'   columns `Sample Name` and `include_final_imaging_cohort` (logical);
#'   `reason_exclude_final_imaging_cohort` is carried for documentation
#'   (see the `image_annotations` sheet of Table S4). Only images with
#'   `include_final_imaging_cohort == TRUE` are retained. If `NULL`, all
#'   images are retained.
#' @param tissue_category Tissue category to retain (default `"Stroma"`).
#' @param min_stroma_area_pct,min_tumor_area_pct,min_combined_area_um2,max_object_density,max_n_obj
#'   Image-level QC thresholds. Default [cfg$aggregation$min_stroma_area_percent],
#'   [cfg$aggregation$min_tumor_area_percent], [cfg$aggregation$min_combined_tissue_area_um2],
#'   [cfg$aggregation$max_obj_density_per_mm2], [cfg$aggregation$max_n_objs].
#' @param additional_summary_cols,additional_summary Optional extra
#'   columns to retain (`additional_summary_cols`) and extra
#'   `dplyr::summarize()` expressions to compute per patient
#'   (`additional_summary`, as a named list spliced in via `!!!`).
#' @param cols_keep Optional explicit character vector of columns to
#'   retain from the merged image-level table; if `NULL`, a built-in
#'   default column set is used.
#'
#' @return A per-patient (`PA_number` x `Tissue Category`) data frame with
#'   `final_obj_area_percent`, `final_obj_density_per_mm2`,
#'   `sum_n_obj_in_subcluster_N`/`final_obj_area_percent_subcluster_N`,
#'   `proportion_obj_in_subcluster_N`, and merged clinical columns.
summarize_summarized_msi_data_per_patient_v2 <- function(
  obj_df,
  tissue_seg_df,
  clin_df_matched,
  clin_df_compact,
  image_annotations = NULL,
  tissue_category = "Stroma",

  min_stroma_area_pct = cfg$aggregation$min_stroma_area_percent,
  min_tumor_area_pct = cfg$aggregation$min_tumor_area_percent,
  min_combined_area_um2 = cfg$aggregation$min_combined_tissue_area_um2,
  max_object_density = cfg$aggregation$max_obj_density_per_mm2,
  max_n_obj = cfg$aggregation$max_n_objs,

  additional_summary_cols = "none",
  additional_summary = NULL,

  cols_keep = NULL
) {
  library(dplyr)
  library(stringr)

  # Merge object-level data with per-image clinical/ID annotations.
  dat2 <- merge(
    x = obj_df[which(obj_df$`Tissue Category` %in% tissue_category), ],
    y = (clin_df_matched %>% dplyr::select(im3_filename, PA_number, StudyID, MSI_tag, Cdj_tag)),
    by.x = "Sample Name", by.y = "im3_filename", all.x = TRUE
  ) %>% as.data.frame()

  # Merge with the tissue-segmentation summary (adds stroma/tumor area, etc.)
  dat2 <- merge(x = dat2, y = tissue_seg_df,
               by = c("Sample Name", "Tissue Category"), all.x = TRUE) %>% as.data.frame()

  # Restrict to images annotated for inclusion in the final imaging cohort.
  if (!is.null(image_annotations)) {
    stopifnot(all(c("Sample Name", "include_final_imaging_cohort") %in% names(image_annotations)),
              !anyDuplicated(image_annotations$`Sample Name`))
    imgs <- image_annotations$`Sample Name`[image_annotations$include_final_imaging_cohort %in% TRUE]
    dat2 <- dat2[which(dat2$`Sample Name` %in% imgs), ]
  }
  message("Number of image-level rows after cohort filtering: ", nrow(dat2))

  dat2 <- dat2 %>% dplyr::ungroup() %>% as.data.frame()
  names(dat2) <- gsub("percent", "%", names(dat2))
  names(dat2) <- gsub("per square mm", "per mm^2", names(dat2))
  names(dat2) <- gsub("square microns", "um^2", names(dat2))
  names(dat2) <- gsub("\\.y", "", names(dat2))

  if (is.null(cols_keep)) {
    cols_keep <- c("Sample Name", "StudyID", "Cdj_tag", "MSI_tag", "PA_number", "punch_id",
                   "BCG_failure", "Sex", "Progression.x", "Tumor_focality", "Location",
                   "Criteria_failure", "LVI", "BRS.merge", "Immune_infiltration_score",
                   "Tissue Category",
                   grep("Area", names(dat2), value = TRUE), "Tumor Area (um^2)", "Stroma Area (um^2)",
                   "Tumor+Stroma Area (um^2)", "Stroma Area (%)", "Tumor Area (%)", "Tumor:Stroma Ratio",
                   "Total Objects", "Object Density (per mm^2)", "Object Area (um^2)",
                   "Object Area (%)",
                   names(obj_df),
                   "Smoking", "consensusClass", "LumP_score", "LumU", "Ba.Sq", "LumNS",
                   "Stroma.rich", "NE.like", "Lund.subtype", "TCGA.subtype", "Uromolclass",
                   additional_summary_cols)
  }

  dat2 <- dat2 %>% dplyr::select(dplyr::any_of(cols_keep))

  # Apply image-level QC thresholds, consolidate duplicate scans of the
  # same tissue punch, then aggregate to the patient (PA_number) level.
  dat_compressed2 <- dat2 %>%
    dplyr::filter(
      `Tissue Category` == tissue_category,
      `Stroma Area (%).tissue_seg_data_summary` > min_stroma_area_pct,
      `Tumor Area (%).tissue_seg_data_summary` > min_tumor_area_pct,
      `Tumor+Stroma Area (um^2).tissue_seg_data_summary` > min_combined_area_um2,
      tissue_filtered_n_total_obj < max_n_obj,
      tissue_filtered_obj_density_per_mm2 < max_object_density
    ) %>%
    dplyr::mutate(core_id = paste0(stringr::str_split_i(pattern = "_\\[", string = `Sample Name`, i = 1))) %>%
    dplyr::group_by(core_id, `Tissue Category`) %>%
    dplyr::summarize(
      `Sample Name` = dplyr::first(`Sample Name`),
      tissue_filtered_obj_density_per_mm2 = mean(tissue_filtered_obj_density_per_mm2),
      tissue_filtered_n_total_obj = mean(tissue_filtered_n_total_obj),
      tissue_filtered_obj_area_um2 = mean(tissue_filtered_obj_area_um2),
      tissue_area_um2 = mean(tissue_area_um2),
      tumor_area_um2 = mean(`Tumor Area (um^2).tissue_seg_data_summary`, na.rm = TRUE),
      MSI_tag = dplyr::first(MSI_tag),
      Cdj_tag = dplyr::first(Cdj_tag),
      PA_number = dplyr::first(PA_number),
      StudyID = dplyr::first(StudyID),
      dplyr::across(dplyr::contains("subcluster"), ~ dplyr::first(.x)),
      .groups = "drop") %>%
    dplyr::ungroup() %>%
    dplyr::group_by(`Sample Name`, `Tissue Category`) %>%
    dplyr::sample_n(size = 1) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(PA_number, `Tissue Category`) %>%
    dplyr::summarize(
      mean_obj_density_per_mm2 = mean(tissue_filtered_obj_density_per_mm2, na.rm = TRUE),
      mean_n_obj = mean(tissue_filtered_n_total_obj, na.rm = TRUE),
      sum_filtered_obj_area_um2 = sum(tissue_filtered_obj_area_um2, na.rm = TRUE),
      sum_n_obj = sum(tissue_filtered_n_total_obj, na.rm = TRUE),
      sum_tissue_area_um2 = sum(tissue_area_um2, na.rm = TRUE),
      sum_tumor_area_um2 = sum(tumor_area_um2, na.rm = TRUE),
      n_consolidated = n(),
      dplyr::across(dplyr::starts_with("tissue_filtered_n_obj"), ~ sum(.x, na.rm = TRUE), .names = "{.col}"),
      dplyr::across(dplyr::starts_with("tissue_filtered_obj_area_um2_"), ~ sum(.x, na.rm = TRUE), .names = "{.col}"),
      StudyID = dplyr::first(StudyID),
      `Tissue Category` = dplyr::first(`Tissue Category`),
      consolidated_im3s = stringr::str_c(`Sample Name`, collapse = ", "),
      consolidated_MSI_tags = stringr::str_c(MSI_tag, collapse = ", "),
      consolidated_Cdj_tags = stringr::str_c(Cdj_tag, collapse = ", "),
      !!!additional_summary,
      .groups = "drop") %>%
    dplyr::rename_with(~ sub("^tissue_filtered", "sum", .x),
                       c(dplyr::starts_with("tissue_filtered_n_"),
                         dplyr::starts_with("tissue_filtered_obj_area_um2_"))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      final_obj_density_per_mm2 = (sum_n_obj / sum_tissue_area_um2) * 1e6,
      final_obj_area_percent = (sum_filtered_obj_area_um2 / sum_tissue_area_um2) * 100,
      sum_stroma_area_um2 = sum_tissue_area_um2,
      sum_total_area_um2 = sum_tissue_area_um2 + sum_tumor_area_um2,
      stroma_percent_of_total_area = (sum_stroma_area_um2 / sum_total_area_um2) * 100,
      tumor_percent_of_total_area = (sum_tumor_area_um2 / sum_total_area_um2) * 100,
      dplyr::across(dplyr::starts_with("sum_obj_area_um2_obj_in"), ~ (.x / sum_stroma_area_um2) * 100,
                    .names = "{sub('sum_obj_area_um2_obj_in_','final_obj_area_percent_',col)}"),
      dplyr::across(dplyr::starts_with("sum_n_obj_in_"), ~ (.x / sum_n_obj),
                    .names = "{sub('sum_n_obj_in_', 'proportion_obj_in_', col)}")
    ) %>%
    unique.data.frame()

  # Merge with compact (one-row-per-patient) clinical data.
  dat_compressed2 <- merge(x = clin_df_compact, y = dat_compressed2,
                           by = c("PA_number", "StudyID"), all.y = TRUE) %>%
    dplyr::ungroup() %>%
    unique.data.frame() %>%
    dplyr::filter(!is.na(PA_number)) %>%
    # clin_df_compact is not always strictly one row per patient -- a patient
    # whose tissue spans multiple TMAs has one clin_df_compact row per TMA
    # (differing only in TMA/TMA_number). unique.data.frame() above only
    # drops byte-identical rows, so it doesn't catch these; the nerve-quant
    # aggregate itself is identical regardless of which TMA row it merges
    # onto, so collapse to one row per (PA_number, StudyID) here.
    dplyr::distinct(PA_number, StudyID, .keep_all = TRUE)

  dat_compressed2
}


#' Filter patients by the cutoffs in `cfg$aggregation$patient_filters`
#'
#' Removes patients that do not meet the configured cutoffs. Each cutoff is
#' evaluated independently (a patient may fail several) and the number
#' removed per cutoff is reported in a message. Missing values fail a
#' cutoff. No columns are added.
#'
#' @param patient_df Per-patient table with `sum_stroma_area_um2`,
#'   `sum_tumor_area_um2` and, if `adequate_bcg` is set, `Adequate_BCG`.
#' @param filters Named list of cutoffs; see [cfg$aggregation$patient_filters]. `NULL`
#'   applies none.
#'
#' @return `patient_df` restricted to patients meeting every enabled cutoff.
filter_patients <- function(patient_df, filters = cfg$aggregation$patient_filters) {
  if (is.null(filters)) return(patient_df)
  unknown <- setdiff(names(filters), names(cfg$aggregation$patient_filters))
  if (length(unknown) > 0) stop("Unknown patient filter(s): ", paste(unknown, collapse = ", "))
  need <- function(col) {
    if (!col %in% names(patient_df)) stop("Patient filter requires column `", col, "`.")
  }

  pass <- list()
  min_cut <- function(name, x) {
    cut <- filters[[name]]
    if (!is.null(cut)) pass[[name]] <<- !is.na(x) & x > cut
  }
  if (!is.null(filters[["min_sum_stroma_area_um2"]])) need("sum_stroma_area_um2")
  if (!is.null(filters[["min_sum_tumor_area_um2"]])) need("sum_tumor_area_um2")
  if (!is.null(filters[["min_sum_tissue_area_um2"]])) { need("sum_stroma_area_um2"); need("sum_tumor_area_um2") }
  min_cut("min_sum_stroma_area_um2", patient_df$sum_stroma_area_um2)
  min_cut("min_sum_tumor_area_um2", patient_df$sum_tumor_area_um2)
  min_cut("min_sum_tissue_area_um2", patient_df$sum_stroma_area_um2 + patient_df$sum_tumor_area_um2)
  if (!is.null(filters[["adequate_bcg"]])) {
    need("Adequate_BCG")
    pass[["adequate_bcg"]] <- patient_df$Adequate_BCG %in% filters[["adequate_bcg"]]
  }

  for (nm in names(pass)) message("Patient filter ", nm, ": removes ", sum(!pass[[nm]]), " patient(s)")
  keep <- Reduce(`&`, pass, init = rep(TRUE, nrow(patient_df)))
  patient_df[keep, , drop = FALSE]
}

#' Full pipeline: filter, aggregate, and derive final column names
#'
#' Chains [filter_obj_region_tissue_level_data()] and
#' [summarize_summarized_msi_data_per_patient_v2()], then renames
#' `final_*` columns to the `nerve_*` names used downstream.
#'
#' @param ocd Unfiltered per-object marker/count table (object-count data).
#' @param tsd Per-region tissue segmentation table.
#' @param tsds Per-image tissue segmentation summary table.
#' @param nerve_obj_table Clustered per-object table (see
#'   [filter_obj_region_tissue_level_data()]).
#' @param clin_df_matched One row per image, with
#'   `im3_filename`/`PA_number`/`StudyID`/`MSI_tag`/`Cdj_tag`.
#' @param clin_df_compact One row per patient.
#' @param image_annotations Optional per-image annotation table with
#'   `Sample Name` and `include_final_imaging_cohort`; see
#'   [summarize_summarized_msi_data_per_patient_v2()]. If `NULL`, all
#'   images are retained.
#' @param regions_to_omit Optional data frame of regions to exclude,
#'   passed through to [filter_obj_region_tissue_level_data()].
#' @param tissue_category Tissue category to retain (default `"Stroma"`).
#' @param nerve_subclusters Numeric vector of subcluster indices to
#'   retain. Defaults to [cfg$aggregation$nerve_clusters_final]; see
#'   [filter_obj_region_tissue_level_data()] for the k-specificity caveat.
#' @param patient_filters Patient-level cutoffs applied to the aggregated
#'   table (default [cfg$aggregation$patient_filters]); see [filter_patients()]. `NULL`
#'   applies none.
#'
#' @return A per-patient data frame with `nerve_obj_area_percent`,
#'   `nerve_obj_density_per_mm2`, `nerve_obj_area_percent_subcluster_N`,
#'   and merged clinical columns, restricted to patients meeting
#'   `patient_filters`. Patient-level inclusion flags and reasons are
#'   provided with the supplementary tables, not generated here.
#'
#' @examples
#' \dontrun{
#' patient_df <- run_patient_nerve_quant_pipeline(
#'   ocd = obj_count_data,
#'   tsd = tissue_seg_data,
#'   tsds = tissue_seg_data_summary,
#'   nerve_obj_table = kmeans_result$nerve_obj_table,
#'   regions_to_omit = kmeans_result$regions_to_omit,
#'   clin_df_matched = read.delim("clin_msi_matched.tsv"),
#'   clin_df_compact = read.delim("clin_compact.tsv"),
#'   image_annotations = table_s4$image_annotations
#' )
#' }
run_patient_nerve_quant_pipeline <- function(ocd,
                                             tsd,
                                             tsds,
                                             nerve_obj_table,
                                             clin_df_matched,
                                             clin_df_compact,
                                             image_annotations = NULL,
                                             regions_to_omit = NULL,
                                             tissue_category = "Stroma",
                                             nerve_subclusters = cfg$aggregation$nerve_clusters_final,
                                             patient_filters = cfg$aggregation$patient_filters) {

  # Stage 1: filter and reconstruct region/tissue-level nerve-object summaries.
  filter_res <- filter_obj_region_tissue_level_data(
    ocd = ocd, tsd = tsd, tsds = tsds,
    nerve_obj_table = nerve_obj_table,
    tissue_category = tissue_category,
    nerve_subclusters = nerve_subclusters,
    min_region_area_um2 = cfg$aggregation$min_region_area_um2,
    max_region_obj_density = cfg$aggregation$max_region_obj_density,
    min_obj_area_um2 = cfg$aggregation$min_obj_area_um2,
    max_obj_area_um2 = cfg$aggregation$max_obj_area_um2,
    regions_to_omit = regions_to_omit,
    include_obj_from_NAregions = FALSE)

  object_df <- filter_res$filtered_tissue_summary %>% unique.data.frame()

  # Stage 2: aggregate to per-patient nerve quantification.
  df <- summarize_summarized_msi_data_per_patient_v2(
    obj_df = object_df,
    tissue_seg_df = tsds,
    clin_df_matched = clin_df_matched,
    clin_df_compact = clin_df_compact,
    image_annotations = image_annotations,
    tissue_category = tissue_category,
    min_stroma_area_pct = cfg$aggregation$min_stroma_area_percent,
    min_tumor_area_pct = cfg$aggregation$min_tumor_area_percent,
    min_combined_area_um2 = cfg$aggregation$min_combined_tissue_area_um2,
    max_object_density = cfg$aggregation$max_obj_density_per_mm2,
    max_n_obj = cfg$aggregation$max_n_objs,
    additional_summary_cols = "none",
    additional_summary = NULL,
    cols_keep = NULL
  )

  message(length(unique(df$PA_number)), " patients pass QC")

  # Rename to final nerve-quantification column names.
  df <- df %>%
    dplyr::rename(nerve_obj_area_percent = final_obj_area_percent,
                  nerve_obj_density_per_mm2 = final_obj_density_per_mm2) %>%
    dplyr::rename_with(~ sub("^final_obj_area_percent_subcluster_",
                             "nerve_obj_area_percent_subcluster_", .x),
                       dplyr::starts_with("final_obj_area_percent_subcluster_"))

  df <- filter_patients(df, patient_filters)
  message(length(unique(df$PA_number)), " patients after patient-level filters")

  df
}
