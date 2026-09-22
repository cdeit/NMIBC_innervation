#' read_table_s4.R
#'
#' Reads Table S4 (segmented nerve objects) and returns the inputs used by
#' `nerve_clustering.R` and `nerve_aggregation.R`, so the analysis
#' can be re-run from the supplementary table.
#'
#' @section Sheets used:
#' \describe{
#'   \item{S4D_Preliminary_obj_clustering}{All Stage 1 objects (marker means, `object_tag`).}
#'   \item{S4C_Candidate_nerve_objs_all}{Archived Stage 2 assignments.}
#'   \item{S4B_Final_cohort_with_cluster3}{Full per-object quantification of the
#'     nerve-candidate objects.}
#'   \item{S4E_Tissue_regions}{Region geometry and area for every image.}
#'   \item{S4F_Image_annotations}{`include_final_imaging_cohort` and
#'     `reason_exclude_final_imaging_cohort` for every image.}
#' }
#'
#' Patient-level clinical tables (`clin_df_matched`, `clin_df_compact`) are
#' not part of Table S4 and are supplied separately to
#' [run_patient_nerve_quant_pipeline()].
NULL

library(dplyr)
library(readxl)
library(stringr)

#' Load Table S4 as pipeline inputs
#'
#' @param path Path to `Table_S4_segmented_PGP95_SYP_objects.xlsx`.
#'
#' @return A list with:
#' \describe{
#'   \item{object_count_df}{Stage 1 object table for
#'     [run_nerve_object_kmeans_pipeline()] (`object_count_df`).}
#'   \item{tissue_seg_df}{Region-level table (`tissue_seg_df`/`tsd`).}
#'   \item{publication_reference}{Archived Stage 2 assignments
#'     (`publication_reference`).}
#'   \item{nerve_object_df}{Per-object quantification of nerve-candidate
#'     objects (`ocd` in [run_patient_nerve_quant_pipeline()]).}
#'   \item{tissue_seg_summary}{Image-level tissue areas (`tsds`).}
#'   \item{S4F_Image_annotations}{Per-image cohort annotations.}
#' }
load_table_s4 <- function(path) {
  required_sheets <- c(
    "S4D_Preliminary_obj_clustering", "S4C_Candidate_nerve_objs_all",
    "S4B_Final_cohort_with_cluster3", "S4E_Tissue_regions", "S4F_Image_annotations"
  )
  missing_sheets <- setdiff(required_sheets, readxl::excel_sheets(path))
  if (length(missing_sheets)) {
    stop("Table S4 is missing required sheets: ", paste(missing_sheets, collapse = ", "),
      ". Supply the complete publication input; see docs/REPRODUCIBILITY.md.",
      call. = FALSE
    )
  }
  rd <- function(sheet) as.data.frame(readxl::read_excel(path, sheet = sheet, .name_repair = "minimal"),
    check.names = FALSE
  )

  # Stage 1 objects: region and object identifiers are recovered from
  # object_tag (<image>_<category>_region<id>_object<id>_<row>).
  prelim <- rd("S4D_Preliminary_obj_clustering")
  ids <- stringr::str_match(prelim$object_tag, "_region(\\d+)_object(\\d+)_\\d+$")
  object_count_df <- prelim
  object_count_df$`Category Region ID.obj_count_data` <- as.numeric(ids[, 2])
  object_count_df$`Object ID.obj_count_data` <- as.numeric(ids[, 3])

  cand <- rd("S4C_Candidate_nerve_objs_all")
  publication_reference <- cand[, intersect(
    c("object_tag", "subcluster_num", "phenotype", "cluster_phenotype", "cell_type"), names(cand)
  )]

  # Quantification columns for nerve-candidate objects; the stored cluster
  # columns are dropped so they are re-derived from the pipeline output.
  nerve_object_df <- rd("S4B_Final_cohort_with_cluster3")
  nerve_object_df <- nerve_object_df[, setdiff(names(nerve_object_df), c(".subcluster_num", "cluster_phenotype"))]

  tissue_seg_df <- rd("S4E_Tissue_regions")
  tissue_seg_df <- tissue_seg_df[, setdiff(
    names(tissue_seg_df),
    c("passes_stage1_region_area_filter", "include_final_imaging_cohort")
  )]

  list(
    object_count_df = object_count_df,
    tissue_seg_df = tissue_seg_df,
    publication_reference = publication_reference,
    nerve_object_df = nerve_object_df,
    tissue_seg_summary = summarize_tissue_areas(tissue_seg_df),
    image_annotations = rd("S4F_Image_annotations")
  )
}

#' Image-level tissue areas from region-level data
#'
#' Tumor and stroma area per image, with each expressed as a percentage of
#' total tissue (tumor + stroma). One row per image and tissue category, as
#' expected for `tsds` by `nerve_aggregation.R`.
#'
#' @param tissue_seg_df Region-level table with `Sample Name`,
#'   `Tissue Category`, and `Region Area (square microns).tissue_seg_data`.
#'
#' @return Data frame of image-level areas.
summarize_tissue_areas <- function(tissue_seg_df) {
  area <- tissue_seg_df$`Region Area (square microns).tissue_seg_data`
  per_image <- tissue_seg_df %>%
    dplyr::mutate(area = area) %>%
    dplyr::group_by(`Sample Name`) %>%
    dplyr::summarize(
      tumor = sum(area[`Tissue Category` == "Tumor"]),
      stroma = sum(area[`Tissue Category` == "Stroma"]), .groups = "drop"
    ) %>%
    dplyr::mutate(
      `Tumor Area (um^2).tissue_seg_data_summary` = tumor,
      `Stroma Area (um^2).tissue_seg_data_summary` = stroma,
      `Tumor+Stroma Area (um^2).tissue_seg_data_summary` = tumor + stroma,
      `Tumor Area (%).tissue_seg_data_summary` = 100 * tumor / (tumor + stroma),
      `Stroma Area (%).tissue_seg_data_summary` = 100 * stroma / (tumor + stroma)
    ) %>%
    dplyr::select(-tumor, -stroma)
  tidyr::crossing(per_image, `Tissue Category` = c("Stroma", "Tumor")) %>%
    dplyr::select(`Sample Name`, `Tissue Category`, dplyr::everything()) %>%
    as.data.frame(check.names = FALSE)
}
