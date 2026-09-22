# Patient-level publication data. Preserve source columns and add explicit aliases
# required by the existing analysis functions; do not impose cohort eligibility here.
read_table_s1 <- function(path) {
  data <- as.data.frame(readxl::read_excel(
    path, sheet = "Table_S1_Cohort_demographics", .name_repair = "minimal"
  ))
  if (anyNA(data$StudyID) || anyDuplicated(data$StudyID)) {
    stop("Table S1 requires one nonmissing StudyID per patient.")
  }
  data$Cohort <- data$BRS_RNAseq_cohort
  data$Cohort[data$Cohort == "Not sequenced"] <- NA_character_
  data$nerve_obj_area_percent <- data$total_STaN_pos_area_percent
  for (k in c(1, 2, 4, 5, 6, 7)) {
    data[[paste0("nerve_obj_area_percent_subcluster_", k)]] <-
      data[[paste0("STaN_pos_area_percent_cluster_", k)]]
    data[[paste0("proportion_obj_in_subcluster_", k)]] <-
      data[[paste0("proportion_obj_in_STaN_cluster_", k)]]
  }
  data
}

# S4F maps source image filenames to patients. Keep the broader mapping intact:
# run_spatial_tweedie() applies the published image flags and clinical eligibility.
read_table_s4_image_annotations <- function(path) {
  data <- as.data.frame(readxl::read_excel(
    path, sheet = "S4F_Image_annotations", .name_repair = "minimal"
  ))
  if (anyNA(data[["Sample Name"]]) || anyDuplicated(data[["Sample Name"]])) {
    stop("Table S4F requires one nonmissing Sample Name per image.")
  }
  data
}
