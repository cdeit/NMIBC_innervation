# Relabel connected components of vessel_filtered.tif after the object size floor.
# Instance IDs here are newly assigned and need not equal upstream CD31 IDs.
# No extra glass-contact cleanup is applied to this already filtered mask.

label_vessel_objects <- function(image_id, cfg) {
  paths <- layer_paths(image_id, cfg)
  m <- enr_read_binary(paths$vessel)
  if (sum(m) == 0) {
    return(list(objects = tibble::tibble(), boundary_points = tibble::tibble()))
  }

  lab <- EBImage::bwlabel(EBImage::Image(m * 1.0, colormode = "Grayscale"))
  lab_mat <- EBImage::imageData(lab)
  n_raw <- max(lab_mat)
  if (n_raw == 0) {
    return(list(objects = tibble::tibble(), boundary_points = tibble::tibble()))
  }

  sizes_raw <- tabulate(lab_mat[lab_mat > 0], nbins = n_raw)
  keep <- which(sizes_raw >= cfg$min_vessel_object_area_px)
  if (length(keep) == 0) {
    return(list(objects = tibble::tibble(), boundary_points = tibble::tibble()))
  }

  filtered_mask <- matrix(lab_mat %in% keep, nrow(lab_mat), ncol(lab_mat))
  lab2 <- EBImage::bwlabel(EBImage::Image(filtered_mask * 1.0, colormode = "Grayscale"))
  lab_mat2 <- EBImage::imageData(lab2)
  n <- max(lab_mat2)

  ids <- seq_len(n)
  areas_px <- tabulate(lab_mat2[lab_mat2 > 0], nbins = n)

  cent <- t(vapply(ids, function(id) {
    idx <- which(lab_mat2 == id, arr.ind = TRUE)
    c(mean(idx[, 2]), mean(idx[, 1]))
  }, numeric(2)))

  vessel_object_id <- paste0(image_id, "_vessel_", ids)

  objects <- tibble::tibble(
    vessel_object_id = vessel_object_id,
    image_id           = image_id,
    local_id             = ids,
    centroid_x_px          = cent[, 1],
    centroid_y_px             = cent[, 2],
    area_px                     = areas_px,
    area_um2                      = areas_px * cfg$um_per_px^2
  )

  oc <- EBImage::ocontour(lab2)
  boundary_list <- list()
  for (id in ids) {
    if (id > length(oc) || is.null(oc[[id]]) || nrow(oc[[id]]) == 0) next
    pts <- oc[[id]]
    boundary_list[[length(boundary_list) + 1]] <- tibble::tibble(
      vessel_object_id = vessel_object_id[id],
      row_px              = pts[, 1] + 1L,
      col_px                = pts[, 2] + 1L
    )
  }

  list(
    objects         = objects,
    boundary_points = dplyr::bind_rows(boundary_list),
    lab_mat         = lab_mat2   ## kept for R/vessel_morphology.R's shape_features_df()
  )
}
