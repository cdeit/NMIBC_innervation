# Label each subtype independently after glass/tissue cleanup.
# ocontour returns zero-based matrix coordinates; +1 converts to R indices.
# Centroids are descriptive columns; distances use contour points.

label_nerve_objects <- function(image_id, cfg, glass_mask, phenotyped) {

  objects_list  <- list()
  boundary_list <- list()
  filled_masks  <- list()   ## post-debris-removal filled masks, per subtype --

  for (code in cfg$nerve_codes) {
    code_chr <- as.character(code)
    # The coded mask is the exact partition of the six supplied subtype masks.
    m <- phenotyped == code
    m <- remove_debris(m, glass_mask, cfg$max_glass_contact_frac, cfg$min_tissue_contact_frac,
                        min_area_px = cfg$min_nerve_object_area_px)
    filled_masks[[code_chr]] <- m
    if (sum(m) == 0) next

    lab <- EBImage::bwlabel(EBImage::Image(m * 1.0, colormode = "Grayscale"))
    lab_mat <- EBImage::imageData(lab)
    n <- max(lab_mat)
    if (n == 0) next

    ids <- seq_len(n)
    areas_px <- tabulate(lab_mat[lab_mat > 0], nbins = n)

    cent <- t(vapply(ids, function(id) {
      idx <- which(lab_mat == id, arr.ind = TRUE)
      c(mean(idx[, 2]), mean(idx[, 1]))   # x_px = col, y_px = row (1-indexed)
    }, numeric(2)))

    nerve_object_id <- paste0(image_id, "_nerve_sc", code, "_", ids)

    objects_list[[code_chr]] <- tibble::tibble(
      nerve_object_id = nerve_object_id,
      image_id         = image_id,
      subtype_code      = as.integer(code),
      local_id           = ids,
      centroid_x_px        = cent[, 1],
      centroid_y_px          = cent[, 2],
      area_px                  = areas_px,
      area_um2                   = areas_px * cfg$um_per_px^2
    )

    oc <- EBImage::ocontour(lab)
    for (id in ids) {
      if (id > length(oc) || is.null(oc[[id]]) || nrow(oc[[id]]) == 0) next
      pts <- oc[[id]]
      boundary_list[[length(boundary_list) + 1]] <- tibble::tibble(
        nerve_object_id = nerve_object_id[id],
        subtype_code     = as.integer(code),
        row_px             = pts[, 1] + 1L,
        col_px               = pts[, 2] + 1L
      )
    }
  }

  list(
    objects         = dplyr::bind_rows(objects_list),
    boundary_points = dplyr::bind_rows(boundary_list),
    filled_masks    = filled_masks
  )
}
