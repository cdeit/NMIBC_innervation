# Vessel shape measurements on newly labeled vessel objects.
# Morphological closing is used as a hull approximation in named legacy features.
# Pixel-unit feature columns are retained to match the historical table contract.

shape_features_df <- function(labeled_eb, instance_mat = NULL) {
  sh  <- as.data.frame(computeFeatures.shape(labeled_eb))
  mom <- as.data.frame(computeFeatures.moment(labeled_eb))
  sh$circularity <- (4 * pi * sh$s.area) / (sh$s.perimeter ^ 2)
  sh$pa_ratio <- sh$s.perimeter / sh$s.area
  sh$radius_cv   <- sh$s.radius.sd / sh$s.radius.mean
  sh$roundness   <- (4 * sh$s.area) / (pi * mom$m.majoraxis ^ 2)
  sh$eccentricity <- mom$m.eccentricity
  sh$pa_ecc_ratio <- sh$pa_ratio / sh$eccentricity


  if (!is.null(instance_mat)) {
    obj_ids <- sort(unique(as.vector(instance_mat)))
    obj_ids <- obj_ids[obj_ids > 0]
    sh$object_id <- obj_ids   # actual TIF instance IDs, not seq_len

    get_bbox <- function(i, pad = 5L) {
      idx <- which(instance_mat == i, arr.ind = TRUE)
      list(
        rmin = max(1L,                  min(idx[,1]) - pad),
        rmax = min(nrow(instance_mat),  max(idx[,1]) + pad),
        cmin = max(1L,                  min(idx[,2]) - pad),
        cmax = min(ncol(instance_mat),  max(idx[,2]) + pad)
      )
    }

    sh$edge_granularity <- vapply(obj_ids, function(i) {
      bb       <- get_bbox(i, pad = 5L)
      ob       <- (instance_mat[bb$rmin:bb$rmax, bb$cmin:bb$cmax] == i) * 1L
      hull     <- mmand::closing(ob, mmand::shapeKernel(c(25, 25), type = "disc"))
      defects  <- hull - ob   # pixels in hull but not in object
      if (sum(defects) == 0) return(0)
      defect_labeled <- bwlabel(Image(defects * 1.0, colormode = "Grayscale"))
      defect_sizes   <- tabulate(imageData(defect_labeled)[imageData(defect_labeled) > 0])
      mean(defect_sizes)   # average size of chunky protrusions
    }, numeric(1))

    sh$smoothness <- vapply(obj_ids, function(i) {
      bb        <- get_bbox(i, pad = 5L)
      ob        <- (instance_mat[bb$rmin:bb$rmax, bb$cmin:bb$cmax] == i) * 1L
      hull      <- mmand::closing(ob, mmand::shapeKernel(c(25, 25), type = "disc"))
      hull_img  <- Image(hull * 1.0, colormode = "Grayscale")
      hull_perim <- computeFeatures.shape(bwlabel(hull_img))[1, "s.perimeter"]
      actual_perim <- sh$s.perimeter[sh$object_id == i]
      if (is.null(hull_perim) || hull_perim == 0) return(NA_real_)
      hull_perim / actual_perim
    }, numeric(1))


    sh$branch_count <- vapply(obj_ids, function(i) {
      bb  <- get_bbox(i, pad = 3L)
      ob  <- (instance_mat[bb$rmin:bb$rmax, bb$cmin:bb$cmax] == i) * 1L

      skel <- mmand::skeletonise(ob, mmand::shapeKernel(c(3,3), type = "box"))

      if (sum(skel) == 0) return(0L)

      skel_img   <- Image(skel * 1.0, colormode = "Grayscale")
      neighbour_count <- filter2(skel_img, matrix(1, 3, 3)) - skel_img  # exclude self
      branch_pts <- sum(imageData(neighbour_count) >= 3 & skel > 0)
      as.integer(branch_pts)
    }, integer(1))


    sh$lumen_fraction <- vapply(obj_ids, function(i) {
      bb  <- get_bbox(i)
      ob  <- instance_mat[bb$rmin:bb$rmax, bb$cmin:bb$cmax] == i
      ofi <- fillHull(Image(ob * 1.0, colormode = "Grayscale"))
      af  <- sum(imageData(ofi) > 0.5)
      if (af > 0) 1 - sum(ob) / af else 0
    }, numeric(1))

    sh$solidity <- vapply(obj_ids, function(i) {
      bb        <- get_bbox(i, pad = 10L)
      ob        <- (instance_mat[bb$rmin:bb$rmax, bb$cmin:bb$cmax] == i) * 1L
      hull      <- mmand::closing(ob, mmand::shapeKernel(c(7, 7), type = "disc"))
      hull_area <- sum(hull > 0)
      if (hull_area > 0) sum(ob > 0) / hull_area else 0
    }, numeric(1))

  } else {
    sh$object_id <- seq_len(nrow(sh))
  }

  sh <- sh %>% dplyr::relocate(object_id, .before = s.area)

  sh
}

compute_vessel_shape_features <- function(vessel_labeling, image_id) {
  lab_mat <- vessel_labeling$lab_mat
  if (is.null(lab_mat) || max(lab_mat) == 0) return(tibble::tibble())

  feats <- shape_features_df(lab_mat, instance_mat = lab_mat)
  feats <- tibble::as_tibble(feats)
  feats$vessel_object_id <- paste0(image_id, "_vessel_", feats$object_id)
  feats
}
