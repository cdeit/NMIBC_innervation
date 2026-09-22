# Stage A functions from the supplied segmentation-preprocessing script.
# Bodies preserved verbatim. See docs/PROVENANCE.md.

write_label_tif <- function(code_mat, path) {
  h <- nrow(code_mat); w <- ncol(code_mat)
  store_mat <- matrix(as.integer(-code_mat), h, w)
  pixel_bytes <- writeBin(as.integer(t(store_mat)), raw(), size = 4L, endian = "little")

  strip_offset <- 8L
  strip_bytes  <- length(pixel_bytes)
  ifd_offset   <- strip_offset + strip_bytes

  ## TIFF baseline tags for one uncompressed, single-strip, single-sample
  ## 32-bit signed-integer image: width/length/bits-per-sample/compression/
  ## photometric/strip-offset/samples-per-pixel/rows-per-strip/strip-byte-
  ## count/sample-format. Listed here in ascending tag-number order, which
  ## the TIFF spec requires.
  tags <- list(
    c(tag = 256L, type = 4L, count = 1L, value = w),             # ImageWidth
    c(tag = 257L, type = 4L, count = 1L, value = h),             # ImageLength
    c(tag = 258L, type = 3L, count = 1L, value = 32L),           # BitsPerSample
    c(tag = 259L, type = 3L, count = 1L, value = 1L),            # Compression: none
    c(tag = 262L, type = 3L, count = 1L, value = 1L),            # Photometric: black-is-zero
    c(tag = 273L, type = 4L, count = 1L, value = strip_offset),  # StripOffsets
    c(tag = 277L, type = 3L, count = 1L, value = 1L),            # SamplesPerPixel
    c(tag = 278L, type = 4L, count = 1L, value = h),             # RowsPerStrip
    c(tag = 279L, type = 4L, count = 1L, value = strip_bytes),   # StripByteCounts
    c(tag = 339L, type = 3L, count = 1L, value = 2L)             # SampleFormat: signed int
  )

  con <- file(path, "wb")
  on.exit(close(con))

  writeBin(charToRaw("II"), con)                                        # little-endian marker
  writeBin(42L,        con, size = 2L, endian = "little")                # TIFF magic number
  writeBin(ifd_offset, con, size = 4L, endian = "little")

  writeBin(pixel_bytes, con)

  writeBin(length(tags), con, size = 2L, endian = "little")
  for (t in tags) {
    writeBin(t[["tag"]],   con, size = 2L, endian = "little")
    writeBin(t[["type"]],  con, size = 2L, endian = "little")
    writeBin(t[["count"]], con, size = 4L, endian = "little")
    if (t[["type"]] == 3L) {
      ## SHORT values occupy the first 2 of the 4 value bytes.
      writeBin(t[["value"]], con, size = 2L, endian = "little")
      writeBin(0L,           con, size = 2L, endian = "little")
    } else {
      writeBin(t[["value"]], con, size = 4L, endian = "little")
    }
  }
  writeBin(0L, con, size = 4L, endian = "little")  # no further IFDs

  invisible(path)
}

read_tif_layer <- function(path, layer = 1L) {
  pages <- readTIFF(path, all = TRUE, as.is = TRUE)
  if (!is.list(pages)) pages <- list(pages)
  if (layer > length(pages)) {
    stop(sprintf("Layer %d requested but %s only has %d layer(s)", layer, basename(path), length(pages)))
  }
  m <- pages[[layer]]
  storage.mode(m) <- "integer"
  m
}

shape_features_df <- function(labeled_eb, instance_mat = NULL) {
  sh  <- as.data.frame(computeFeatures.shape(labeled_eb))
  mom <- as.data.frame(computeFeatures.moment(labeled_eb))
  sh$circularity  <- (4 * pi * sh$s.area) / (sh$s.perimeter ^ 2)
  sh$pa_ratio     <- sh$s.perimeter / sh$s.area
  sh$radius_cv    <- sh$s.radius.sd / sh$s.radius.mean
  sh$roundness    <- (4 * sh$s.area) / (pi * mom$m.majoraxis ^ 2)
  sh$eccentricity <- mom$m.eccentricity
  sh$pa_ecc_ratio <- sh$pa_ratio / sh$eccentricity

  if (!is.null(instance_mat)) {
    obj_ids <- sort(unique(as.vector(instance_mat)))
    obj_ids <- obj_ids[obj_ids > 0]
    sh$object_id <- obj_ids

    get_bbox <- function(i, pad = 5L) {
      idx <- which(instance_mat == i, arr.ind = TRUE)
      list(
        rmin = max(1L,                 min(idx[, 1]) - pad),
        rmax = min(nrow(instance_mat), max(idx[, 1]) + pad),
        cmin = max(1L,                 min(idx[, 2]) - pad),
        cmax = min(ncol(instance_mat), max(idx[, 2]) + pad)
      )
    }

    ## Convexity defects: pixels inside a closed hull but outside the object
    ## itself. Their mean size is a crude granularity score -- a smooth tube
    ## has a small number of tiny defects, a ragged blob has large ones.
    sh$edge_granularity <- vapply(obj_ids, function(i) {
      bb      <- get_bbox(i, pad = 5L)
      ob      <- (instance_mat[bb$rmin:bb$rmax, bb$cmin:bb$cmax] == i) * 1L
      hull    <- mmand::closing(ob, mmand::shapeKernel(c(25, 25), type = "disc"))
      defects <- hull - ob
      if (sum(defects) == 0) return(0)
      defect_labeled <- bwlabel(Image(defects * 1.0, colormode = "Grayscale"))
      defect_sizes   <- tabulate(imageData(defect_labeled)[imageData(defect_labeled) > 0])
      mean(defect_sizes)
    }, numeric(1))

    ## Hull perimeter vs. actual perimeter -- close to 1 for a smooth outline,
    ## higher for a jagged one.
    sh$smoothness <- vapply(obj_ids, function(i) {
      bb           <- get_bbox(i, pad = 5L)
      ob           <- (instance_mat[bb$rmin:bb$rmax, bb$cmin:bb$cmax] == i) * 1L
      hull         <- mmand::closing(ob, mmand::shapeKernel(c(25, 25), type = "disc"))
      hull_img     <- Image(hull * 1.0, colormode = "Grayscale")
      hull_perim   <- computeFeatures.shape(bwlabel(hull_img))[1, "s.perimeter"]
      actual_perim <- sh$s.perimeter[sh$object_id == i]
      if (is.null(hull_perim) || hull_perim == 0) return(NA_real_)
      hull_perim / actual_perim
    }, numeric(1))

    ## Skeleton branch points (3+ skeleton neighbors) -- distinguishes a
    ## single tube (0-2 branches) from a branchy or sheet-like blob.
    sh$branch_count <- vapply(obj_ids, function(i) {
      bb   <- get_bbox(i, pad = 3L)
      ob   <- (instance_mat[bb$rmin:bb$rmax, bb$cmin:bb$cmax] == i) * 1L
      skel <- mmand::skeletonise(ob, mmand::shapeKernel(c(3, 3), type = "box"))
      if (sum(skel) == 0) return(0L)
      skel_img        <- Image(skel * 1.0, colormode = "Grayscale")
      neighbour_count <- filter2(skel_img, matrix(1, 3, 3)) - skel_img
      as.integer(sum(imageData(neighbour_count) >= 3 & skel > 0))
    }, integer(1))

    ## Fraction of the filled-hull area that is NOT occupied by the object
    ## itself -- an open lumen (a vessel cross-section) scores high here.
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

  dplyr::relocate(sh, object_id, .before = s.area)
}

filter_vessels <- function(vessel_instance, image_id) {
  ## Shape-classification thresholds, tuned once against real vessel data and
  ## not treated as a user-configurable knob.
  min_area_px                <- 400
  max_area_px                <- 100000

  min_circularity             <- 0.2
  min_lumen_frac              <- 0.035
  min_lumen_frac_confident    <- 0.1

  min_solidity_nolumen        <- 0.93
  min_area_px_nolumen_round   <- 400
  max_round_nolumen_branches  <- 100

  min_elongated_area_px       <- 400
  max_radius_cv_elongated     <- 0.6
  min_solidity_elongated      <- 0.85
  min_circularity_elongated   <- 0.10
  min_eccentricity_elongated  <- 0.85
  max_pa_ecc_ratio             <- 0.2

  max_roundness_cshape        <- 0.4
  max_radius_cv_cshape        <- 0.5
  min_solidity_cshape         <- 0.87
  min_area_cshape              <- 500

  min_solidity_auto            <- 0.97

  min_area_linear              <- 100
  max_roundness_linear         <- 0.60
  max_radius_cv_linear         <- 0.75
  min_solidity_linear          <- 0.70
  min_eccentricity_linear      <- 0.85

  n_obj <- max(vessel_instance)
  if (n_obj == 0L) {
    return(list(mask = matrix(0L, nrow(vessel_instance), ncol(vessel_instance)), feats = NULL))
  }

  labeled <- Image(vessel_instance * 1.0, colormode = "Grayscale")
  lab_mat <- imageData(labeled)

  feats <- tryCatch(shape_features_df(labeled, lab_mat), error = function(e) {
    message("  [vessel] shape_features_df failed for ", image_id, ": ", e$message)
    NULL
  })
  if (is.null(feats)) return(list(mask = matrix(0L, nrow(vessel_instance), ncol(vessel_instance)), feats = NULL))
  feats$image_id <- image_id

  ## Round, with a confident lumen.
  feats$b1 <- feats$lumen_fraction >= min_lumen_frac &
    (feats$lumen_fraction >= min_lumen_frac_confident | feats$circularity >= min_circularity)

  ## Round, no resolvable lumen, but smooth and not too small; excludes
  ## branchy blobs that would otherwise pass on solidity alone.
  feats$b2 <- feats$circularity >= min_circularity &
    feats$lumen_fraction < min_lumen_frac &
    feats$solidity >= min_solidity_nolumen &
    feats$s.area >= min_area_px_nolumen_round &
    feats$pa_ecc_ratio <= 0.2
  feats$b2 <- feats$b2 & !(feats$lumen_fraction < 0.02 & feats$branch_count >= max_round_nolumen_branches)

  ## Elongated tube, or a near-perfectly linear object regardless of the
  ## other shape thresholds.
  feats$b3 <- feats$s.area >= min_elongated_area_px &
    ((feats$circularity >= min_circularity_elongated &
        feats$radius_cv <= max_radius_cv_elongated &
        feats$solidity >= min_solidity_elongated &
        feats$eccentricity >= min_eccentricity_elongated &
        feats$pa_ecc_ratio <= max_pa_ecc_ratio) |
       feats$eccentricity >= 0.98)

  ## C-shaped / partial cross-section: not round enough for b1/b2, but smooth
  ## and compact.
  feats$b4 <- feats$s.area >= min_area_cshape &
    feats$roundness <= max_roundness_cshape &
    feats$radius_cv <= max_radius_cv_cshape &
    feats$solidity >= min_solidity_cshape &
    feats$circularity < min_circularity

  ## Catch-all: extremely smooth outline regardless of the other metrics.
  feats$b5 <- feats$solidity >= min_solidity_auto

  ## Small linear fragments -- admitted through a lower area floor than b3.
  feats$b6 <- feats$s.area >= min_area_linear &
    feats$roundness <= max_roundness_linear &
    feats$radius_cv <= max_radius_cv_linear &
    feats$solidity >= min_solidity_linear &
    feats$eccentricity >= min_eccentricity_linear &
    feats$pa_ecc_ratio <= max_pa_ecc_ratio

  feats$is_vessel <- (feats$s.area >= min_area_px & feats$s.area <= max_area_px &
                         (feats$b1 | feats$b2 | feats$b3 | feats$b4 | feats$b5)) |
    (feats$s.area >= min_area_linear & feats$b6)

  vessel_ids <- feats$object_id[feats$is_vessel]
  out <- matrix(0L, nrow(vessel_instance), ncol(vessel_instance))
  out[lab_mat %in% vessel_ids] <- 1L

  list(mask = out, feats = feats)
}

map_nerve_phenotypes <- function(nerve_instance, nerve_df, valid_codes) {
  empty <- matrix(0L, nrow(nerve_instance), ncol(nerve_instance))
  n_obj <- max(nerve_instance)
  if (nrow(nerve_df) == 0L || n_obj == 0L) return(list(mask = empty))

  codes <- as.integer(gsub("subcluster_", "", nerve_df$phenotype))
  in_range <- !is.na(nerve_df$objID) & nerve_df$objID >= 1L & nerve_df$objID <= n_obj
  codes[!codes %in% valid_codes] <- 0L

  ## lookup[1] is background (raw instance label 0); lookup[k + 1] is the
  ## code for raw instance label k.
  lookup <- integer(n_obj + 1L)
  lookup[nerve_df$objID[in_range] + 1L] <- codes[in_range]

  list(mask = matrix(lookup[nerve_instance + 1L], nrow(nerve_instance), ncol(nerve_instance)))
}

