# Connected-component removal using glass/tissue contact and conditional area floors.
# Floating-point filter2 results feed a strict glass_contact > 0 comparison.
# Preserve these calculations and the validated runtime for historical agreement.

remove_debris <- function(bin_mat, exclusion_mat,
                           max_glass_contact_frac,
                           min_tissue_contact_frac,
                           min_area_px = NULL,
                           tissue_mat = NULL) {

  if (sum(bin_mat > 0, na.rm = TRUE) == 0) return(bin_mat > 0)

  bin_mat <- bin_mat > 0
  exclusion_mat <- exclusion_mat > 0

  if (is.null(tissue_mat)) {
    tissue_mat <- !exclusion_mat
  } else {
    tissue_mat <- tissue_mat > 0
  }
  tissue_mat <- tissue_mat & !exclusion_mat

  # Label full connected components before evaluating their boundaries.
  lab_img <- EBImage::bwlabel(EBImage::Image(bin_mat * 1, colormode = "Grayscale"))
  lab_mat <- EBImage::imageData(lab_img)
  n_comp <- max(lab_mat)
  if (n_comp == 0) return(bin_mat)

  glass_env  <- exclusion_mat & !bin_mat
  tissue_env <- tissue_mat

  # Count the eight neighbors; retain FFT arithmetic without rounding.
  kern <- matrix(1, 3, 3)
  kern[2, 2] <- 0

  glass_nb  <- EBImage::imageData(EBImage::filter2(EBImage::Image(glass_env * 1, colormode = "Grayscale"), kern))
  tissue_nb <- EBImage::imageData(EBImage::filter2(EBImage::Image(tissue_env * 1, colormode = "Grayscale"), kern))

  # Component boundary = foreground removed by the size-3 erosion.
  bin_eroded <- EBImage::imageData(EBImage::erode(
    EBImage::Image(bin_mat * 1, colormode = "Grayscale"),
    EBImage::makeBrush(3, shape = "disc")
  )) > 0.5

  border_pixels <- lab_mat > 0 & !bin_eroded
  comp_idx   <- which(border_pixels)
  comp_label <- lab_mat[comp_idx]

  n_glass  <- tapply(glass_nb[comp_idx],  comp_label, sum)
  n_tissue <- tapply(tissue_nb[comp_idx], comp_label, sum)

  glass_contact  <- numeric(n_comp)
  tissue_contact <- numeric(n_comp)
  glass_contact[as.integer(names(n_glass))]   <- as.numeric(n_glass)
  tissue_contact[as.integer(names(n_tissue))] <- as.numeric(n_tissue)

  # Fractions summarize all component-boundary neighbor contacts.
  denom <- pmax(glass_contact + tissue_contact, 1)
  glass_frac  <- glass_contact  / denom
  tissue_frac <- tissue_contact / denom

  # Separate complete-surround criterion using a size-7 glass dilation.
  glass_surround <- EBImage::imageData(EBImage::dilate(
    EBImage::Image(glass_env * 1, colormode = "Grayscale"),
    EBImage::makeBrush(7, shape = "disc")
  )) > 0.5

  fully_surrounded <- vapply(seq_len(n_comp), function(id) {
    idx <- lab_mat == id
    if (!any(idx)) return(FALSE)
    all(glass_surround[idx])
  }, logical(1))

  comp_sizes <- tabulate(lab_mat[lab_mat > 0], nbins = n_comp)
  # Crucial historical condition: the size floor is not universal.
  too_small <- if (!is.null(min_area_px)) {
    (comp_sizes < min_area_px) & (glass_contact > 0)
  } else {
    rep(FALSE, n_comp)
  }

  remove_labels <- which(
    fully_surrounded |
      glass_frac >= max_glass_contact_frac |
      tissue_frac < min_tissue_contact_frac |
      too_small
  )

  lab_mat[lab_mat %in% remove_labels] <- 0L
  lab_mat > 0
}
