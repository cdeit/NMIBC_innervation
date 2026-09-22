# Raster border construction and approximate Euclidean distance transforms.
# makeBrush uses size (not radius); distance sampling is twofold by default.
# Segment length is pixel count times scale, not a subpixel arc-length estimate.

compute_border <- function(inside_mat, outside_mat, um_per_px,
                            dilate_px = 3L, valid_mat = NULL) {
  h <- nrow(inside_mat); w <- ncol(inside_mat)
  if (sum(inside_mat) == 0 || sum(outside_mat) == 0) {
    return(list(border = matrix(FALSE, h, w), border_length_um = 0))
  }
  # dilate_px is makeBrush(size), not a radius.
  outside_dil <- EBImage::imageData(EBImage::dilate(
    EBImage::Image(outside_mat * 1.0, colormode = "Grayscale"),
    EBImage::makeBrush(dilate_px, shape = "disc")
  )) > 0.5

  border <- inside_mat & outside_dil
  if (!is.null(valid_mat)) border <- border & valid_mat

  if (sum(border) == 0) {
    return(list(border = border, border_length_um = 0))
  }
  list(border = border, border_length_um = sum(border) * um_per_px)
}

compute_tumor_stroma_border <- function(tumor, stroma, glass, um_per_px, dilate_px = 3L) {
  compute_border(tumor, stroma, um_per_px, dilate_px, valid_mat = !glass)
}

compute_tumor_glass_border <- function(tumor, glass, um_per_px, dilate_px = 3L) {
  compute_border(tumor, glass, um_per_px, dilate_px, valid_mat = NULL)
}

compute_stroma_glass_border <- function(stroma, glass, um_per_px, dilate_px = 3L) {
  compute_border(stroma, glass, um_per_px, dilate_px, valid_mat = NULL)
}

label_border_segments <- function(border_mask, um_per_px) {
  h <- nrow(border_mask); w <- ncol(border_mask)
  if (sum(border_mask) == 0) {
    return(list(lab_mat = matrix(0L, h, w), lengths_um = numeric(0),
                longest_id = NA_integer_, longest_mask = border_mask,
                total_length_um = 0))
  }
  lab <- EBImage::bwlabel(EBImage::Image(border_mask * 1.0, colormode = "Grayscale"))
  lab_mat <- EBImage::imageData(lab)
  n <- max(lab_mat)
  sizes <- tabulate(lab_mat[lab_mat > 0], nbins = n)
  lengths_um <- sizes * um_per_px
  longest_id <- which.max(sizes)
  longest_mask <- lab_mat == longest_id
  list(lab_mat = lab_mat, lengths_um = lengths_um, longest_id = longest_id,
       longest_mask = longest_mask, total_length_um = sum(lengths_um))
}

fast_distmap_um <- function(bin_mat, um_per_px, factor = 2L) {
  h <- nrow(bin_mat); w <- ncol(bin_mat)
  # Preserve the sample origin and index rounding used by the original run.
  down <- bin_mat[seq(1, h, factor), seq(1, w, factor)]
  dt <- EBImage::imageData(EBImage::distmap(
    EBImage::Image((!down) * 1.0, colormode = "Grayscale")
  )) * um_per_px * factor
  dt[round(seq(1, nrow(dt), length.out = h)),
     round(seq(1, ncol(dt), length.out = w))]
}
