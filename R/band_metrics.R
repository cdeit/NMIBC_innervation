## band_metrics.R -- generalized band tabulation (Parts 1 & 2) and radius-windowed
## enrichment. Replaces the prior nerve-only compute_enrichment()/calc_radius_enrichment().

## Tallies area/density per distance band, for one or more mask "layers"
## (nerve subtypes, vessel) against one reference surface (tumor-stroma border
## or vessel mask) in a single image.
##
##   distmap_um    - numeric matrix, um distance to the reference surface
##   band_breaks_um- numeric vector, band edges (length n_bands + 1)
##   band_labels   - character vector, length n_bands, aligned with band_breaks_um
##   eligible_mask - logical matrix: pixels assignable to a band at all
##   denom_mask    - logical matrix: density denominator (stroma). Kept separate
##                   from eligible_mask so a future denominator change doesn't
##                   require touching the tallying logic, though in this
##                   pipeline the two are currently identical.
##   layer_masks   - named list of matrices: binary (logical) or coded (integer)
##   layer_kind    - named list matching layer_masks: "binary" | "coded"
##   layer_codes   - named list, REQUIRED for "coded" layers: the full set of
##                   codes that must always get a row (even if absent = 0 area)
##                   so downstream joins never confuse "not measured" with
##                   "measured, zero".
##
## Returns one row per (image_id, reference_name, band_index, layer_name, subtype_code).
compute_band_metrics <- function(distmap_um, band_breaks_um, band_labels,
                                  eligible_mask, denom_mask,
                                  layer_masks, layer_kind, layer_codes = list(),
                                  image_id, reference_name, um_per_px) {

  n_bands <- length(band_breaks_um) - 1L
  stopifnot(length(band_labels) == n_bands)

  eligible_idx <- which(as.vector(eligible_mask))
  band_idx_eligible <- cut(as.vector(distmap_um)[eligible_idx],
                            breaks = band_breaks_um, labels = FALSE, include.lowest = TRUE)

  denom_idx <- which(as.vector(denom_mask))
  band_idx_denom <- cut(as.vector(distmap_um)[denom_idx],
                         breaks = band_breaks_um, labels = FALSE, include.lowest = TRUE)
  denom_area_px <- tabulate(band_idx_denom, nbins = n_bands)
  denom_area_um2 <- denom_area_px * um_per_px^2

  make_rows <- function(area_px_band, layer_name, subtype_code) {
    tibble::tibble(
      image_id       = image_id,
      reference_name = reference_name,
      band_index     = seq_len(n_bands),
      band_label     = band_labels,
      band_inner_um  = band_breaks_um[-length(band_breaks_um)],
      band_outer_um  = band_breaks_um[-1],
      layer_name     = layer_name,
      subtype_code   = subtype_code,
      area_px        = area_px_band,
      area_um2       = area_px_band * um_per_px^2,
      stroma_area_um2 = denom_area_um2,
      band_density   = dplyr::if_else(denom_area_um2 > 0,
                                       (area_px_band * um_per_px^2) / denom_area_um2,
                                       NA_real_)
    )
  }

  results <- list()
  for (layer_name in names(layer_masks)) {
    mask <- layer_masks[[layer_name]]
    kind <- layer_kind[[layer_name]]

    if (kind == "binary") {
      vals <- as.vector(mask)[eligible_idx]
      area_px_band <- tabulate(band_idx_eligible[vals], nbins = n_bands)
      results[[length(results) + 1]] <- make_rows(area_px_band, layer_name, NA_integer_)

    } else if (kind == "coded") {
      codes <- layer_codes[[layer_name]]
      if (is.null(codes)) stop("layer_codes[['", layer_name, "']] is required for coded layers")
      code_vals <- as.vector(mask)[eligible_idx]
      for (code in codes) {
        present <- code_vals == code
        area_px_band <- tabulate(band_idx_eligible[present], nbins = n_bands)
        results[[length(results) + 1]] <- make_rows(area_px_band, layer_name, as.integer(code))
      }
    } else {
      stop("Unknown layer_kind '", kind, "' for layer '", layer_name, "'")
    }
  }
  dplyr::bind_rows(results)
}

