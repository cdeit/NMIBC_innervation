# Object-to-border and object-to-object distances.
# Border distances use the downsampled raster approximation in geometry.R.

object_to_border_distances_suffixed <- function(boundary_points, id_col, border_mask,
                                                  suffix, um_per_px, cfg) {
  out <- object_to_border_distances(boundary_points, id_col, border_mask, um_per_px, cfg)
  rename_map <- c(
    dist_nearest_um           = paste0("dist_nearest_", suffix, "_um"),
    nearest_segment_length_um = paste0("nearest_", suffix, "_segment_length_um"),
    dist_longest_um           = paste0("dist_longest_", suffix, "_um"),
    longest_segment_length_um = paste0("longest_", suffix, "_segment_length_um")
  )
  names(out)[match(names(rename_map), names(out))] <- rename_map
  out
}

object_to_border_distances <- function(boundary_points, id_col, border_mask, um_per_px, cfg) {
  ids <- unique(boundary_points[[id_col]])

  empty_result <- function() {
    tibble::tibble(
      !!id_col := ids,
      dist_nearest_um = NA_real_, nearest_segment_length_um = NA_real_,
      dist_longest_um = NA_real_, longest_segment_length_um = NA_real_
    )
  }
  if (sum(border_mask) == 0 || length(ids) == 0) return(empty_result())

  seg <- label_border_segments(border_mask, um_per_px)
  n_seg <- length(seg$lengths_um)
  if (n_seg == 0) return(empty_result())

  seg_dists <- lapply(seq_len(n_seg), function(s) {
    fast_distmap_um(seg$lab_mat == s, um_per_px, cfg$distmap_downsample_factor)
  })

  longest_id <- seg$longest_id
  longest_length_um <- seg$lengths_um[longest_id]

  rows <- lapply(ids, function(oid) {
    pts <- boundary_points[boundary_points[[id_col]] == oid, c("row_px", "col_px")]
    if (nrow(pts) == 0) {
      return(tibble::tibble(dist_nearest_um = NA_real_, nearest_segment_length_um = NA_real_,
                             dist_longest_um = NA_real_, longest_segment_length_um = longest_length_um))
    }
    idx <- cbind(pts$row_px, pts$col_px)
    per_seg_min <- vapply(seq_len(n_seg), function(s) min(seg_dists[[s]][idx]), numeric(1))
    nearest_seg <- which.min(per_seg_min)
    tibble::tibble(
      dist_nearest_um           = per_seg_min[nearest_seg],
      nearest_segment_length_um = seg$lengths_um[nearest_seg],
      dist_longest_um           = per_seg_min[longest_id],
      longest_segment_length_um = longest_length_um
    )
  })

  out <- dplyr::bind_rows(rows)
  out[[id_col]] <- ids
  dplyr::relocate(out, !!id_col)
}

