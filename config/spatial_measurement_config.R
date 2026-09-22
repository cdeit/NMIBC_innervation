# Parameters for tumor-stroma border bands/enrichment and object-border distances.
# Values match the successful Alpine probe. Sourcing this file runs no analysis.
cfg <- list(
  # Original image scale; one pixel occupies 0.25 square micrometers.
  um_per_px = 0.5,
  file_tumor = 'tumor.tif',
  file_stroma = 'stroma.tif',
  file_glass = 'glass.tif',
  file_vessel = 'vessel_filtered.tif',
  file_nerve_phenotyped = 'nerve_phenotyped.tif',
  nerve_codes = c(1L, 2L, 4L, 5L, 6L, 7L),

  # Component cleanup; size floors apply ONLY when computed glass_contact > 0.
  # R/debris.R retains the original size-3 erosion, eight-neighbor convolution
  # and size-7 dilation of glass. Do not round FFT results or add a tolerance.
  max_glass_contact_frac = 0.9,
  min_tissue_contact_frac = 0.1,
  min_tumor_area_px = 25000L,
  min_stroma_area_px = 25000L,

  # makeBrush SIZE (not radius). Raster distances are sampled every second
  # row/column, scaled by 0.5*2 um, then resampled using rounded indices.
  border_dilate_px = 3L,
  distmap_downsample_factor = 2L,

  # Actual intervals: [0,25], (25,50], ..., (300,Inf]. Updated-goal display labels
  # are retained for CSV agreement; numerical breaks determine membership.
  border_breaks = c(0,25,50,100,200,300,Inf),
  border_labels = c('0-25 um','26-50 um','51-100 um','101-200 um','201-300 um','>300 um'),
  # Nerve minimum is conditional on glass contact; vessel minimum is
  # unconditional, applied to connected components of vessel_filtered.tif.
  min_nerve_object_area_px = 5L,
  min_vessel_object_area_px = 30L,

  # Display labels in updated_format_goal_spatial_output, not inferred phenotypes.
  nerve_labels = c('1'='Cluster 1 (SP+/VGLUT1+)', '2'='Cluster 2 (TH-mid)',
    '4'='Cluster 4 (VACHT+)', '5'='Cluster 5 (TH-high)',
    '6'='Cluster 6 (SYP-high)', '7'='Cluster 7 (Indeterminate)'),
  # Output controls only; neither setting changes measurement calculations.
  export_geometry_tiffs = TRUE,
  export_mask_previews = TRUE
)
