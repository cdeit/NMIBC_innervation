# Exact updated-table measurement definitions

- Calibration is 0.5 um/pixel; one pixel is 0.25 um2.
- Tissue cleanup labels connected components, identifies boundary pixels by
  size-3 disc erosion and accumulates eight-neighbor glass/nonglass contacts.
  Nonglass contact includes same-class tissue. Remove components for glass
  fraction >=0.9, tissue fraction <0.1, complete coverage by size-7 dilated glass,
  or area <25,000 pixels AND computed glass contact >0. The historical FFT
  arithmetic and strict positivity test are preserved without rounding.
- Tumor–stroma border = cleaned tumor AND size-3 dilated cleaned stroma AND NOT
  glass. `makeBrush(3)` specifies size, not radius. Stroma–glass uses the same
  border builder with the corresponding compartments.
- Distance rasters sample every second row/column, starting at 1, apply distmap
  to the complement, scale by 1 um and expand using rounded index sequences.
  These are approximate raster distances. No correction for a segment lost
  during downsampling is inserted into the historical calculation.
- Band denominators use valid cleaned stroma INCLUDING vessel-occupied pixels.
  Signal areas intersect it with the original coded nerve or filtered vessel
  mask. Object-specific nerve cleanup does not alter band signal areas.
- Intervals are [0,25], (25,50], (50,100], (100,200], (200,300], (300,Inf] um.
  Target labels, including `201-300 um`, are preserved as display labels only.
  `proportion_object_area_per_band` = signal area / stromal band area; zero
  denominator gives NA. No pseudocount or separate enrichment-ratio calculation
  is required for these three target tables.
- Nerves are labeled by subtype after glass/contact cleanup with a conditional
  five-pixel floor. Vessel components are labeled from the prefiltered mask with
  an unconditional 30-pixel floor. Contours use ocontour with +1 to convert its
  zero-based coordinates to R matrix indices. Centroids are descriptive only.
- Per-object border distances take the minimum raster value at contour points.
  Both nearest and longest connected-segment results and segment lengths are
  saved for tumor–stroma and stroma–glass borders. Length is pixel count * 0.5 um.
- Vessel shape features retain the original EBImage shape/moment and mmand
  closing/skeletonization calculations. Hull-related features use the original
  morphological closing approximations, not substituted geometric convex hulls.
- Final labels, StudyID and inclusion/exclusion annotations are joined from the
  supplied annotation file. No measurements are imported from goal tables.
  Clinical/manual inclusion is not inferable from pixels. All requested image
  rows remain in output regardless of their inclusion flags, as in the targets.
- No 300 um distance cutoff is applied to object-distance tables.
