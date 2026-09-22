# Calculate only measurements needed by the three updated-format goal CSVs.
# Original TIFFs are read-only. No references or archived geometry enter here.
process_one_image <- function(image_id, cfg) {
  p<-layer_paths(image_id,cfg)
  tumor<-enr_read_binary(p$tumor);stroma<-enr_read_binary(p$stroma)
  glass<-enr_read_binary(p$glass);vessel<-enr_read_binary(p$vessel)
  phenotyped<-enr_read_coded(p$phenotyped,c(0L,cfg$nerve_codes))
  dims<-lapply(list(tumor,stroma,glass,vessel,phenotyped),dim)
  stopifnot(all(vapply(dims,identical,logical(1),dims[[1]])))

  # Apply tissue cleanup in memory, preserving the original FFT/contact rules.
  tc<-remove_debris(tumor,glass,cfg$max_glass_contact_frac,cfg$min_tissue_contact_frac,cfg$min_tumor_area_px)
  sc<-remove_debris(stroma,glass,cfg$max_glass_contact_frac,cfg$min_tissue_contact_frac,cfg$min_stroma_area_px)
  valid_stroma<-sc & !glass
  ts<-compute_tumor_stroma_border(tc,sc,glass,cfg$um_per_px,cfg$border_dilate_px)$border
  sg<-compute_stroma_glass_border(sc,glass,cfg$um_per_px,cfg$border_dilate_px)$border
  dt<-fast_distmap_um(ts,cfg$um_per_px,cfg$distmap_downsample_factor)

  # Signal and denominator include valid cleaned stroma. Object cleanup below
  # does NOT modify the coded nerve mask used for band-area measurement.
  bands<-compute_band_metrics(dt,cfg$border_breaks,cfg$border_labels,
    valid_stroma,valid_stroma,list(nerve=phenotyped,vessel=vessel),
    list(nerve='coded',vessel='binary'),list(nerve=cfg$nerve_codes),
    image_id,'tumor_border',cfg$um_per_px)

  # Derive subtype masks from their exact coded-mask partition, then label
  # nerves/vessels using the original, distinct object-cleanup rules.
  nerve<-label_nerve_objects(image_id,cfg,glass,phenotyped)
  ves<-label_vessel_objects(image_id,cfg)
  nd<-nerve$objects;vd<-ves$objects
  if(nrow(nd)>0L){
    nd<-dplyr::left_join(nd,object_to_border_distances_suffixed(
      nerve$boundary_points,'nerve_object_id',ts,'tumor_stroma_border',cfg$um_per_px,cfg),by='nerve_object_id')
    nd<-dplyr::left_join(nd,object_to_border_distances_suffixed(
      nerve$boundary_points,'nerve_object_id',sg,'stroma_glass_border',cfg$um_per_px,cfg),by='nerve_object_id')
  }
  # Updated vessel goal includes morphology and stroma-glass distances.
  if(nrow(vd)>0L){
    vd<-dplyr::left_join(vd,object_to_border_distances_suffixed(
      ves$boundary_points,'vessel_object_id',ts,'tumor_stroma_border',cfg$um_per_px,cfg),by='vessel_object_id')
    vd<-dplyr::left_join(vd,object_to_border_distances_suffixed(
      ves$boundary_points,'vessel_object_id',sg,'stroma_glass_border',cfg$um_per_px,cfg),by='vessel_object_id')
    vd<-dplyr::left_join(vd,compute_vessel_shape_features(ves,image_id),by='vessel_object_id')
  }
  export_geometry(image_id,cfg,tc,valid_stroma,ts,sg,vessel,nerve,ves)
  list(band_metrics_border=bands,nerve_object_distances=nd,vessel_object_distances=vd)
}
