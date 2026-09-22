# Export the actual geometry used by the measurements. Signed TIFFs preserve
# the input convention; PNG previews use ordinary foreground=white encoding.
export_geometry <- function(image_id, cfg, tumor, stroma, ts, sg, vessel, nerve, ves) {
  folder<-file.path(cfg$out_dir,'masks',image_id)
  dir.create(folder,recursive=TRUE,showWarnings=FALSE)
  coded_clean<-matrix(0L,nrow(tumor),ncol(tumor))
  for(code in names(nerve$filled_masks))coded_clean[nerve$filled_masks[[code]]]<-as.integer(code)
  vessel_retained<-if(is.null(ves$lab_mat))matrix(FALSE,nrow(tumor),ncol(tumor)) else ves$lab_mat>0
  masks<-list(tumor_clean=tumor,stroma_clean=stroma,tumor_stroma_border=ts,
    stroma_glass_border=sg,placement_valid=stroma & !vessel,
    nerve_phenotyped_clean=coded_clean,vessel_objects_retained=vessel_retained)
  if(isTRUE(cfg$export_geometry_tiffs))for(nm in names(masks)){
    write_label_tif(masks[[nm]],file.path(folder,paste0(nm,'.tif')))
  }
  if(isTRUE(cfg$export_mask_previews))for(nm in names(masks)){
    # Coded nerve preview shows the combined footprint; subtype codes are in TIFF.
    m<-masks[[nm]]>0
    png::writePNG(matrix(as.numeric(m),nrow(m),ncol(m)),file.path(folder,paste0(nm,'.png')))
  }
  saveRDS(list(tumor_clean=tumor,stroma_clean=stroma,border=ts,
    stroma_glass_border=sg,placement=stroma & !vessel),
    file.path(cfg$out_dir,'geometry',paste0(image_id,'_geometry.rds')))
}
