# Annotation joins and exact updated-format column order.
# Annotations are supplied clinical/manual decisions, not inferred measurements.
read_annotations <- function(path, ids) {
  x<-read.csv(path,stringsAsFactors=FALSE,check.names=FALSE)
  fields<-c('image_id','StudyID','keep_or_exclude','include_image_final_spatial_cohort','reason_exclude_final_spatial_cohort')
  if(!all(fields%in%names(x))||anyDuplicated(x$image_id)||any(!ids%in%x$image_id))stop('Annotations require a unique complete row for every requested image.')
  x<-x[match(ids,x$image_id),fields,drop=FALSE]
  if(!is.logical(x$include_image_final_spatial_cohort)||anyNA(x$include_image_final_spatial_cohort))stop('Final cohort flags must be TRUE/FALSE.')
  if(anyNA(x$StudyID)||any(!x$keep_or_exclude%in%c('keep','omit')))stop('Invalid StudyID or keep_or_exclude annotation.')
  x
}
format_goal_tables <- function(measured, annotations, cfg) {
  bands<-measured$band_metrics_border
  bands$band_label<-cfg$border_labels[bands$band_index]
  bands$subtype_label<-ifelse(bands$layer_name=='vessel','Vessel (CD31+)',unname(cfg$nerve_labels[as.character(bands$subtype_code)]))
  names(bands)[names(bands)=='band_density']<-'proportion_object_area_per_band'
  flags<-c('image_id','include_image_final_spatial_cohort','reason_exclude_final_spatial_cohort')
  bands<-dplyr::left_join(bands,annotations[,flags],by='image_id')
  out<-list(nerve_vessel_area_by_distance_to_tumor_stroma_border=bands)
  for(kind in c('nerve','vessel')){
    target<-paste0(kind,'_object_distances_to_tumor_stroma_border')
    x<-measured[[paste0(kind,'_object_distances')]]
    if(nrow(x)==0L){out[[target]]<-output_schema[[target]];next}
    x$subtype_label<-if(kind=='vessel')rep('Vessel (CD31+)',nrow(x)) else unname(cfg$nerve_labels[as.character(x$subtype_code)])
    x<-dplyr::left_join(x,annotations,by='image_id')
    out[[target]]<-x
  }
  for(nm in names(out))out[[nm]]<-conform_table(out[[nm]],nm)
  out[names(output_schema)]
}
