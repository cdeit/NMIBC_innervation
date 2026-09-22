#!/usr/bin/env Rscript
# Run from any working directory. A new output directory is mandatory.
local({
  args<-commandArgs(trailingOnly=TRUE)
  if(!length(args)%in%c(3L,4L,5L))stop('Usage: Rscript R/compute_spatial_measurements.R MASK_DIR IMAGE_IDS.txt NEW_OUTPUT_DIR [spatial_measurement_config.R] [ANNOTATIONS.csv]')
  script<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])
  root<-normalizePath(file.path(dirname(script),'..'))
  for(f in c('spatial_io','runtime','output_schema','mask_io','preprocess_segmentation','debris','geometry','band_metrics','nerve_objects','vessel_objects','border_distances','vessel_morphology','process_image','format_outputs','export_masks'))source(file.path(root,'R',paste0(f,'.R')),local=TRUE)
  runtime_ok<-check_runtime(Sys.getenv('SPATIAL_ALLOW_UNVALIDATED')=='1')
  suppressPackageStartupMessages({library(EBImage);library(mmand);library(tiff);library(dplyr)})
  config_path<-if(length(args)>=4L)args[4] else file.path(root,'config/spatial_measurement_config.R')
  source(config_path,local=TRUE)
  ids<-read_image_ids(args[2])
  annotation_path<-if(length(args)==5L)args[5] else file.path(root,'data','metadata','spatial','cohort_annotations.csv')
  annotations<-read_annotations(annotation_path,ids)
  cfg$in_dir<-normalizePath(args[1],mustWork=TRUE)
  if(file.exists(args[3]))stop('Output already exists. Choose a NEW directory; existing runs are never overwritten.')
  # Validate every input before creating output. Measurement functions read only these masks.
  input_paths<-unlist(lapply(ids,function(id){p<-layer_paths(id,cfg);unlist(p[c('tumor','stroma','glass','vessel','phenotyped')])}))
  if(any(!file.exists(input_paths)))stop('Missing inputs: ',paste(input_paths[!file.exists(input_paths)],collapse='\n'))
  dir.create(args[3],recursive=TRUE);cfg$out_dir<-normalizePath(args[3])
  for(d in c('checkpoints','provenance','masks','geometry'))dir.create(file.path(cfg$out_dir,d))
  writeLines(ids,file.path(cfg$out_dir,'image_ids.txt'))
  file.copy(annotation_path,file.path(cfg$out_dir,'provenance','cohort_annotations.csv'))
  saveRDS(cfg,file.path(cfg$out_dir,'provenance','resolved_config.rds'))
  file.copy(config_path,file.path(cfg$out_dir,'provenance','spatial_measurement_config.R'))
  record_runtime(file.path(cfg$out_dir,'provenance'))
  code<-c(list.files(file.path(root,'R'),pattern='[.]R$',full.names=TRUE),list.files(file.path(root,'scripts'),pattern='[.](R|sh)$',full.names=TRUE),normalizePath(config_path))
  # Keep an executable source snapshot alongside the effective config.
  for(path in code){
    relative<-substring(path,nchar(root)+2L)
    target<-file.path(cfg$out_dir,'provenance','source',relative)
    dir.create(dirname(target),recursive=TRUE,showWarnings=FALSE)
    stopifnot(file.copy(path,target))
  }
  files<-c(input_paths,normalizePath(args[2]),normalizePath(annotation_path),normalizePath(config_path),code)
  write.csv(data.frame(file=files,md5=checked_md5(files)),file.path(cfg$out_dir,'provenance','input_and_code_md5.csv'),row.names=FALSE)
  writeLines(if(runtime_ok)'Recorded Alpine band runtime checks passed; output validation still required.' else 'UNVALIDATED DEVELOPMENT RUNTIME',file.path(cfg$out_dir,'provenance','runtime_status.txt'))
  for(id in ids){
    message('Processing ',id)
    result<-process_one_image(id,cfg)
    saveRDS(result,file.path(cfg$out_dir,'checkpoints',paste0(id,'.rds')))
    rm(result);gc(verbose=FALSE)
  }
  measured<-setNames(lapply(c('band_metrics_border','nerve_object_distances','vessel_object_distances'),function(nm){
    dplyr::bind_rows(lapply(ids,function(id)readRDS(file.path(cfg$out_dir,'checkpoints',paste0(id,'.rds')))[[nm]]))
  }),c('band_metrics_border','nerve_object_distances','vessel_object_distances'))
  formatted<-format_goal_tables(measured,annotations,cfg)
  for(nm in names(formatted))write_table(formatted[[nm]],nm,cfg$out_dir)
  writeLines(c(paste('Execution completed for',length(ids),'images.'),'Measurements completed; reference comparison is performed separately by R/validate_outputs.R.'),file.path(cfg$out_dir,'COMPLETE.txt'))
})
