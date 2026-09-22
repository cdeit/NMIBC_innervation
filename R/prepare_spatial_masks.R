#!/usr/bin/env Rscript
# Input is already segmented TIFFs, not unsegmented fluorescence microscopy.
local({
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=5L)stop('Usage: Rscript R/prepare_spatial_masks.R RAW_TISSUE_NERVE_DIR RAW_CD31_DIR PHENOTYPE_CSV IMAGE_IDS.txt NEW_MASK_DIR')
script<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]);root<-normalizePath(file.path(dirname(script),'..'))
source(file.path(root,'R/spatial_io.R'),local=TRUE)
source(file.path(root,'R/runtime.R'),local=TRUE)
check_runtime(Sys.getenv('SPATIAL_ALLOW_UNVALIDATED')=='1')
require_packages(c('EBImage','mmand','tiff','dplyr','tibble','tidyr','fftwtools'))
suppressPackageStartupMessages({library(EBImage);library(mmand);library(tiff);library(dplyr)})
source(file.path(root,'R/preprocess_segmentation.R'),local=TRUE)
ids<-read_image_ids(args[4]);dest<-args[5];if(file.exists(dest))stop('Choose a new output directory; existing masks will not be overwritten.')
p<-read.csv(args[3],stringsAsFactors=FALSE);stopifnot(all(c('image_id','objID','phenotype')%in%names(p)))
p<-p[p$image_id%in%ids,];if(anyDuplicated(paste(p$image_id,p$objID)))stop('Duplicate image/object phenotype assignments.')
if(anyNA(p$objID)||any(p$objID<1|p$objID!=as.integer(p$objID)))stop('objID must contain positive integer raw-instance IDs.')
codes<-c(1L,2L,4L,5L,6L,7L)
if(any(!p$phenotype%in%paste0('subcluster_',codes)))stop('Supply the filtered nerve phenotype table with only subclusters 1,2,4,5,6,7.')
paths<-c(file.path(args[1],paste0(ids,'_binary_seg_maps.tif')),file.path(args[2],paste0(ids,'_binary_seg_maps.tif')))
if(any(!file.exists(paths)))stop('Missing raw segmentation files. Expected normalized image IDs in filenames.')
dir.create(dest,recursive=TRUE);record_runtime(file.path(dest,'provenance'))
writeLines(capture.output(sessionInfo()),file.path(dest,'preprocessing_sessionInfo.txt'))
write.csv(data.frame(file=c(paths,args[3],args[4]),md5=checked_md5(c(paths,args[3],args[4]))),file.path(dest,'preprocessing_input_checksums.csv'),row.names=FALSE)
for(id in ids){
 message('Preparing ',id)
 pages<-tiff::readTIFF(file.path(args[1],paste0(id,'_binary_seg_maps.tif')),all=TRUE,as.is=TRUE)
 if(length(pages)<2L)stop('Raw tissue/nerve TIFF needs semantic page 1 and instance page 2.')
 sem<-pages[[1]];inst<-pages[[2]];storage.mode(sem)<-'integer';storage.mode(inst)<-'integer'
 cd31<-read_tif_layer(file.path(args[2],paste0(id,'_binary_seg_maps.tif')),2L)
 stopifnot(identical(dim(sem),dim(inst)),identical(dim(sem),dim(cd31)),all(sem%in%0:4),all(inst>=0))
 ph<-p[p$image_id==id,];if(any(!ph$objID%in%inst))stop('Phenotype objID missing from raw instance layer: ',id)
 nerve<-map_nerve_phenotypes(inst,ph,codes)$mask
 vessel_result<-filter_vessels(cd31,id)
 if(any(cd31>0)&&is.null(vessel_result$feats))stop("Vessel feature extraction failed: ",id)
 vessel<-vessel_result$mask
 masks<-list(tumor=1L*(sem==0L),stroma=1L*(sem==1L),glass=1L*(sem==2L),nerve_phenotyped=nerve,vessel_filtered=vessel,vessel_unfiltered=1L*(cd31>0))
 # Export the original per-subtype single-channel masks as well as coded nerve.
 subtype_names<-c('1'='nerve_SP_pos_VGLUT1_pos_sc1','2'='nerve_TH_mid_sc2','4'='nerve_VACHT_pos_sc4','5'='nerve_TH_high_sc5','6'='nerve_SYP_high_sc6','7'='nerve_TH_neg_VACHT_neg_SP_neg_VGLUT1_neg_sc7')
 for(code in codes)masks[[subtype_names[as.character(code)]]]<-1L*(nerve==code)
 folder<-file.path(dest,id);dir.create(folder)
 for(nm in names(masks))write_label_tif(masks[[nm]],file.path(folder,paste0(nm,'.tif')))
}
writeLines(paste('Prepared',length(ids),'images. Phenotypes were mapped by image_id + objID, not table row order.'),file.path(dest,'COMPLETE.txt'))
})
