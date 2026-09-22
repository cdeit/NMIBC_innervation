#!/usr/bin/env Rscript
# Goals are read ONLY by this validator, never by the measurement runner.
local({
 args<-commandArgs(trailingOnly=TRUE)
 if(length(args)!=2L)stop('Usage: Rscript R/validate_outputs.R RUN_DIR UPDATED_GOAL_DIR')
 script<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]);root<-normalizePath(file.path(dirname(script),'..'))
 source(file.path(root,'R/output_schema.R'),local=TRUE)
 out<-normalizePath(args[1]);ref<-normalizePath(args[2]);ids<-readLines(file.path(out,'image_ids.txt'))
 dest<-file.path(out,'validation');dir.create(dest,showWarnings=FALSE)
 marker<-file.path(dest,'REFERENCE_MATCH.txt');if(file.exists(marker))unlink(marker)
 tables<-names(output_schema)
 report<-list();summary<-list();failed<-FALSE
 for(nm in tables){
  ap<-file.path(ref,paste0(nm,'.csv'));bp<-file.path(out,paste0(nm,'.csv'))
  if(!file.exists(ap)||!file.exists(bp))stop('Missing required comparison table: ',nm)
  a<-read.csv(ap,check.names=FALSE);b<-read.csv(bp,check.names=FALSE)
  a<-a[a$image_id%in%ids,,drop=FALSE]
  fields<-if(startsWith(nm,'nerve_vessel_area'))c('image_id','reference_name','layer_name','subtype_code','band_index') else if(startsWith(nm,'nerve'))'nerve_object_id' else 'vessel_object_id'
  key<-function(x)do.call(paste,c(x[fields],sep='|'))
  ka<-key(a);kb<-key(b)
  keys_ok<-!anyDuplicated(ka)&&!anyDuplicated(kb)&&setequal(ka,kb)
  cols_ok<-identical(names(b),names(output_schema[[nm]]))&&identical(names(a),names(b))
  if(cols_ok)a<-a[,names(b),drop=FALSE]
  summary[[nm]]<-data.frame(table=nm,reference_rows=nrow(a),candidate_rows=nrow(b),keys_match=keys_ok,columns_match=cols_ok,missing_keys=length(setdiff(ka,kb)),extra_keys=length(setdiff(kb,ka)))
  if(!keys_ok||!cols_ok){failed<-TRUE;next}
  b<-b[match(ka,kb),,drop=FALSE]
  for(col in names(a)){
   x<-a[[col]];y<-b[[col]];same<-is.na(x)&is.na(y);ok<-!is.na(x)&!is.na(y)
   same[ok]<-x[ok]==y[ok]
   # Integer counts, areas, IDs, labels and missingness are exact.
   exact<-grepl('(^|_)(area_px|area_um2|s.area|branch_count)$',col)||col%in%c('stroma_area_um2','radius_area','radius_denom','loo_area','loo_denom','outer_denom','local_id','object_id','subtype_code','band_index')
   finite<-ok&is.finite(if(is.numeric(x))x else rep(NA_real_,length(x)))&is.finite(if(is.numeric(y))y else rep(NA_real_,length(y)))
   if(is.numeric(x)&&is.numeric(y)&&!exact)same[finite]<-abs(x[finite]-y[finite])<=1e-8+1e-10*abs(x[finite])
   delta<-if(is.numeric(x)&&is.numeric(y)&&any(finite))max(abs(x[finite]-y[finite])) else NA_real_
   report[[length(report)+1]]<-data.frame(table=nm,column=col,differing_cells=sum(!same),missingness_changes=sum(is.na(x)!=is.na(y)),max_abs_difference=delta,rule=if(exact)'exact' else 'numeric tolerance / categorical equality')
   if(any(!same))failed<-TRUE
  }
 }
 write.csv(do.call(rbind,summary),file.path(dest,'table_summary.csv'),row.names=FALSE)
 if(length(report))write.csv(do.call(rbind,report),file.path(dest,'column_comparison.csv'),row.names=FALSE)
 if(failed)stop('Reference comparison failed. Inspect validation/table_summary.csv and column_comparison.csv.')
 writeLines(c(paste('PASS for',length(ids),'requested images across all three updated-format goal tables.'),'Exact keys, areas/counts and missingness; other numeric values within 1e-8 + 1e-10 * abs(reference).','Row order and CSV serialization are not compared.'),marker)
 cat(readLines(marker),sep='\n')
})
