# Public orchestration helpers; computational kernels are preserved separately.
read_image_ids <- function(path) {
  ids <- trimws(readLines(path, warn=FALSE)); ids <- ids[nzchar(ids) & !startsWith(ids,'#')]
  if(!length(ids) || anyDuplicated(ids))stop('Image manifest must be nonempty and unique.')
  if(any(grepl('[/\\\\]',ids)) || any(ids %in% c('.','..')))stop('Image IDs must be directory names, not paths.')
  ids
}
require_packages <- function(packages) {
  missing<-packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
  if(length(missing))stop('Missing packages: ',paste(missing,collapse=', '),'. See environment/README.md.')
}
capture_environment <- function(dest) {
  writeLines(capture.output(sessionInfo()),file.path(dest,'sessionInfo.txt'))
  pkgs<-c('EBImage','fftwtools','tiff','dplyr','tibble','tidyr')
  write.csv(data.frame(package=pkgs,version=vapply(pkgs,function(p)as.character(packageVersion(p)),character(1))),file.path(dest,'packages.csv'),row.names=FALSE)
  if(Sys.info()[['sysname']]=='Linux'){
    if(file.exists('/etc/os-release'))file.copy('/etc/os-release',file.path(dest,'os_release.txt'))
    for(cmd in c('lscpu','ldd'))if(nzchar(Sys.which(cmd))){
      args<-if(cmd=='ldd')shQuote(system.file('libs','fftwtools.so',package='fftwtools')) else character()
      writeLines(tryCatch(system2(cmd,args,stdout=TRUE,stderr=TRUE),error=function(e)conditionMessage(e)),file.path(dest,paste0(cmd,'.txt')))
    }
  }
}
