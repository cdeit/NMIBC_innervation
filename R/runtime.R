# Runtime and file-contract checks. These do not alter measurement parameters.
publication_packages <- c('EBImage','fftwtools','tiff','png','mmand','dplyr','tidyr','tibble')
check_runtime <- function(allow_unvalidated = FALSE) {
  require_packages(publication_packages)
  expected <- c(mmand='1.6.3',EBImage='4.48.0',fftwtools='0.9-11',tiff='0.1-12',dplyr='1.1.4',tidyr='1.3.1',tibble='3.2.1')
  installed <- vapply(names(expected),function(p)packageDescription(p)$Version,character(1))
  # packageVersion() normalizes hyphens to dots; compare DESCRIPTION strings.
  matched <- installed == expected
  print(data.frame(package=names(expected),expected=unname(expected),installed=unname(installed),matches=matched),row.names=FALSE)
  fftw <- if(Sys.info()[['sysname']]=='Linux' && nzchar(Sys.which('dpkg-query'))) {
    suppressWarnings(system2('dpkg-query',c('-W','libfftw3-double3'),stdout=TRUE,stderr=TRUE))
  } else ''
  os <- if(file.exists('/etc/os-release'))readLines('/etc/os-release') else ''
  linkage <- if(Sys.info()[['sysname']]=='Linux' && nzchar(Sys.which('ldd'))) {
    system2('ldd',shQuote(system.file('libs','fftwtools.so',package='fftwtools')),stdout=TRUE,stderr=TRUE)
  } else ''
  ok <- all(matched) && as.character(getRversion())=='4.4.1' &&
    grepl('x86_64.*linux',R.version$platform) && any(grepl('3.3.8-2ubuntu8',fftw,fixed=TRUE)) &&
    any(grepl('Ubuntu 22.04.5 LTS',os,fixed=TRUE)) && any(grepl('libfftw3.so.3 =>',linkage,fixed=TRUE))
  if(!ok && !allow_unvalidated)stop('Runtime differs from the validated Alpine band environment. Use the Alpine launcher. SPATIAL_ALLOW_UNVALIDATED=1 is for development only.')
  if(!ok)warning('UNVALIDATED RUNTIME: successful execution does not imply historical agreement.')
  invisible(ok)
}
record_runtime <- function(dest) {
  dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  capture_environment(dest)
  write.csv(data.frame(package=publication_packages,version=vapply(publication_packages,function(p)packageDescription(p)$Version,character(1))),file.path(dest,'packages.csv'),row.names=FALSE)
  writeLines(capture.output(sessionInfo()),file.path(dest,'sessionInfo.txt'))
  binaries<-unlist(lapply(publication_packages,function(p)list.files(system.file('libs',package=p),pattern='\\.(so|dylib)$',full.names=TRUE)))
  if(length(binaries))write.csv(data.frame(file=binaries,md5=unname(tools::md5sum(binaries))),file.path(dest,'package_binary_md5.csv'),row.names=FALSE)
  if(Sys.info()[['sysname']]=='Linux'&&nzchar(Sys.which('dpkg-query')))writeLines(suppressWarnings(system2('dpkg-query',c('-W','libfftw3-double3'),stdout=TRUE,stderr=TRUE)),file.path(dest,'fftw_system_package.txt'))
}
# Hash files individually, retrying transient filesystem read failures once.
checked_md5 <- function(files) {
  out <- vapply(files,function(path) {
    value <- suppressWarnings(unname(tools::md5sum(path)))
    if(is.na(value))value <- unname(tools::md5sum(path))
    if(is.na(value))stop('Unable to checksum input/source file: ',path)
    value
  },character(1))
  unname(out)
}

# Empty-object images still produce properly headed CSVs.
conform_table <- function(x, name) {
  template<-output_schema[[name]]
  if(is.null(template))return(as.data.frame(x))
  if(nrow(x)>0L && !setequal(names(x),names(template)))stop('Unexpected schema for ',name)
  if(nrow(x)==0L)return(template)
  as.data.frame(x[,names(template),drop=FALSE])
}
write_table <- function(x,name,dest) {
  x<-conform_table(x,name)
  write.csv(x,file.path(dest,paste0(name,'.csv')),row.names=FALSE,na='NA')
  invisible(x)
}
