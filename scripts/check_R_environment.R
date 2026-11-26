# Explicit, offline dependency preparation. Never sourced by analysis entry points.
# Default mode only checks the installed library against the exported inventory.
args <- commandArgs(trailingOnly=TRUE)
value <- function(key,default=NULL) {
  i<-match(key,args);if(is.na(i))return(default)
  if(i==length(args)||startsWith(args[i+1L],'--'))stop('Missing value for ',key)
  args[i+1L]
}
filearg <- grep('^--file=',commandArgs(),value=TRUE)
root <- normalizePath(file.path(dirname(sub('^--file=','',filearg[1L])),'..'),mustWork=TRUE)
inventory <- file.path(root,'environment','R-dependencies.csv')
if(!file.exists(inventory))stop('Missing environment/R-dependencies.csv')
pins <- read.csv(inventory,stringsAsFactors=FALSE,check.names=FALSE)
if(!all(c('package','version','dependencies') %in% names(pins))||anyDuplicated(pins$package)) stop('Invalid dependency inventory')
lib <- value('--library')
if(!is.null(lib)) {
  if('--install' %in% args)dir.create(lib,recursive=TRUE,showWarnings=FALSE)
  lib<-normalizePath(lib,mustWork=TRUE);.libPaths(c(lib,.libPaths()))
}
installed_version <- function(p) {
  d<-tryCatch(utils::packageDescription(p),error=function(e)NULL)
  if(is.null(d)||!is.list(d)||is.null(d$Version))return(NA_character_)
  as.character(d$Version)
}
matches <- function(p,v) {actual<-installed_version(p);!is.na(actual)&&utils::compareVersion(actual,v)==0L}
if('--download' %in% args) {
  source_dir<-value('--source-dir');if(is.null(source_dir))stop('--download requires --source-dir')
  dir.create(source_dir,recursive=TRUE,showWarnings=FALSE)
  for(i in seq_len(nrow(pins))) {
    p<-pins$package[i];v<-pins$version[i];name<-paste0(p,'_',v,'.tar.gz')
    dest<-file.path(source_dir,name)
    if(file.exists(dest))next
    urls<-c(paste0('https://cran.r-project.org/src/contrib/',name),
            paste0('https://cran.r-project.org/src/contrib/Archive/',p,'/',name))
    if(p%in%c('BiocGenerics','graph','RBGL'))
      urls<-paste0('https://bioconductor.org/packages/3.17/bioc/src/contrib/',name)
    ok<-FALSE
    for(url in urls) {
      success<-tryCatch(suppressWarnings(download.file(url,dest,mode='wb',quiet=TRUE))==0L,error=function(e)FALSE)
      if(success){ok<-TRUE;break}
      if(file.exists(dest))unlink(dest)
    }
    if(!ok)stop('Pinned source unavailable: ',name,'. Obtain the original archive before continuing.')
  }
  cat('Pinned source archives downloaded. Review checksums and upstream licenses before installation.\n')
  if(!('--install'%in%args))quit(status=0L)
}
if('--install' %in% args) {
  if(is.null(lib))stop('--install requires an explicit --library directory')
  source_dir<-value('--source-dir');if(is.null(source_dir))stop('--install requires --source-dir containing manually obtained package_version.tar.gz files')
  source_dir<-normalizePath(source_dir,mustWork=TRUE)
  pending<-pins$package;done<-character()
  while(length(pending)) {
    ready<-pending[vapply(pending,function(p){d<-strsplit(pins$dependencies[match(p,pins$package)],';',fixed=TRUE)[[1L]]
      !length(intersect(d,pending))},TRUE)]
    if(!length(ready))stop('Dependency cycle or malformed inventory; no automatic fallback')
    for(p in ready) {
      v<-pins$version[match(p,pins$package)]
      if(!matches(p,v)) {
        archive<-file.path(source_dir,paste0(p,'_',v,'.tar.gz'))
        if(!file.exists(archive))stop('Missing pinned source archive: ',basename(archive),'; obtain it explicitly. Availability has not been verified.')
        utils::install.packages(archive,repos=NULL,type='source',lib=lib)
        if(!matches(p,v))stop('Installed version does not match requested pin: ',p,' ',v)
      }
    }
    done<-c(done,ready);pending<-setdiff(pending,ready)
  }
}
actual<-vapply(pins$package,installed_version,'')
ok<-mapply(matches,pins$package,pins$version)
print(data.frame(package=pins$package,pinned=pins$version,installed=actual,match=ok),row.names=FALSE)
if(!all(ok))quit(status=1L)
cat('All recorded dependency versions match. This is not an OS/toolchain or historical-environment guarantee.\n')
