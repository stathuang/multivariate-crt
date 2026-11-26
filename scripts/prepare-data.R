# Explicit import or preparation; input directories are never written.
args<-commandArgs(TRUE)
stopifnot(length(args)==2L,args[1L]%in%c('import','raw'))
script<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1L])
root<-normalizePath(file.path(dirname(script),'..'))
input<-normalizePath(args[2L],mustWork=TRUE)
local<-file.path(root,'data/local');dir.create(local,recursive=TRUE,showWarnings=FALSE)
source(file.path(root,'experiments/readable-experiments/R/run_tools.R'))
source(file.path(root,'experiments/readable-experiments/R/comparators.R'))
source(file.path(root,'experiments/readable-experiments/R/data.R'))
source(file.path(root,'experiments/readable-experiments/R/stock_data.R'))
source(file.path(root,'experiments/readable-experiments/R/bcr_data.R'))
RNGkind('Mersenne-Twister','Inversion','Rejection')
if(args[1L]=='raw') {
  needed<-c('Stock/DataForExp.RData','Stock/Pre-X.RData','Stock/list_stocks.csv',
            'Stock/DIA.json','Breast/GSE2034_series_matrix.txt','Breast/GSE2034-PI.csv','Breast/ImpVar.csv')
  stopifnot(all(file.exists(file.path(input,needed))))
  expected<-jsonlite::fromJSON(file.path(root,'data/required-inputs.json'))
  for(i in seq_len(nrow(expected))) {
    actual<-digest::digest(file=file.path(input,expected$file[i]),algo='sha256',serialize=FALSE)
    if(actual!=expected$sha256[i])stop('Input snapshot hash differs: ',expected$file[i])
  }
  stock<-prepare_stock(file.path(input,'Stock'))
  # Replace the incomplete historical final DIA window using the full 461-date snapshot.
  chart<-jsonlite::fromJSON(file.path(input,'Stock/DIA.json'),simplifyVector=FALSE)$chart$result[[1L]]
  prices<-as.numeric(unlist(chart$indicators$adjclose[[1L]]$adjclose))
  dates<-as.Date(as.POSIXct(as.numeric(unlist(chart$timestamp)),origin='1970-01-01',tz='UTC'))
  stopifnot(length(prices)==461L,length(dates)==461L,!anyDuplicated(dates),all(is.finite(prices)),all(prices>0))
  returns<-exp(diff(log(prices)))-1
  groups<-(seq_along(returns)-1L)%/%5L+1L
  stopifnot(all(tabulate(groups)==5L))
  stock$y<-as.numeric(tapply(returns,groups,mean))*100
  stock$nulls<-lapply(stock$targets,function(target)with_seed(2023L,{
    outside<-setdiff(seq_len(ncol(stock$X)),target)
    fit<-glmnet::cv.glmnet(stock$X[,outside,drop=FALSE],stock$y)
    predicted<-predict(fit,newx=stock$X[,outside,drop=FALSE],s=fit$lambda.min)
    list(coef=as.numeric(coef(fit,s=fit$lambda.min)),sd=sd(stock$y-predicted),lambda=fit$lambda.min,outside=outside,seed=2023L)
  }))
  stock$starting_point<-'Original 103-stock matrix; corrected 92 complete five-day DIA windows.'
  bcr<-prepare_bcr(file.path(input,'Breast'),selection_seed=2023L,null_seed=2023L)
} else {
  stock<-readRDS(file.path(input,'stock.rds'))
  bcr<-readRDS(file.path(input,'bcr_erplus60.rds'))
}
stopifnot(identical(dim(stock$X),c(92L,103L)),length(stock$y)==92L,
          identical(dim(bcr$X),c(209L,21L)),length(bcr$target)==18L,length(bcr$condition_probes)==3L,
          identical(bcr$panel,'wang_erplus60'))
objects<-list(stock.rds=stock,bcr_erplus60.rds=bcr)
destinations<-list(
  stock.rds=c(file.path(local,'stock.rds'),file.path(root,'experiments/readable-experiments/data/prepared/stock.rds'),
              file.path(root,'experiments/stock-rf1000/engine/data/external-cache/prepared/stock.rds')),
  bcr_erplus60.rds=c(file.path(local,'bcr_erplus60.rds'),file.path(root,'experiments/readable-experiments/data/prepared/bcr_erplus60.rds')))
for(name in names(objects)) for(path in destinations[[name]]) {
  if(file.exists(path)) stop('Refusing to overwrite: ',path)
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE);saveRDS(objects[[name]],path,version=3L)
}
manifest<-data.frame(name=names(objects),input_mode=args[1L],sha256=vapply(names(objects),function(n)
  digest::digest(file=file.path(local,n),algo='sha256',serialize=FALSE),''))
write.csv(manifest,file.path(local,'IMPORT.csv'),row.names=FALSE)
cat('Prepared stock (92 x 103) and ER-positive BCR (209 x 21; 3 background, 18 targets).\n')
