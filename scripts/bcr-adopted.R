# Reproduce the adopted original RNG path, then apply the certified LRT and RF repair.
args <- commandArgs(TRUE)
stopifnot(length(args)>=2L,args[1L]%in%c('bcr-observed','bcr-null'))
root <- Sys.getenv('REPLICATION_ROOT'); stopifnot(nzchar(root))
out <- normalizePath(args[2L],mustWork=TRUE)
old <- new.env(parent=globalenv()); current <- new.env(parent=globalenv())
for(p in c('jsonlite','digest','MASS','Matrix','glmnet','randomForest','rags2ridges',
           'VGAM','lpSolve','gmp','Rmpfr','CVXR','ECOSolveR','nnet'))
  stopifnot(requireNamespace(p,quietly=TRUE))
for(f in c('core.R','statistics.R','baselines.R','applications.R'))
  sys.source(file.path(root,'experiments/bcr-original/R',f),old)
for(f in c('sampler.R','data.R','linear_tests.R','gaussian_tests.R','multinomial_likelihood.R','bcr_tests.R'))
  sys.source(file.path(root,'experiments/readable-experiments/R',f),current)
observed <- args[1L]=='bcr-observed'
manifest <- jsonlite::fromJSON(file.path(root,'reference/bcr',if(observed)'observed-manifest.json' else 'original-manifest.json'),simplifyVector=FALSE)
cfg <- manifest$provenance$config
for(k in c('statistics','baselines','methods','graph_ids')) cfg[[k]]<-unlist(cfg[[k]])
cfg$data_file <- 'bcr_erplus60.rds'
data_dir <- file.path(root,'data/local'); stopifnot(file.exists(file.path(data_dir,cfg$data_file)))
raw_gcrt <- old$gcrt_test
saved <- NULL
old$gcrt_test <- function(...) {
  a<-list(...)
  saved <<- list(arguments=a,rng=.Random.seed)
  do.call(raw_gcrt,a)
}
ids <- if(observed)1L else if('--trial'%in%args)c(1L,74L) else 1:400
results <- list()
RNGkind('Mersenne-Twister','Inversion','Rejection')
for(i in ids) {
  set.seed(if(observed)cfg$observed_seed else cfg$seed_base+i)
  legacy <- old$application_replicate(cfg,i,data_dir)
  stopifnot(!is.null(saved))
  # Legacy calls above consume exactly their original random stream, including old failures.
  a<-saved$arguments
  fit<-current$bcr_lrt_test(a[[1L]],a[[2L]],a[[3L]])
  repaired<-legacy
  k<-which(repaired$method=='LRT')
  repaired$pvalue[k]<-fit$pvalue;repaired$status[k]<-'ok';repaired$error[k]<-''
  k<-which(repaired$method=='RF')
  if(repaired$status[k]!='ok') {
    assign('.Random.seed',saved$rng,envir=.GlobalEnv)
    fixed<-current$gcrt_rf_test(a[[1L]],a[[2L]],a[[3L]],a[[4L]],'RF',cfg$p_rule,response_family='multinomial')
    repaired$pvalue[k]<-fixed$pvalue;repaired$status[k]<-'ok';repaired$error[k]<-''
  }
  stopifnot(all(repaired$status=='ok'),all(is.finite(repaired$pvalue)))
  expected<-read.csv(file.path(root,'reference/bcr',if(observed)'observed.csv' else 'UNRESTRICTED_RESULTS.csv'))
  expected<-expected[expected$replicate==i,]
  j<-match(repaired$method,expected$method);stopifnot(!anyNA(j))
  error<-max(abs(repaired$pvalue-expected$pvalue[j]))
  saveRDS(list(legacy=legacy,adopted=repaired,maximum_reference_difference=error,
               lrt_diagnostics=fit$diagnostics),file.path(out,sprintf('replicate-%05d.rds',i)))
  cat('Replication',i,'maximum adopted-reference p-value difference',error,'\n');flush.console()
  stopifnot(error<1e-10)
  results[[as.character(i)]]<-repaired
}
write.csv(do.call(rbind,results),file.path(out,'results.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
