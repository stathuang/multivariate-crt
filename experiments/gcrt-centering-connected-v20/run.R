# Fixed 2500-replicate run. No adaptive stopping or scenario changes.
options(warn=2,digits=17)
stopifnot(basename(normalizePath('.'))=='gcrt-centering-connected-v20',
          !file.exists('run/COMPLETE.json'))
cfg <- jsonlite::fromJSON('run/CONFIG.json')
stopifnot(cfg$replicates==2500, cfg$n==80,cfg$t==8,cfg$s==112,cfg$M==99,
  jsonlite::fromJSON('run/IMPLEMENTATION_CHECK.json')$passed)
RNGkind(cfg$rng[1],cfg$rng[2],cfg$rng[3])
source('core.R')
dir.create('run/records',showWarnings=FALSE)
capture.output(list(time=as.character(Sys.time()),session=sessionInfo(),RNG=RNGkind(),
  threads=Sys.getenv(c('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS'))),
  file='run/environment.txt')
tests <- list(); checks <- list(); start <- proc.time()[['elapsed']]
for (batch in 1:25) {
  path <- sprintf('run/records/batch-%02d.rds',batch)
  stopifnot(!file.exists(path))
  records <- list()
  for (k in ((batch-1)*100+1):(batch*100)) {
    ans <- tryCatch(one_replicate(k),error=function(e)e)
    if (inherits(ans,'error')) {
      write.csv(data.frame(replicate=k,error=conditionMessage(ans)), 'run/failures.csv',row.names=FALSE)
      saveRDS(records,sprintf('run/records/failed-batch-%02d.rds',batch))
      stop('Production stopped; failure retained; no denominator reduction.')
    }
    records[[as.character(k)]] <- ans
    tests[[k]] <- test_rows(ans)
    checks[[k]] <- data.frame(replicate=k,data_seed=ans$data_seed,sampler_seed=ans$sampler_seed,
      invariant_error=ans$invariant_error,min_pool_eigenvalue=ans$min_pool_eigenvalue,seconds=ans$seconds)
  }
  saveRDS(records,path,compress=TRUE)
  cat(sprintf('%s: %d / 2500 independent datasets completed; elapsed %.1f seconds\n',
    format(Sys.time(),'%H:%M:%S'),batch*100,proc.time()[['elapsed']]-start));flush.console()
}
tests <- do.call(rbind,tests); checks <- do.call(rbind,checks)
stopifnot(nrow(tests)==60000,nrow(checks)==2500,
  !anyDuplicated(tests[c('replicate','case','method')]),all(table(tests$case,tests$method)==2500))
write.csv(tests,'run/tests.csv',row.names=FALSE);write.csv(checks,'run/checks.csv',row.names=FALSE)
write.csv(data.frame(replicate=integer(),error=character()),'run/failures.csv',row.names=FALSE)
summary <- do.call(rbind,lapply(cells$case,function(ca) do.call(rbind,lapply(methods,function(me) {
  d <- tests[tests$case==ca & tests$method==me,]; ci <- binom.test(sum(d$reject),nrow(d))$conf.int
  data.frame(case=ca,b=d$b[1],r=d$r[1],method=me,replicates=nrow(d),rejects=sum(d$reject),
    rate=mean(d$reject),mcse=sqrt(mean(d$reject)*(1-mean(d$reject))/nrow(d)),ci_lower=ci[1],ci_upper=ci[2])
}))))
comparisons <- list(c('centered_pool','original'),c('centered_identity','original'),
                    c('centered_pool','centered_identity'))
paired <- do.call(rbind,lapply(cells$case,function(ca) do.call(rbind,lapply(comparisons,function(pair) {
  a <- tests[tests$case==ca & tests$method==pair[1],]; b <- tests[tests$case==ca & tests$method==pair[2],]
  z <- merge(a,b,by='replicate',suffixes=c('_a','_b')); d <- as.numeric(z$reject_a)-as.numeric(z$reject_b)
  diff <- mean(d); se <- sd(d)/sqrt(length(d)); half <- qnorm(.975)*se
  data.frame(case=ca,b=a$b[1],r=a$r[1],method_a=pair[1],method_b=pair[2],pairs=length(d),
    difference=diff,paired_mcse=se,ci_lower=diff-half,ci_upper=diff+half,
    a_only=sum(d==1),b_only=sum(d==-1))
}))))
write.csv(summary,'run/summary.csv',row.names=FALSE);write.csv(paired,'run/paired_differences.csv',row.names=FALSE)
jsonlite::write_json(list(completed_at=as.character(Sys.time()),independent_datasets=2500,
  pvalues=nrow(tests),summary_rows=nrow(summary),paired_rows=nrow(paired),failures=0,
  elapsed_seconds=proc.time()[['elapsed']]-start,max_invariant_error=max(checks$invariant_error),
  minimum_pooled_eigenvalue=min(checks$min_pool_eigenvalue)),
  'run/COMPLETE.json',pretty=TRUE,auto_unbox=TRUE,digits=17)
cat('All prespecified datasets and cells completed.\n')
