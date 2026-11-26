# Run only from this directory: Rscript --vanilla run.R
# The fixed data generator and actual hub-and-spoke sampler are adapted from
# ../v19-direction-minimal-20261003-230434/run.R. That file is never sourced.
options(warn=2, digits=17)
Sys.setenv(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1')
stopifnot(basename(normalizePath(getwd()))=='gcrt-symmetric-centering-20261004',
          Sys.readlink(getwd())=='', !file.exists('records.rds'),
          requireNamespace('jsonlite', quietly=TRUE))
cfg <- jsonlite::fromJSON('CONFIG.json')
RNGkind(cfg$rng[1], cfg$rng[2], cfg$rng[3])
M <- as.integer(cfg$M); reps <- as.integer(cfg$repetitions_per_n)
stopifnot(M==19L,reps==100L,cfg$L==1L,
          cfg$background_penalty_coefficient==1,cfg$target_penalty==1,
          identical(as.numeric(cfg$sample_sizes),c(40,120,360)))
cases <- names(cfg$cases); methods <- cfg$methods
center <- function(A) sweep(A,2,colMeans(A),'-')
now <- function() proc.time()[['elapsed']]
rank_p <- function(x) (1+sum(x[-1L]>=x[1L]))/length(x)
center_scores <- function(C,G) {
  delta <- sweep(C,1,rowMeans(C),'-')
  colSums(delta*solve(G,delta))
}
relative_error <- function(x,y) max(abs(x-y))/max(1,max(abs(x)),max(abs(y)))
records <- list()
failures <- data.frame(n=integer(),replicate=integer(),seed=integer(),error=character())
started <- now()
capture.output(list(started=as.character(Sys.time()),system=Sys.info(),
  threads=Sys.getenv(c('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS')),
  RNG=RNGkind(),session=sessionInfo()),file='environment.txt')
one <- function(n,b,seed) {
  set.seed(seed); time0 <- now()
  # This block and its RNG order match the reference exactly.
  raw <- matrix(rnorm(4L*n),n,4L); centered <- center(raw)
  z <- centered[,1]; E <- centered[,2:3,drop=FALSE]; eps <- centered[,4]
  X <- E+z; zz <- sum(z*z); h <- n/(zz+n)
  H <- function(A) A-(1-h)*tcrossprod(z,as.vector(crossprod(z,A))/zz)
  Yold <- sapply(c(null=0,positive=.25,negative=-.25),function(beta) z+beta*rowSums(E)+eps)
  Y <- cbind(Yold,orthogonal=z+n^(-1/3)*(E[,1]-E[,2])/sqrt(2)+eps)
  rotate <- function(A,order) {
    for(j in order) {
      fit <- z*sum(z*A[,j])/zz; old <- A[,j]-fit
      fresh <- rnorm(n); fresh <- fresh-mean(fresh)
      fresh <- fresh-z*sum(z*fresh)/zz
      A[,j] <- fit+fresh*sqrt(sum(old*old)/sum(fresh*fresh))
    }
    A
  }
  hub <- rotate(X,1:2)
  candidates <- c(list(X),lapply(seq_len(M),function(m) rotate(hub,2:1)))
  rng_after_sampler <- .Random.seed
  base_cross <- as.vector(crossprod(z,X)); base_norm <- colSums(X*X)
  invariant_error <- max(vapply(candidates,function(A)
    max(abs(colMeans(A)),abs(as.vector(crossprod(z,A))-base_cross)/n,
        abs(colSums(A*A)-base_norm)/n),0.))
  stopifnot(is.finite(invariant_error),invariant_error<1e-10)
  HY <- H(Y); denom <- colSums(HY*HY)/n
  numerator <- matrix(NA_real_,M+1L,length(cases),dimnames=list(0:M,cases))
  C <- array(NA_real_,c(2L,M+1L,length(cases)),dimnames=list(NULL,0:M,cases))
  D <- vector('list',M+1L)
  for(j in seq_along(candidates)) {
    HX <- H(candidates[[j]]); G <- crossprod(HX)/n+diag(1/n,2)
    cv <- crossprod(HX,HY)/n
    numerator[j,] <- colSums(cv*solve(G,cv))
    C[,j,] <- cv; D[[j]] <- G
  }
  original <- sweep(numerator,2,pmax(denom,1e-30/n),'/')
  Dpool <- Reduce('+',D)/(M+1L)
  min_pool_eigenvalue <- min(eigen(Dpool,symmetric=TRUE,only.values=TRUE)$values)
  stopifnot(min_pool_eigenvalue>0)
  chol(Dpool)
  scores <- array(NA_real_,c(M+1L,length(cases),length(methods)),
                  dimnames=list(0:M,cases,methods))
  scores[,,'original'] <- original
  scores[,,'scalar_twosided'] <- sweep(original,2,colMeans(original),'-')^2
  permutation <- (M+1L):1L  # deterministic; consumes no random numbers
  Dperm <- Reduce('+',D[permutation])/(M+1L)
  check <- c(permutation_error=0,translation_error=0,m1_equidistance_error=0,
             permutation_rank_mismatches=0,translation_rank_mismatches=0,
             closest_observed_rank_gap=Inf)
  for(j in seq_along(cases)) {
    for(metric in c('centered_pool','centered_identity')) {
      G <- if(metric=='centered_pool') Dpool else diag(2)
      gp <- if(metric=='centered_pool') Dperm else diag(2)
      a <- center_scores(C[,,j],G)
      ap <- center_scores(C[,permutation,j],gp)
      at <- center_scores(sweep(C[,,j],1,c(.7,-.2),'+'),G)
      d1 <- Reduce('+',D[1:2])/2
      am1 <- center_scores(C[,1:2,j],if(metric=='centered_pool')d1 else diag(2))
      check['permutation_error'] <- max(check['permutation_error'],relative_error(ap,a[permutation]))
      check['translation_error'] <- max(check['translation_error'],relative_error(at,a))
      check['m1_equidistance_error'] <- max(check['m1_equidistance_error'],relative_error(am1[1],am1[2]))
      check['permutation_rank_mismatches'] <- check['permutation_rank_mismatches']+
        as.integer(!identical(rank_p(ap[order(permutation)]),rank_p(a)))
      check['translation_rank_mismatches'] <- check['translation_rank_mismatches']+
        as.integer(!identical(rank_p(at),rank_p(a)))
      check['closest_observed_rank_gap'] <- min(check['closest_observed_rank_gap'],min(abs(a[-1]-a[1])))
      scores[,j,metric] <- a
    }
  }
  stopifnot(all(is.finite(scores)),all(scores>=-1e-12),
            max(check[c('permutation_error','translation_error','m1_equidistance_error')])<1e-12,
            check['permutation_rank_mismatches']==0,check['translation_rank_mismatches']==0,
            identical(.Random.seed,rng_after_sampler))
  pvals <- apply(scores,c(2,3),rank_p)
  stopifnot(identical(as.numeric(pvals[,'original']),as.numeric(apply(numerator,2,rank_p))))
  list(n=n,replicate=b,seed=seed,scores=scores,c=C,D=D,Dpool=Dpool,
       numerator=numerator,denominator=denom,pvalue=pvals,reject=pvals<=cfg$alpha,
       invariant_error=invariant_error,min_pool_eigenvalue=min_pool_eigenvalue,
       checks=check,seconds=now()-time0,rng_after_sampler=rng_after_sampler)
}
for(k in seq_along(cfg$sample_sizes)) {
  n <- as.integer(cfg$sample_sizes[k]); block_start <- now()
  for(b in seq_len(reps)) {
    seed <- as.integer(cfg$seed_base+10000*k+b)
    ans <- tryCatch(one(n,b,seed),error=function(e)e)
    if(inherits(ans,'error')) {
      failures <- rbind(failures,data.frame(n=n,replicate=b,seed=seed,error=conditionMessage(ans)))
      write.csv(failures,'failures.csv',row.names=FALSE)
    } else records[[sprintf('n%d-r%04d',n,b)]] <- ans
  }
  saveRDS(records,sprintf('checkpoint-n%d.rds',n),version=3)
  cat(sprintf('Completed n=%d: %d planned replicates; cumulative failures=%d; elapsed=%.3fs\n',n,reps,nrow(failures),now()-block_start));flush.console()
}
saveRDS(records,'records.rds',version=3);write.csv(failures,'failures.csv',row.names=FALSE)
tests <- do.call(rbind,lapply(names(records),function(id) {
  r <- records[[id]]
  do.call(rbind,lapply(methods,function(me) data.frame(dataset_id=id,n=r$n,replicate=r$replicate,
    seed=r$seed,case=cases,method=me,pvalue=r$pvalue[,me],reject=r$reject[,me],
    T0=r$scores[1,,me],Tbar=colMeans(r$scores[-1,,me]),
    gap=r$scores[1,,me]-colMeans(r$scores[-1,,me]))))
}))
write.csv(tests,'tests.csv',row.names=FALSE)
scores <- do.call(rbind,lapply(names(records),function(id) {
  r <- records[[id]]
  do.call(rbind,lapply(methods,function(me)do.call(rbind,lapply(cases,function(ca)
    data.frame(dataset_id=id,n=r$n,replicate=r$replicate,seed=r$seed,case=ca,method=me,
               candidate=0:M,score=r$scores[,ca,me])))))
}))
write.csv(scores,'scores.csv',row.names=FALSE)
checks <- do.call(rbind,lapply(names(records),function(id) {
  r <- records[[id]]
  cbind(data.frame(dataset_id=id,n=r$n,replicate=r$replicate,seed=r$seed,
      seconds=r$seconds,invariant_error=r$invariant_error,min_pool_eigenvalue=r$min_pool_eigenvalue),
      as.data.frame(as.list(r$checks)))
}))
write.csv(checks,'checks.csv',row.names=FALSE)
# Exact replay of all three old p-values for the predeclared first 100 seeds.
old <- read.csv(file.path(Sys.getenv('REPLICATION_ROOT'),'reference/direction-prefix/tests.csv'))
old <- old[old$replicate<=reps,]
replay <- merge(old,tests[tests$method=='original' & tests$case!='orthogonal',],
  by=c('dataset_id','n','replicate','seed','case'),suffixes=c('_old','_new'))
replay$pvalue_equal <- replay$pvalue_old==replay$pvalue_new
replay$reject_equal <- replay$reject_old==replay$reject_new
write.csv(replay,'reference_replay.csv',row.names=FALSE)
stopifnot(nrow(failures)==0L,length(records)==length(cfg$sample_sizes)*reps,
  nrow(tests)==length(records)*length(cases)*length(methods),nrow(replay)==900L,
  all(replay$pvalue_equal),all(replay$reject_equal))
summary <- do.call(rbind,lapply(cfg$sample_sizes,function(n)
  do.call(rbind,lapply(cases,function(ca)do.call(rbind,lapply(methods,function(me) {
    z <- tests[tests$n==n & tests$case==ca & tests$method==me,]
    ci <- binom.test(sum(z$reject),reps)$conf.int
    data.frame(n=n,case=ca,method=me,replicates=reps,rejects=sum(z$reject),
      rate=mean(z$reject),ci_lower=ci[1],ci_upper=ci[2],mcse=sqrt(mean(z$reject)*(1-mean(z$reject))/reps),
      mean_gap=mean(z$gap),gap_mcse=sd(z$gap)/sqrt(reps))
  }))))))
write.csv(summary,'summary.csv',row.names=FALSE)
paired <- do.call(rbind,lapply(cfg$sample_sizes,function(n)
  do.call(rbind,lapply(cases,function(ca)do.call(rbind,lapply(methods[-1],function(me) {
    old <- tests[tests$n==n & tests$case==ca & tests$method=='original',]
    new <- tests[tests$n==n & tests$case==ca & tests$method==me,]
    stopifnot(identical(old$dataset_id,new$dataset_id))
    d <- as.numeric(new$reject)-as.numeric(old$reject)
    data.frame(n=n,case=ca,method=me,pairs=reps,new_minus_original=mean(d),
      paired_mcse=sd(d)/sqrt(reps),new_only=sum(new$reject & !old$reject),
      original_only=sum(old$reject & !new$reject))
  }))))))
write.csv(paired,'paired_differences.csv',row.names=FALSE)
validation <- list(datasets=length(records),failures=nrow(failures),tests=nrow(tests),
  reference_pvalues_compared=nrow(replay),reference_pvalue_mismatches=sum(!replay$pvalue_equal),
  reference_max_T0_difference=max(abs(replay$T0_new-replay$T0_old)),
  max_sampler_invariant_error=max(checks$invariant_error),
  minimum_pool_eigenvalue=min(checks$min_pool_eigenvalue),
  max_permutation_error=max(checks$permutation_error),
  max_translation_error=max(checks$translation_error),
  max_m1_equidistance_error=max(checks$m1_equidistance_error),
  permutation_rank_mismatches=sum(checks$permutation_rank_mismatches),
  translation_rank_mismatches=sum(checks$translation_rank_mismatches),
  closest_observed_centered_rank_gap=min(checks$closest_observed_rank_gap),
  score_rounding='none',m1_note='Only mathematical equidistance checked up to floating error; M=1 is excluded from the tested procedure.',
  seconds=now()-started)
jsonlite::write_json(validation,'VALIDATION.json',pretty=TRUE,auto_unbox=TRUE,digits=17)
writeLines(sprintf('Total wall seconds including checks and output: %.6f',now()-started),'elapsed.txt')
print(summary,row.names=FALSE)
cat('Completed exactly the predeclared settings; no replacements or additional experiments.\n')
