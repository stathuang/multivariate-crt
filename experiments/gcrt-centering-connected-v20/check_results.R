# Recompute all scores, ranks, intervals and paired uncertainty from saved components.
options(warn=2,digits=17)
stopifnot(basename(normalizePath('.'))=='gcrt-centering-connected-v20',file.exists('run/COMPLETE.json'))
cfg <- jsonlite::fromJSON('run/CONFIG.json')
tests <- read.csv('run/tests.csv'); rates <- read.csv('run/summary.csv')
paired <- read.csv('run/paired_differences.csv')
calc <- list(); maximum_score_error <- 0; pool_error <- 0; raw_rank_mismatch <- 0
files <- list.files('run/records',pattern='^batch-[0-9]+[.]rds$',full.names=TRUE)
stopifnot(length(files)==25)
for (file in files) {
  records <- readRDS(file)
  for (a in records) {
    K <- dim(a$c)[2]; q <- dim(a$c)[1]
    common_gram <- matrix(0,q,q)
    for (j in 1:K) common_gram <- common_gram + a$D[,,j]/K
    pool_error <- max(pool_error,max(abs(common_gram-a$Dpool)))
    for (cc in seq_len(dim(a$c)[3])) {
      values <- matrix(NA_real_,K,3,dimnames=list(NULL,c('original','centered_pool','centered_identity')))
      centered <- t(a$c[,,cc])-matrix(rowMeans(a$c[,,cc]),K,q,byrow=TRUE)
      values[,2] <- diag(centered%*%solve(common_gram)%*%t(centered))
      values[,3] <- rowSums(centered^2)
      for (j in 1:K) {
        z <- a$c[,j,cc]
        values[j,1] <- drop(t(z)%*%solve(a$D[,,j])%*%z)
      }
      maximum_score_error <- max(maximum_score_error,max(abs(values-a$scores[,cc,])))
      for (me in colnames(values)) {
        pv <- sum(values[,me]>=values[1,me])/K
        raw_rank_mismatch <- raw_rank_mismatch + as.integer(pv != a$pvalue[cc,me])
        calc[[length(calc)+1L]] <- data.frame(replicate=a$replicate,
          case=dimnames(a$c)[[3]][cc],method=me,pvalue=pv,reject=pv<=.05)
      }
    }
  }
}
calc <- do.call(rbind,calc)
key <- c('replicate','case','method')
stopifnot(nrow(calc)==60000,!anyDuplicated(calc[key]),!anyDuplicated(tests[key]),
          identical(sort(unique(calc$replicate)),1:2500),all(table(calc$case,calc$method)==2500))
z <- merge(calc,tests,by=key,all=TRUE,suffixes=c('_calc','_saved'))
stopifnot(!anyNA(z),nrow(z)==60000)
summary_error <- paired_error <- 0
for (i in 1:nrow(rates)) {
  a <- rates[i,]; y <- calc$reject[calc$case==a$case & calc$method==a$method]
  ci <- binom.test(sum(y),length(y))$conf.int
  summary_error <- max(summary_error,abs(sum(y)-a$rejects),abs(mean(y)-a$rate),
    abs(sqrt(mean(y)*(1-mean(y))/length(y))-a$mcse),abs(ci-c(a$ci_lower,a$ci_upper)))
}
for (i in 1:nrow(paired)) {
  a <- paired[i,]; s <- calc[calc$case==a$case,]
  x <- s[s$method==a$method_a,]; y <- s[s$method==a$method_b,]
  xy <- merge(x,y,by='replicate'); d <- as.numeric(xy$reject.x)-as.numeric(xy$reject.y)
  se <- sd(d)/sqrt(length(d)); delta <- mean(d)
  paired_error <- max(paired_error,abs(delta-a$difference),abs(se-a$paired_mcse),
    abs(delta-qnorm(.975)*se-a$ci_lower),abs(delta+qnorm(.975)*se-a$ci_upper),
    abs(sum(d==1)-a$a_only),abs(sum(d==-1)-a$b_only))
}
result <- list(checked_at=as.character(Sys.time()),replicates=2500,pvalues=nrow(calc),
  pvalue_mismatches=sum(z$pvalue_calc!=z$pvalue_saved),reject_mismatches=sum(z$reject_calc!=z$reject_saved),
  record_rank_mismatches=raw_rank_mismatch,maximum_score_error=maximum_score_error,
  maximum_pool_error=pool_error,summary_rows=nrow(rates),paired_rows=nrow(paired),
  maximum_summary_error=summary_error,maximum_paired_error=paired_error,
  scope='Separate deterministic recomputation from stored score vectors and candidate Gram matrices; no new sampling.')
stopifnot(result$pvalue_mismatches==0,result$reject_mismatches==0,raw_rank_mismatch==0,
          maximum_score_error<1e-10,summary_error<1e-12,paired_error<1e-12)
jsonlite::write_json(result,'run/RESULTS_CHECK.json',pretty=TRUE,auto_unbox=TRUE,digits=17)
cat(jsonlite::toJSON(result,pretty=TRUE,auto_unbox=TRUE,digits=17),'\n')
