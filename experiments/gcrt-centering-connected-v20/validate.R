# Separate algebraic implementation checks; no outcomes are used to select settings.
options(warn=2,digits=17)
stopifnot(basename(normalizePath('.'))=='gcrt-centering-connected-v20')
cfg <- jsonlite::fromJSON('run/CONFIG.json')
RNGkind(cfg$rng[1],cfg$rng[2],cfg$rng[3])
source('core.R')
time0 <- proc.time()[['elapsed']]
x <- one_replicate(1,keep_data=TRUE)
dat <- x$data; X <- dat$X; S <- dat$S; Y <- dat$Y
# Reference local update via explicit projection and a common deterministic vector.
fresh <- sin(seq_len(cfg$n))
update_error <- 0
for (j in 1:8) {
  A <- cbind(1,S[,1],X[,neighbours[[j]],drop=FALSE])
  P <- A%*%solve(crossprod(A),t(A)); R <- diag(cfg$n)-P
  old <- R%*%X[,j]; z <- R%*%fresh
  ref <- P%*%X[,j]+z*sqrt(sum(old^2)/sum(z^2))
  update_error <- max(update_error,abs(update_target(X,S[,1],j,fresh)[,j]-ref))
}
# Independent primal ridge solves, rather than the production dual residual operator.
Sc <- center(S); Yc <- center(Y)
reg <- crossprod(Sc)+diag(.2*cfg$n,cfg$s)
RY <- Yc - Sc%*%solve(reg,crossprod(Sc,Yc))
c_ref <- x$c; D_ref <- x$D; original_ref <- x$scores[,,'original']
for (j in 1:(cfg$M+1)) {
  A <- center(dat$candidates[[j]])
  RX <- A - Sc%*%solve(reg,crossprod(Sc,A))
  c_ref[,j,] <- crossprod(RX,RY)/cfg$n
  D_ref[,,j] <- crossprod(RX)/cfg$n+diag(1/cfg$n,cfg$t)
  coef <- solve(crossprod(RX)+diag(cfg$t),crossprod(RX,RY))
  original_ref[j,] <- (colSums(RY^2)-colSums((RY-RX%*%coef)^2)-colSums(coef^2))/cfg$n
}
pool_ref <- Reduce('+',lapply(1:100,function(j) D_ref[,,j]))/100
pool_scores <- unit_scores <- original_ref
for (ca in 1:8) {
  delta <- t(c_ref[,,ca])-matrix(rowMeans(c_ref[,,ca]),100,8,byrow=TRUE)
  pool_scores[,ca] <- diag(delta%*%solve(pool_ref)%*%t(delta))
  unit_scores[,ca] <- rowSums(delta^2)
}
permutation <- 100:1
permuted <- score_all(x$c[,permutation,],x$D[,,permutation])$scores
translated <- score_all(sweep(x$c,1,seq_len(8)/10,'+'),x$D)$scores
repeat_x <- one_replicate(1)
errors <- c(update=max(update_error),score_vector=max(abs(c_ref-x$c)),
  candidate_Gram=max(abs(D_ref-x$D)),original=max(abs(original_ref-x$scores[,,'original'])),
  pooled=max(abs(pool_scores-x$scores[,,'centered_pool'])),
  unit=max(abs(unit_scores-x$scores[,,'centered_identity'])),
  permutation=max(abs(permuted-x$scores[permutation,,])),
  translation=max(abs(translated[,,2:3]-x$scores[,,2:3])))
rank_mismatch <- sum(apply(original_ref,2,rank_p)!=x$pvalue[,'original']) +
  sum(apply(pool_scores,2,rank_p)!=x$pvalue[,'centered_pool']) +
  sum(apply(unit_scores,2,rank_p)!=x$pvalue[,'centered_identity'])
stopifnot(max(errors)<1e-10,rank_mismatch==0,
  identical(repeat_x$scores,x$scores),identical(repeat_x$rng_after_sampler,x$rng_after_sampler))
elapsed <- proc.time()[['elapsed']]-time0
jsonlite::write_json(list(passed=TRUE,checked_at=as.character(Sys.time()),
  errors=as.list(errors),rank_mismatch=rank_mismatch,
  sufficient_statistics_error=x$invariant_error,seed_replay_exact=TRUE,
  one_production_replicate_seconds=repeat_x$seconds,
  estimated_2500_replicates_seconds=2500*repeat_x$seconds,
  validation_seconds=elapsed,replicate_used=1,
  scope='Checks reuse the first prespecified seed; rejection rates were not inspected or used for design changes.'),
  'run/IMPLEMENTATION_CHECK.json',pretty=TRUE,auto_unbox=TRUE,digits=17)
cat('Implementation checks passed. Seconds per replicate:',repeat_x$seconds,'\n')
