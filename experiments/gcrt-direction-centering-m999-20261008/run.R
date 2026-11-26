options(warn=2, digits=17)
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
stopifnot(basename(normalizePath(getwd()))=="gcrt-direction-centering-m999-20261008",
          Sys.readlink(getwd())=="", !file.exists("RUN_STARTED.json"),
          requireNamespace("jsonlite",quietly=TRUE))
cfg <- jsonlite::fromJSON("CONFIG.json")
package_root <- Sys.getenv("REPLICATION_ROOT")
stopifnot(nzchar(package_root))
cfg$direction_reference <- file.path(package_root,"reference/direction-prefix")
cfg$centering_reference <- file.path(package_root,"reference/centering-prefix")
stopifnot(cfg$M==999L,cfg$L==1L,cfg$direction_repetitions_per_n==500L,
          cfg$centering_repetitions_per_n==100L,cfg$background_penalty_coefficient==1,
          cfg$target_penalty==1,identical(as.numeric(cfg$sample_sizes),c(40,120,360)))
M <- cfg$M; N <- M+1L; methods <- cfg$methods; cases <- cfg$cases
RNGkind(cfg$rng[1],cfg$rng[2],cfg$rng[3])
old_direction <- readRDS(file.path(cfg$direction_reference,"records.rds"))
old_center <- readRDS(file.path(cfg$centering_reference,"records.rds"))
started <- proc.time()[["elapsed"]]
jsonlite::write_json(list(started=as.character(Sys.time()),pid=Sys.getpid()),
                    "RUN_STARTED.json",pretty=TRUE,auto_unbox=TRUE)
capture.output(list(system=Sys.info(),RNG=RNGkind(),session=sessionInfo()),
               file="environment.txt")
dir.create("records"); dir.create("binary")
rank_p <- function(x) (1+sum(x[-1L]>=x[1L]))/length(x)
center <- function(A) sweep(A,2,colMeans(A),"-")
make_scores <- function(C,D,original) {
  nv <- nrow(original); nc <- ncol(original)
  ans <- array(NA_real_,c(nv,nc,4L),dimnames=list(0:(nv-1L),colnames(original),methods))
  ans[,,"original"] <- original
  ans[,,"scalar_twosided"] <- sweep(original,2,colMeans(original),"-")^2
  pool <- Reduce("+",D)/nv
  for(j in seq_len(nc)) {
    delta <- sweep(C[,,j],1,rowMeans(C[,,j]),"-")
    ans[,j,"centered_pool"] <- colSums(delta*solve(pool,delta))
    ans[,j,"centered_identity"] <- colSums(delta^2)
  }
  ans
}
one <- function(n,b,seed,id) {
  t0 <- proc.time()[["elapsed"]]; set.seed(seed)
  raw <- matrix(rnorm(4L*n),n,4L); centered <- center(raw)
  z <- centered[,1]; E <- centered[,2:3,drop=FALSE]; eps <- centered[,4]
  X <- E+z; zz <- sum(z*z); h <- n/(zz+n)
  H <- function(A) A-(1-h)*tcrossprod(z,as.vector(crossprod(z,A))/zz)
  Y <- sapply(c(null=0,positive=.25,negative=-.25),
              function(beta) z+beta*rowSums(E)+eps)
  if(b<=cfg$centering_repetitions_per_n)
    Y <- cbind(Y,orthogonal=z+n^(-1/3)*(E[,1]-E[,2])/sqrt(2)+eps)
  rotate <- function(A,order) {
    for(j in order) {
      fit <- z*sum(z*A[,j])/zz; old <- A[,j]-fit
      fresh <- rnorm(n); fresh <- fresh-mean(fresh)
      fresh <- fresh-z*sum(z*fresh)/zz
      A[,j] <- fit+fresh*sqrt(sum(old*old)/sum(fresh*fresh))
    }
    A
  }
  hub <- rotate(X,1:2); candidates <- vector("list",N); candidates[[1L]] <- X
  for(m in seq_len(M)) {
    candidates[[m+1L]] <- rotate(hub,2:1)
    if(m==19L) prefix_rng <- .Random.seed
  }
  final_rng <- .Random.seed
  base_cross <- as.vector(crossprod(z,X)); base_norm <- colSums(X*X)
  invariant <- max(vapply(candidates,function(A)
    max(abs(colMeans(A)),abs(as.vector(crossprod(z,A))-base_cross)/n,
        abs(colSums(A*A)-base_norm)/n),0.))
  stopifnot(is.finite(invariant),invariant<1e-10)
  HY <- H(Y); denom <- colSums(HY*HY)/n; nc <- ncol(Y)
  num <- matrix(NA_real_,N,nc,dimnames=list(0:M,colnames(Y)))
  C <- array(NA_real_,c(2L,N,nc)); D <- vector("list",N); direct_error <- 0
  for(v in seq_len(N)) {
    HX <- H(candidates[[v]]); G <- crossprod(HX)/n+diag(1/n,2)
    cv <- crossprod(HX,HY)/n; C[,v,] <- cv; D[[v]] <- G
    num[v,] <- colSums(cv*solve(G,cv))
    if(b==1L) {
      betahat <- solve(crossprod(HX)+diag(2),crossprod(HX,HY))
      direct <- (colSums(HY^2)-colSums((HY-HX%*%betahat)^2)-
                 colSums(betahat^2))/pmax(colSums(HY^2),1e-30)
      direct_error <- max(direct_error,max(abs(direct-num[v,]/pmax(denom,1e-30/n))))
    }
  }
  original <- sweep(num,2,pmax(denom,1e-30/n),"/")
  stopifnot(all(is.finite(original)),all(original>=-1e-12),direct_error<1e-12)
  p <- apply(original,2,rank_p)
  stopifnot(identical(as.numeric(p),as.numeric(apply(num,2,rank_p))))
  ref <- old_direction[[id]]; stopifnot(ref$seed==seed)
  prefix_error <- max(abs(original[1:20,1:3]-ref$scores))
  prefix_p <- apply(original[1:20,1:3],2,rank_p)
  stopifnot(prefix_error<1e-12,identical(as.numeric(prefix_p),as.numeric(ref$pvalue)))
  scores <- prefix_scores <- NULL; prefix_center_error <- 0
  if(b<=cfg$centering_repetitions_per_n) {
    scores <- make_scores(C,D,original)
    prefix_scores <- make_scores(C[,1:20,,drop=FALSE],D[1:20],original[1:20,,drop=FALSE])
    cref <- old_center[[id]]; prefix_center_error <- max(abs(prefix_scores-cref$scores))
    stopifnot(prefix_center_error<1e-12,identical(prefix_rng,cref$rng_after_sampler),
      identical(as.numeric(apply(prefix_scores,c(2,3),rank_p)),as.numeric(cref$pvalue)))
    stopifnot(all(is.finite(scores)),all(scores>=-1e-12))
  }
  stopifnot(identical(.Random.seed,final_rng))
  list(n=n,replicate=b,seed=seed,original=original[,1:3,drop=FALSE],
       numerator=num[,1:3,drop=FALSE],denominator=denom[1:3],scores=scores,
       c=if(!is.null(scores))C else NULL,D=if(!is.null(scores))D else NULL,
       denominator_all=if(!is.null(scores))denom else NULL,
       data=list(z=z,E=E,epsilon=eps,hub=hub),rng_after_sampler=final_rng,
       check=list(invariant_error=invariant,direct_error=direct_error,
                  prefix_original_error=prefix_error,prefix_center_error=prefix_center_error),
       seconds=proc.time()[["elapsed"]]-t0)
}
direction_rows <- center_rows <- check_rows <- list()
failures <- data.frame(n=integer(),replicate=integer(),seed=integer(),error=character())
for(k in seq_along(cfg$sample_sizes)) {
  n <- cfg$sample_sizes[k]; rec <- list()
  connections <- lapply(c("direction","centering","geometry"),function(kind)
    file(sprintf("binary/%s-n%d.f64",kind,n),"wb"))
  for(b in seq_len(cfg$direction_repetitions_per_n)) {
    seed <- as.integer(cfg$seed_base+10000*k+b); id <- sprintf("n%d-r%04d",n,b)
    r <- tryCatch(one(n,b,seed,id),error=function(e)e)
    if(inherits(r,"error")) {
      failures <- rbind(failures,data.frame(n=n,replicate=b,seed=seed,error=conditionMessage(r)))
      write.csv(failures,"failures.csv",row.names=FALSE)
      cat("FAILED",id,conditionMessage(r),"\n"); flush.console(); next
    }
    rec[[id]] <- r
    writeBin(as.double(r$original),connections[[1]],size=8,endian="little")
    op <- apply(r$original,2,rank_p)
    direction_rows[[id]] <- data.frame(dataset_id=id,n=n,replicate=b,seed=seed,
      case=cases[1:3],pvalue=op,reject=op<=cfg$alpha,T0=r$original[1,],
      Tbar=colMeans(r$original[-1,,drop=FALSE]),
      gap=r$original[1,]-colMeans(r$original[-1,,drop=FALSE]))
    if(!is.null(r$scores)) {
      writeBin(as.double(r$scores),connections[[2]],size=8,endian="little")
      writeBin(c(as.double(r$c),as.double(unlist(r$D)),r$denominator_all),
               connections[[3]],size=8,endian="little")
      center_rows[[id]] <- do.call(rbind,lapply(methods,function(me) {
        q <- apply(r$scores[,,me],2,rank_p)
        data.frame(dataset_id=id,n=n,replicate=b,seed=seed,case=cases,method=me,
          pvalue=q,reject=q<=cfg$alpha,T0=r$scores[1,,me],
          Tbar=colMeans(r$scores[-1,,me]),gap=r$scores[1,,me]-colMeans(r$scores[-1,,me]))
      }))
    }
    check_rows[[id]] <- cbind(data.frame(dataset_id=id,n=n,replicate=b,seed=seed,
      seconds=r$seconds),as.data.frame(r$check))
    if(b%%100L==0L) {cat(sprintf("n=%d completed %d/500, failures=%d, elapsed=%.1fs\n",
      n,b,nrow(failures),proc.time()[["elapsed"]]-started));flush.console()}
  }
  lapply(connections,close)
  saveRDS(rec,sprintf("records/n%d.rds",n),compress=FALSE,version=3)
}
direction <- do.call(rbind,direction_rows); centering <- do.call(rbind,center_rows)
checks <- do.call(rbind,check_rows)
write.csv(direction,"direction-tests.csv",row.names=FALSE)
write.csv(centering,"centering-tests.csv",row.names=FALSE)
write.csv(checks,"checks.csv",row.names=FALSE);write.csv(failures,"failures.csv",row.names=FALSE)
stopifnot(nrow(failures)==0L,nrow(direction)==4500L,nrow(centering)==4800L,nrow(checks)==1500L)
summary_rows <- function(dat,groups) {
  keys <- unique(dat[groups])
  do.call(rbind,lapply(seq_len(nrow(keys)),function(i) {
    use <- rep(TRUE,nrow(dat))
    for(g in groups) use <- use & dat[[g]]==keys[[g]][i]
    z <- dat[use,];ci <- binom.test(sum(z$reject),nrow(z))$conf.int
    cbind(keys[i,,drop=FALSE],data.frame(replicates=nrow(z),rejects=sum(z$reject),
      rate=mean(z$reject),ci_lower=ci[1],ci_upper=ci[2],
      mean_T0=mean(z$T0),mean_Tbar=mean(z$Tbar),mean_gap=mean(z$gap),
      gap_mcse=sd(z$gap)/sqrt(nrow(z))))
  }))
}
ds <- summary_rows(direction,c("n","case"))
cs <- summary_rows(centering,c("n","case","method"))
write.csv(ds,"direction-summary.csv",row.names=FALSE)
write.csv(cs,"centering-summary.csv",row.names=FALSE)
jsonlite::write_json(list(datasets=1500,direction_test_rows=4500,centering_test_rows=4800,
  shared_test_rows=900,failures=0,max_invariant_error=max(checks$invariant_error),
  max_direct_error=max(checks$direct_error),prefix_original_scores_compared=90000,
  prefix_center_scores_compared=96000,prefix_pvalue_mismatches=0,
  max_prefix_original_error=max(checks$prefix_original_error),
  max_prefix_center_error=max(checks$prefix_center_error),
  computation_seconds=proc.time()[["elapsed"]]-started),
  "RUN_COMPLETE.json",pretty=TRUE,auto_unbox=TRUE,digits=17)
print(ds,row.names=FALSE)
cat("Completed fixed larger-copy experiment. Original results and manuscript unchanged.\n")
