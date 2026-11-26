# Connected-target data generator, sampler, and quadratic statistics; caller supplies cfg.
methods <- c('original', 'centered_pool', 'centered_identity')
cells <- expand.grid(r = cfg$signal, b = cfg$background_sign)
cells$case <- sprintf('b%s_r%02d', ifelse(cells$b < 0, 'neg', 'pos'), round(100*cells$r))
Th <- diag(cfg$t)
Th[abs(row(Th) - col(Th)) == 1L] <- -.2
V <- solve(Th); cholV <- chol(V)
v <- sqrt(2/9) * sin((1:8)*pi/9)
lambda_star <- 1 - .4*cos(pi/9)
neighbours <- lapply(1:8, function(j) intersect(c(j-1L,j+1L), 1:8))
center <- function(x) sweep(x, 2L, colMeans(x), '-')
rank_p <- function(s) (1 + sum(s[-1L] >= s[1L]))/length(s)
update_target <- function(X, S1, j, fresh = rnorm(nrow(X))) {
  A <- cbind(1, S1, X[, neighbours[[j]], drop=FALSE])
  qrA <- qr(A)
  stopifnot(qrA$rank == ncol(A))
  old <- qr.resid(qrA, X[,j]); fitted <- X[,j] - old
  new <- qr.resid(qrA, fresh)
  denom <- sum(new^2); oldnorm <- sum(old^2)
  stopifnot(is.finite(denom), denom > 0, oldnorm > 0)
  X[,j] <- fitted + new*sqrt(oldnorm/denom)
  X
}
scan_targets <- function(X, S1, order) {
  for (j in order) X <- update_target(X,S1,j)
  X
}
sufficient <- function(X,S1) {
  G <- crossprod(X)
  c(colMeans(X), diag(G)/nrow(X), G[cbind(1:7,2:8)]/nrow(X),
    as.numeric(crossprod(S1,X))/nrow(X))
}
score_all <- function(cv,D) {
  K <- dim(cv)[2L]; q <- dim(cv)[1L]; nc <- dim(cv)[3L]
  pool <- apply(D,c(1,2),mean)
  score <- array(NA_real_,c(K,nc,3L),dimnames=list(0:(K-1),cells$case,methods))
  for (j in seq_len(K)) {
    a <- cv[,j,]
    score[j,,'original'] <- colSums(a*solve(D[,,j],a))
  }
  for (ca in seq_len(nc)) {
    delta <- sweep(cv[,,ca],1,rowMeans(cv[,,ca]),'-')
    score[,ca,'centered_pool'] <- colSums(delta*solve(pool,delta))
    score[,ca,'centered_identity'] <- colSums(delta^2)
  }
  list(scores=score,pool=pool)
}
one_replicate <- function(k, keep_data=FALSE) {
  time0 <- proc.time()[['elapsed']]
  n <- cfg$n; t <- cfg$t; M <- cfg$M
  data_seed <- as.integer(cfg$data_seed_base + k)
  sampler_seed <- as.integer(cfg$sampler_seed_base + k)
  set.seed(data_seed)
  S <- matrix(rnorm(n*cfg$s),n,cfg$s)
  E <- matrix(rnorm(n*t),n,t) %*% cholV
  eps <- rnorm(n)
  X <- E + .5*tcrossprod(S[,1],v)
  Y <- outer(S[,1],cells$b) + outer(as.numeric(E%*%v)*sqrt(lambda_star),cells$r) + eps
  colnames(Y) <- cells$case
  Sc <- center(S); Yc <- center(Y)
  penalty <- cfg$background_penalty_coefficient*n
  H <- penalty*chol2inv(chol(tcrossprod(Sc) + diag(penalty,n)))
  HY <- H%*%Yc
  set.seed(sampler_seed)
  hub <- scan_targets(X,S[,1],seq_len(t))
  candidates <- c(list(X),lapply(seq_len(M),function(m) scan_targets(hub,S[,1],t:1)))
  rng_after <- .Random.seed
  baseline <- sufficient(X,S[,1])
  invariant_error <- max(vapply(candidates,function(A)
    max(abs(sufficient(A,S[,1])-baseline)),0.))
  cvec <- array(NA_real_,c(t,M+1L,nrow(cells)),dimnames=list(NULL,0:M,cells$case))
  D <- array(NA_real_,c(t,t,M+1L))
  for (j in seq_along(candidates)) {
    HX <- H%*%center(candidates[[j]])
    cvec[,j,] <- crossprod(HX,HY)/n
    D[,,j] <- crossprod(HX)/n + diag(cfg$target_penalty/n,t)
  }
  ans <- score_all(cvec,D)
  pv <- apply(ans$scores,c(2,3),rank_p)
  min_eig <- min(eigen(ans$pool,symmetric=TRUE,only.values=TRUE)$values)
  stopifnot(invariant_error < 1e-9, min_eig > 0, all(is.finite(ans$scores)),
            all(ans$scores >= -1e-12), identical(rng_after,.Random.seed))
  ret <- list(replicate=k,data_seed=data_seed,sampler_seed=sampler_seed,
    c=cvec,D=D,Dpool=ans$pool,scores=ans$scores,pvalue=pv,reject=pv<=cfg$alpha,
    invariant_error=invariant_error,min_pool_eigenvalue=min_eig,
    seconds=proc.time()[['elapsed']]-time0,rng_after_sampler=rng_after)
  if (keep_data) ret$data <- list(S=S,E=E,eps=eps,X=X,Y=Y,H=H,candidates=candidates)
  ret
}
test_rows <- function(x) {
  z <- expand.grid(case=cells$case,method=methods,stringsAsFactors=FALSE)
  z$replicate <- x$replicate; z$data_seed <- x$data_seed; z$sampler_seed <- x$sampler_seed
  z$b <- cells$b[match(z$case,cells$case)]; z$r <- cells$r[match(z$case,cells$case)]
  z$pvalue <- as.vector(x$pvalue); z$reject <- as.vector(x$reject)
  z
}
