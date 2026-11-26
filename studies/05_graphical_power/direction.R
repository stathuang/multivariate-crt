rank_p <- function(x) (1+sum(x[-1L]>=x[1L]))/length(x)
center <- function(A) sweep(A,2,colMeans(A),"-")
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
  list(n=n,replicate=b,seed=seed,original=original[,1:3,drop=FALSE],
       numerator=num[,1:3,drop=FALSE],denominator=denom[1:3])
}
