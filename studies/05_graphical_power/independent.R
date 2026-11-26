center <- function(A) sweep(A,2,colMeans(A),'-')
now <- function() proc.time()[['elapsed']]
rank_p <- function(x) (1+sum(x[-1L]>=x[1L]))/length(x)
center_scores <- function(C,G) {
  delta <- sweep(C,1,rowMeans(C),'-')
  colSums(delta*solve(G,delta))
}
relative_error <- function(x,y) max(abs(x-y))/max(1,max(abs(x)),max(abs(y)))
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
