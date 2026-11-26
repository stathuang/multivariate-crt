# Enabled only by explicit Gaussian cache policy in newly versioned campaigns.
# Gaussian hdi nodewise fits are design-only; non-Gaussian transforms may use y.
# See docs/CACHE_CANDIDATE.md before considering integration.
hdi_gaussian_cache <- function(X) {
  .baseline_require('hdi')
  if (!is.matrix(X) || !is.numeric(X) || any(!is.finite(X))) stop('Invalid X')
  if (any(apply(X,2,stats::sd) <= 0)) stop('Constant design column')
  # This matches prepare.data(..., family='gaussian', standardize=TRUE).
  # The internal builder is version-recorded; public lasso.proj consumes Z.
  standardized <- scale(scale(X,center=TRUE,scale=TRUE),center=TRUE,scale=FALSE)
  calculate <- getFromNamespace('calculate.Z','hdi')
  fitted <- calculate(x=standardized,parallel=FALSE,ncores=1L,verbose=FALSE,Z=NULL,do.ZnZ=FALSE)
  structure(list(X=X,Z=fitted$Z,hdi_version=as.character(utils::packageVersion('hdi')),
                 standardize=TRUE,family='gaussian',policy='one-nodewise-fit-per-design-v1'),
            class='hdi_gaussian_cache')
}

baseline_test_cached_candidate <- function(X,y,target,method,response_family='gaussian',cache,
                                          suppress_unused_group_testing=TRUE) {
  if (!inherits(cache,'hdi_gaussian_cache') || !identical(X,cache$X))
    stop('Cache design mismatch: use precisely the original matrix and column order')
  if (!identical(cache$hdi_version,as.character(utils::packageVersion('hdi'))))
    stop('hdi version differs from cache builder')
  if (method!='De-Lasso' || response_family!='gaussian')
    stop('Candidate cache supports Gaussian De-Lasso only')
  if (!is.numeric(y) || length(y)!=nrow(X) || any(!is.finite(y))) stop('Invalid response')
  if (!is.numeric(target) || !length(target) || anyNA(target) || anyDuplicated(target) ||
      any(target!=as.integer(target)) || any(target<1 | target>ncol(X))) stop('Invalid target')
  fit <- hdi::lasso.proj(X,y,family='gaussian',standardize=TRUE,
                        multiplecorr.method='holm',Z=cache$Z,
                        suppress.grouptesting=suppress_unused_group_testing)
  .baseline_bonferroni(fit$pval[target])
}
