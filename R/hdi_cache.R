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

