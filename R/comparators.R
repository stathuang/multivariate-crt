# Keep the accepted Gaussian De-Lasso cache and its two independent RNG streams.
# Restoring the caller's state prevents these fits from changing later methods.
with_seed <- function(seed, expr) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) previous <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) assign(".Random.seed", previous, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)
  force(expr)
}

de_lasso_seed <- function(config, replicate, stage, theta_index = 0L) {
  offset <- switch(stage, nodewise = 1357911,
                   response = 2468022 + 104729 * theta_index,
                   stop("Unknown De-Lasso RNG stage"))
  as.integer((as.double(config$seed_base) + replicate + offset) %% 2147483646 + 1)
}

# Build once per X, before sampling copies, as in the accepted reconstruction.
prepare_de_lasso <- function(config, replicate, X) {
  policy <- config$de_lasso_policy
  if (is.null(policy) || policy == "uncached-ordered-v1") return(list(enabled = FALSE))
  if (policy != "one-nodewise-fit-per-design-v1") stop("Unknown De-Lasso cache policy")
  if (config$response_family != "gaussian" || !"De-Lasso" %in% config$baselines)
    stop("De-Lasso cache policy requires Gaussian De-Lasso in the method list")
  if (!identical(config$de_lasso_seed_policy, "separate-nodewise-response-seeds-v1"))
    stop("The Gaussian cache requires its explicit separate-stream policy")
  suppress_groups <- config$de_lasso_suppress_unused_group_testing
  if (!is.logical(suppress_groups) || length(suppress_groups) != 1L || is.na(suppress_groups))
    stop("Specify the unused-group-testing policy explicitly")
  seed <- de_lasso_seed(config, replicate, "nodewise")
  result <- capture_step(with_seed(seed, hdi_gaussian_cache(X)))
  list(enabled = TRUE, nodewise_seed = seed, result = result)
}

de_lasso_test <- function(X, Y, target, cache = NULL, suppress_groups = TRUE) {
  check_test_inputs(X, Y, target)
  if (is.null(cache)) {
    fit <- hdi::lasso.proj(X, Y, family = "gaussian", multiplecorr.method = "holm")
  } else {
    if (!inherits(cache, "hdi_gaussian_cache") || !identical(X, cache$X))
      stop("Cache design mismatch: use precisely the original matrix and column order")
    if (!identical(cache$hdi_version, as.character(utils::packageVersion("hdi"))))
      stop("hdi version differs from cache builder")
    fit <- hdi::lasso.proj(X, Y, family = "gaussian", standardize = TRUE,
                          multiplecorr.method = "holm", Z = cache$Z,
                          suppress.grouptesting = suppress_groups)
  }
  individual <- fit$pval[target]
  coordinate_z <- abs(fit$bhat[target] / fit$se[target])
  # The group p-value still uses hdi's original p-values, not recomputed tails.
  list(pvalue = .baseline_bonferroni(individual), statistic = max(coordinate_z),
       diagnostics = list(coordinate_z = coordinate_z, coordinate_pvalues = individual,
                          combination = "Bonferroni"))
}

cgm_linear_test <- function(X, Y, target) {
  check_test_inputs(X, Y, target)
  fit <- cgm_lf_diagnostic(X, Y, diag(ncol(X))[, target, drop = FALSE],
                          model = "linear", intercept = TRUE, intercept.loading = FALSE)
  coordinate_z <- abs(fit$est.debias.vec / fit$se.vec)
  individual <- 2 * stats::pnorm(coordinate_z, lower.tail = FALSE)
  list(pvalue = .baseline_bonferroni(individual), statistic = max(coordinate_z),
       diagnostics = list(coordinate_z = coordinate_z, coordinate_pvalues = individual,
                          combination = "Bonferroni"))
}
