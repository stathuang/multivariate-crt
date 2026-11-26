# D.10 covariance recipe from Liu et al., with an estimated intercept.
# Reference: reference/dcrt-author-functions.R, linshrink_cov (lines 1221-1256).
lw_covariance <- function(X) {
  n <- nrow(X)
  p <- ncol(X)
  if (n < 3L || p < 2L || any(!is.finite(X))) stop("Invalid covariate matrix")
  mu <- colMeans(X)
  centered <- sweep(X, 2L, mu, "-")
  effective_n <- n - 1L
  S <- crossprod(centered) / effective_n
  target_scale <- sum(diag(S)) / p
  distance <- sum((S - diag(target_scale, p))^2) / p
  # Exact algebraic expansion of the author's rowwise outer-product sum.
  numerator <- sum(rowSums(centered^2)^2) -
    2 * sum(crossprod(centered) * S) + n * sum(S^2)
  variance <- max(0, numerator / (p * effective_n^2))
  if (!is.finite(target_scale) || target_scale <= 0) stop("Zero covariate variance")
  shrinkage <- if (distance == 0) 1 else min(distance, variance) / distance
  unadjusted <- shrinkage * diag(target_scale, p) + (1 - shrinkage) * S
  precision <- solve(unadjusted)
  # Column j equals X_j - X_{-j} beta_j under the unadjusted covariance.
  residuals <- sweep(centered %*% precision, 2L, diag(precision), "/")
  residual_variance <- colMeans(residuals^2)
  model_variance <- 1 / diag(precision)
  diagonal_scale <- sqrt(residual_variance / model_variance)
  adjusted <- unadjusted * outer(diagonal_scale, diagonal_scale)
  if (any(!is.finite(adjusted)) || any(diagonal_scale <= 0))
    stop("Undefined D.10 variance adjustment")
  chol(adjusted)
  list(mean = mu, covariance = adjusted, unadjusted = unadjusted,
       shrinkage = shrinkage, scale = diagonal_scale,
       residual_variance = residual_variance, model_variance = model_variance)
}

# Apply the literature covariance estimator to a target BLOCK given backgrounds.
# D.10 uses scalar targets; this Gaussian block conditioning is our adaptation.
lw_conditional <- function(X, s) {
  fit <- lw_covariance(X)
  background <- seq_len(s)
  target <- seq.int(s + 1L, ncol(X))
  Sigma <- fit$covariance
  beta <- solve(Sigma[background, background, drop = FALSE],
                Sigma[background, target, drop = FALSE])
  mean <- sweep(sweep(X[, background, drop = FALSE], 2L,
                     fit$mean[background], "-") %*% beta,
                2L, fit$mean[target], "+")
  covariance <- Sigma[target, target, drop = FALSE] -
    Sigma[target, background, drop = FALSE] %*% beta
  covariance <- (covariance + t(covariance)) / 2
  list(mean = mean, covariance = covariance, root = chol(covariance), fit = fit)
}

lw_copy_stats <- function(dat, s, M, lambda, power = FALSE) {
  conditional <- lw_conditional(dat$X, s)
  n <- nrow(dat$X)
  t_dim <- ncol(dat$G)
  context <- statistic_context(dat$PR, dat$rY)
  stats <- if (power) numeric(M) else
    matrix(NA_real_, M, 3L, dimnames = list(NULL, c("raw score", "poly2 score", "ridge gain")))
  for (m in seq_len(M)) {
    Xtilde <- conditional$mean + matrix(rnorm(n * t_dim), n, t_dim) %*% conditional$root
    G <- dat$PR %*% Xtilde
    if (power) stats[m] <- ridge_stat(G, context$rd, lambda) else
      stats[m, ] <- stat_values(G, context, lambda)
  }
  list(statistics = stats, fit = conditional)
}

append_lw <- function(original, cfg, M, lambda, power, keep_example) {
  # Added AFTER the four original methods: their data, p-values and RNG path stay intact.
  rng_after_original <- .Random.seed
  new <- capture_step(lw_copy_stats(original$diagnostics$data, cfg$s, M, lambda, power))
  name <- "LW plug-in CRT (D.10)"
  if (power) {
    p <- if (new$status == "ok") mc_pvalue(original$statistic, new$value$statistics) else NA_real_
    original$pvalue <- c(original$pvalue, setNames(p, name))
  } else {
    p <- if (new$status == "ok") vapply(seq_along(original$statistic), function(j)
      mc_pvalue(original$statistic[j], new$value$statistics[, j]), 0) else rep(NA_real_, 3L)
    original$pvalue <- cbind(original$pvalue, p)
    colnames(original$pvalue)[ncol(original$pvalue)] <- name
  }
  diagnostic <- list(lw_status = new$status, lw_error = new$error, lw_warnings = new$warnings,
                     rng_after_original = rng_after_original)
  if (new$status == "ok") {
    fit <- new$value$fit$fit
    diagnostic$shrinkage <- fit$shrinkage
    diagnostic$scale <- fit$scale
    diagnostic$conditional_variance <- diag(new$value$fit$covariance)
  }
  if (keep_example) {
    original$diagnostics$randomized_statistics[[name]] <- if (new$status == "ok") new$value$statistics else NULL
    original$diagnostics$lw_fit <- if (new$status == "ok") new$value$fit else NULL
    original$diagnostics <- c(original$diagnostics, diagnostic)
  } else original$diagnostics <- diagnostic
  original
}

size_replication_lw <- function(b, cfg, coefs, setting_index, seed_base, M,
                                lambda_ridge, keep_example = FALSE) {
  original <- size_replication(b, cfg, coefs, setting_index, seed_base, M,
                               lambda_ridge, keep_example = TRUE)
  append_lw(original, cfg, M, lambda_ridge, FALSE, keep_example)
}

power_replication_lw <- function(b, cfg, coefs, theta, setting_index, theta_index,
                                 seed_base, M, lambda_ridge, keep_example = FALSE) {
  original <- power_replication(b, cfg, coefs, theta, setting_index, theta_index,
                                seed_base, M, lambda_ridge, keep_example = TRUE)
  append_lw(original, cfg, M, lambda_ridge, TRUE, keep_example)
}
