# Shared Gaussian design: permuted band precision, unit marginal variances, random beta.
# Return X, the true graph, beta, and the shared response noise.
gaussian_design <- function(n, p, bandwidth, precision_offdiag) {
  if (n < 2L || p < 3L || bandwidth < 0L || bandwidth >= p)
    stop("Invalid Gaussian dimensions or graph bandwidth")
  distance <- abs(outer(seq_len(p), seq_len(p), `-`))
  graph <- 1L * (distance > 0 & distance <= bandwidth)
  precision <- diag(p) + precision_offdiag * graph
  covariance <- stats::cov2cor(solve(precision))

  # Preserve the reconstruction's draw order: permutation, X, beta, noise.
  permutation <- sample.int(p)
  covariance <- covariance[permutation, permutation, drop = FALSE]
  X <- MASS::mvrnorm(n, rep(0, p), covariance)
  graph <- graph[permutation, permutation, drop = FALSE]
  beta <- numeric(p)
  active <- seq_len(min(p, 20L))
  # Nonlinear responses ignore beta; retain its draw to preserve the following RNG sequence.
  beta[active] <- stats::runif(length(active), 1, 2) / sqrt(length(active))
  noise_uniform <- stats::runif(n)

  list(X = X, graph = graph, beta = beta, noise_uniform = noise_uniform,
       covariance = covariance, permutation = permutation)
}

# Share beta and noise across signal levels within each replicate.
linear_response <- function(X, target, theta, beta, noise_uniform) {
  beta[target] <- beta[target] * theta
  predictor <- drop(X %*% beta)
  predictor + stats::qnorm(noise_uniform)
}

# Seven nonlinear terms for each ten-variable block, in the accepted order.
nonlinear_basis <- function(X) {
  cbind(X[, 1]^2 / 2, 1 / (1 + X[, 2]^2), cos(pi * X[, 3]),
        sin(pi * X[, 4]), X[, 5] * X[, 6], sin(X[, 7] * X[, 8]),
        sin(pi * X[, 9]) / (4 + X[, 10]^2))
}

# theta scales all seven terms from columns 1:10; the terms from 11:20 have weight 1.
# E05/E06 test columns 1:8 conditional on the rest, including 9:10 in the scaled block.
# The response uses fixed column positions, independent of beta and the test's target setting.
nonlinear_response <- function(X, theta, noise_uniform) {
  if (ncol(X) < 20L) stop("The nonlinear design requires twenty covariates")
  signal_terms <- nonlinear_basis(X[, 1:10, drop = FALSE])
  background_terms <- nonlinear_basis(X[, 11:20, drop = FALSE])
  predictor <- drop(cbind(signal_terms, background_terms) %*%
                      c(rep(theta, 7), rep(1, 7)))
  predictor + stats::qnorm(noise_uniform)
}
