# Core simulation code. This file only defines functions.
cells_for <- function(cfg) {
  data.frame(cell_id = seq_len(12L),
             distribution = rep(cfg$distributions, each = 4L),
             n = rep(cfg$sample_sizes, 3L), stringsAsFactors = FALSE)
}

set_stream <- function(seed) assign(".Random.seed", as.integer(seed), .GlobalEnv)

generate_data <- function(n, distribution, cfg) {
  Z <- matrix(rnorm(n * cfg$s), n, cfg$s)
  Sigma <- cfg$rho^abs(outer(seq_len(cfg$t), seq_len(cfg$t), "-"))
  G <- matrix(rnorm(n * cfg$t), n, cfg$t) %*% chol(Sigma)
  V <- rep(NA_real_, n)
  scale <- rep(1, n)
  if (distribution != "Gaussian") {
    nu <- switch(distribution, t5 = 5, t3 = 3, stop("Unknown distribution"))
    V <- rchisq(n, nu)
    scale <- sqrt((nu - 2) / V)
  }
  e <- sweep(G, 1L, scale, "*")
  mu <- (Z[, seq_len(cfg$t)] + Z[, seq_len(cfg$t) + 1L]) / sqrt(2)
  U <- mu + e
  epsilon <- rnorm(n)
  Y <- rowSums(Z) / sqrt(cfg$s) + (Z[, 1L]^2 - 1) / sqrt(2) + epsilon
  list(Z = Z, U = U, Y = Y, G = G, e = e, V = V,
       scale = scale, epsilon = epsilon, Sigma = Sigma)
}

# The test receives only observed Z, U, and Y, never population parameters.
prepare_statistic <- function(Z, U, Y, lambda = 1) {
  D <- cbind(1, Z)
  qD <- qr(D)
  if (qD$rank != ncol(D)) stop("Background design is rank deficient")
  d <- nrow(D) - qD$rank
  if (d <= ncol(U)) stop("This implementation requires d > t")
  r <- qr.resid(qD, U)
  ry <- as.vector(qr.resid(qD, Y))
  Gram <- crossprod(r)
  eigenvalues <- eigen(Gram, symmetric = TRUE, only.values = TRUE)$values
  if (min(eigenvalues) <= 0) stop("Residual Gram matrix is not positive definite")
  cvec <- as.vector(crossprod(r, ry))
  z <- backsolve(chol(Gram + diag(lambda, ncol(U))), cvec, transpose = TRUE)
  T0 <- sum(z^2)
  list(D = D, qD = qD, r = r, ry = ry, Gram = Gram, d = d,
       eigenvalues = eigenvalues, weights = eigenvalues / (eigenvalues + lambda),
       response_ss = sum(ry^2), T0 = T0)
}

# Exact pushforward of independent uniform Stiefel residual copies for this
# statistic. The distribution is evaluated directly, without an asymptotic approximation.
draw_ridge_scores <- function(context, M) {
  t <- length(context$weights)
  normals2 <- matrix(rnorm(M * t)^2, M, t)
  other_ss <- rchisq(M, context$d - t)
  context$response_ss * as.vector(normals2 %*% context$weights) /
    (rowSums(normals2) + other_ss)
}

rank_pvalue <- function(T0, Tm) (1L + sum(Tm >= T0)) / (length(Tm) + 1L)

one_replication <- function(seed, n, distribution, cfg) {
  set_stream(seed)
  dat <- generate_data(n, distribution, cfg)
  context <- prepare_statistic(dat$Z, dat$U, dat$Y, cfg$lambda)
  Tm <- draw_ridge_scores(context, cfg$M)
  if (any(!is.finite(c(context$T0, Tm)))) stop("Non-finite statistic")
  p <- rank_pvalue(context$T0, Tm)
  list(T0 = context$T0, Tm = Tm, p_value = p,
       exceedances = sum(Tm >= context$T0), ties = sum(Tm == context$T0),
       reject = p <= cfg$alpha, d = context$d, rank_D = context$qD$rank,
       min_eigenvalue = min(context$eigenvalues),
       max_eigenvalue = max(context$eigenvalues), response_ss = context$response_ss)
}

wilson_interval <- function(k, N, confidence = 0.95) {
  if (N == 0) return(c(NA_real_, NA_real_))
  z <- qnorm(1 - (1 - confidence) / 2)
  p <- k / N
  denom <- 1 + z^2 / N
  centre <- (p + z^2 / (2 * N)) / denom
  half <- z * sqrt(p * (1 - p) / N + z^2 / (4 * N^2)) / denom
  c(centre - half, centre + half)
}

singapore_time <- function() format(Sys.time(), tz = "Asia/Singapore", usetz = TRUE)
