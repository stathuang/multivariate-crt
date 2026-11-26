# Shared mathematical helpers for E01/E02, retained from the corrected benchmark.
# D contains the intercept and conditioning columns. PR has orthonormal rows
# spanning their residual space. G = PR %*% XT represents the target residuals.
# Haar/Stiefel draws keep the approved QR diagonal-sign correction.

residual_basis <- function(D) {
  qrD <- qr(D)
  q <- ncol(D)
  t(qr.Q(qrD, complete = TRUE)[, (q + 1):nrow(D), drop = FALSE])
}

psd_factor <- function(Sigma, tol = 1e-10) {
  Sigma <- (Sigma + t(Sigma)) / 2
  eig <- eigen(Sigma, symmetric = TRUE)
  scale <- max(1, max(abs(eig$values)))
  if (any(eig$values < -tol * scale)) {
    stop("The fitted covariance matrix is not positive semidefinite.")
  }
  keep <- eig$values > tol * scale
  if (!any(keep)) {
    return(matrix(0, nrow(Sigma), 0L))
  }
  sweep(
    eig$vectors[, keep, drop = FALSE],
    2,
    sqrt(pmax(eig$values[keep], 0)),
    `*`
  )
}

spd_quadratic_form <- function(M, v) {
  z <- backsolve(chol(M), v, transpose = TRUE)
  sum(z^2)
}

haar_orthogonal <- function(d) {
  z <- matrix(rnorm(d * d), d, d)
  qrz <- qr(z)
  U <- qr.Q(qrz)
  R <- qr.R(qrz)
  # QR's sign convention means its raw Q need not be Haar-distributed.
  # Flip each column of Q by sign(R[j,j]); the corresponding R then has
  # positive diagonal. For iid standard Gaussian input, this produces
  # Haar-distributed orthogonal matrices.
  # No additional random signs are needed.
  signs <- sign(diag(R))
  # Avoid zeroing a column; this is not a remedy for rank deficiency.
  # Exact zeros occur with probability zero under the full-rank Gaussian model.
  signs[signs == 0] <- 1
  sweep(U, 2, signs, `*`)
}

stiefel_matrix <- function(d, t_dim) {
  z <- matrix(rnorm(d * t_dim), d, t_dim)
  qrz <- qr(z)
  U <- qr.Q(qrz)
  R <- qr.R(qrz)
  # QR's sign convention means its raw Q need not be Haar-distributed.
  # Flip each column of Q by sign(R[j,j]); the corresponding R then has
  # positive diagonal. For iid standard Gaussian input, this produces
  # uniform Stiefel frames, provided d >= t_dim.
  # No additional random signs are needed.
  signs <- sign(diag(R))
  # Avoid zeroing a column; this is not a remedy for rank deficiency.
  # Exact zeros occur with probability zero under the full-rank Gaussian model.
  signs[signs == 0] <- 1
  sweep(U, 2, signs, `*`)
}

mc_pvalue <- function(T0, Tm) {
  (1 + sum(Tm >= T0)) / (length(Tm) + 1)
}
