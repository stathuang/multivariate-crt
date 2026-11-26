# Local reference functions copied from the handoff's oracle_common.R.
# Old experiment files are only read, never sourced or executed.
# Source SHA-256: 1ba74625e08e0fc8c220a1f7a0cfd90ba7c51e3a5a5e6b53494c2ebd95946567
residual_basis <- function(D) {
  qrD <- qr(D)
  q <- ncol(D)
  t(qr.Q(qrD, complete = TRUE)[, (q + 1):nrow(D), drop = FALSE])
}

spd_quadratic_form <- function(M, v) {
  z <- backsolve(chol(M), v, transpose = TRUE)
  sum(z^2)
}

stiefel_matrix <- function(d, t_dim) {
  z <- matrix(rnorm(d * t_dim), d, t_dim)
  qrz <- qr(z)
  U <- qr.Q(qrz)
  R <- qr.R(qrz)
  signs <- sign(diag(R))
  signs[signs == 0] <- 1
  sweep(U, 2, signs, `*`)
}

mc_pvalue <- function(T0, Tm) (1 + sum(Tm >= T0)) / (length(Tm) + 1)

# Ridge-only algebra from oracle_size.R:stat_values, omitting its other scores.
# Source SHA-256: 94ae6d125abff66626be2b2b5c339a8fcc975e3f0c99c5f7d7a79b8dc88fe8ce
reference_ridge <- function(G, rd, lambda = 1) {
  spd_quadratic_form(crossprod(G) + diag(lambda, ncol(G)),
                     as.vector(crossprod(G, rd)))
}
