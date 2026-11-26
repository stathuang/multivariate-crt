# Logistic response with beta and uniform noise shared across signal levels.
logistic_response <- function(X, target, theta, beta, noise_uniform) {
  beta[target] <- beta[target] * theta
  predictor <- drop(X %*% beta)
  as.numeric(noise_uniform < stats::plogis(predictor))
}

# Fixed nonlinear binary model: transform each of the first twenty columns.
# theta scales columns 1:10, while the terms from 11:20 have weight 1.
# E09/E10 test columns 1:8 conditional on the rest, including columns 9:10.
# beta is unused; changing the test's target does not move this signal block.
nonlinear_binary_response <- function(X, theta, noise_uniform) {
  if (ncol(X) < 20L) stop("The nonlinear binary design requires twenty covariates")
  z <- X[, 1:20, drop = FALSE]
  transformed <- ifelse(z >= stats::qnorm(1/3) & z <= stats::qnorm(2/3), -.5, .25)
  predictor <- drop(transformed %*% c(rep(theta, 10), rep(1, 10)))
  as.numeric(noise_uniform < stats::plogis(predictor))
}
