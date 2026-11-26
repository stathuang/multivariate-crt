# E02: ridge-gain power curve using full randomized covariate matrices.

power_coefficients <- function(s, t, setting_index) {
  set.seed(91000L + setting_index)
  q <- s + 1L
  list(
    Bcoef = matrix(rnorm(q * t, sd = 0.7 / sqrt(q)), q, t),
    gamma = rnorm(q, sd = 1 / sqrt(q)),
    signal = rep(1 / sqrt(t), t)
  )
}

power_data <- function(cfg, coefs, theta) {
  XS <- matrix(rnorm(cfg$n * cfg$s), cfg$n, cfg$s)
  D <- cbind(1, XS)
  E <- matrix(rnorm(cfg$n * cfg$t), cfg$n, cfg$t)
  XT <- D %*% coefs$Bcoef + E
  Y <- drop(
    D %*% coefs$gamma +
      theta * E %*% coefs$signal +
      rnorm(cfg$n)
  )
  PR <- residual_basis(D)
  Bhat <- qr.coef(qr(D), XT)
  G <- PR %*% XT
  rY <- drop(Y - D %*% qr.coef(qr(D), Y))
  list(D = D, PR = PR, Bhat = Bhat, G = G, rY = rY, X = cbind(XS, XT), Y = Y, E = E)
}

ridge_stat <- function(G, rd, lambda = 1) {
  v <- drop(crossprod(G, rd))
  d <- nrow(G)
  t_dim <- ncol(G)
  if (t_dim <= d) {
    spd_quadratic_form(
      crossprod(G) + diag(lambda, t_dim),
      v
    )
  } else {
    H <- tcrossprod(G)
    sum(rd^2) - lambda * spd_quadratic_form(
      H + diag(lambda, d),
      rd
    )
  }
}

ridge_copy_stats <- function(
    Gobs,
    D,
    PR,
    Bhat,
    Bcoef,
    rd,
    method,
    denom = NULL,
    M,
    lambda = 1) {
  d <- nrow(Gobs)
  t_dim <- ncol(Gobs)
  out <- numeric(M)

  if (method == "MVN-CRT" && d >= t_dim) {
    C <- chol(crossprod(Gobs))
    for (m in seq_len(M)) {
      Gdraw <- stiefel_matrix(d, t_dim) %*% C
      Xtilde <- D %*% Bhat + t(PR) %*% Gdraw
      G <- PR %*% Xtilde
      out[m] <- ridge_stat(G, rd, lambda)
    }
    return(out)
  }

  if (method == "MVN-CRT") {
    for (m in seq_len(M)) {
      Gdraw <- haar_orthogonal(d) %*% Gobs
      Xtilde <- D %*% Bhat + t(PR) %*% Gdraw
      G <- PR %*% Xtilde
      out[m] <- ridge_stat(G, rd, lambda)
    }
    return(out)
  }

  if (method == "Oracle CRT") {
    for (m in seq_len(M)) {
      Xtilde <- (
        D %*% Bcoef +
          matrix(rnorm(nrow(D) * t_dim), nrow(D), t_dim)
      )
      G <- PR %*% Xtilde
      out[m] <- ridge_stat(G, rd, lambda)
    }
    return(out)
  }

  Sigmahat <- crossprod(Gobs) / denom
  L <- psd_factor(Sigmahat)
  rank_sigma <- ncol(L)
  for (m in seq_len(M)) {
    noise <- if (rank_sigma == 0L) {
      matrix(0, nrow(D), nrow(L))
    } else {
      matrix(
        rnorm(nrow(D) * rank_sigma),
        nrow(D),
        rank_sigma
      ) %*% t(L)
    }
    Xtilde <- D %*% Bhat + noise
    G <- PR %*% Xtilde
    out[m] <- ridge_stat(G, rd, lambda)
  }
  out
}

# The four methods run sequentially on this dataset, in the original order.
power_replication <- function(b, cfg, coefs, theta, setting_index, theta_index,
                               seed_base, M, lambda_ridge, keep_example = FALSE) {
  set.seed(seed_base + 1000000L * setting_index + 10000L * theta_index + b)
  dat <- power_data(cfg, coefs, theta)
  rd <- drop(dat$PR %*% dat$rY)
  T0 <- ridge_stat(dat$G, rd, lambda_ridge)
  randomized <- list(
    "Oracle CRT" = ridge_copy_stats(dat$G, dat$D, dat$PR, dat$Bhat, coefs$Bcoef,
      rd, "Oracle CRT", M = M, lambda = lambda_ridge),
    "Plug-in MLE CRT" = ridge_copy_stats(dat$G, dat$D, dat$PR, dat$Bhat, coefs$Bcoef,
      rd, "Plug-in MLE CRT", denom = cfg$n, M = M, lambda = lambda_ridge),
    "Plug-in CRT (d.f. adjustment)" = ridge_copy_stats(dat$G, dat$D, dat$PR, dat$Bhat,
      coefs$Bcoef, rd, "Plug-in CRT (d.f. adjustment)", denom = cfg$d, M = M, lambda = lambda_ridge),
    "MVN-CRT" = ridge_copy_stats(dat$G, dat$D, dat$PR, dat$Bhat, coefs$Bcoef,
      rd, "MVN-CRT", M = M, lambda = lambda_ridge)
  )
  pvals <- vapply(randomized, function(values) mc_pvalue(T0, values), 0)
  list(pvalue = pvals, statistic = T0,
       diagnostics = if (keep_example) list(data = dat, residual_response = rd,
                          randomized_statistics = randomized) else NULL)
}
