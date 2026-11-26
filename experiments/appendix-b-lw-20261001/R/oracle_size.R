# E01: keep explicit oracle, plug-in and MVN randomization matrices.

size_coefficients <- function(s, t, setting_index) {
  set.seed(91000L + setting_index)
  q <- s + 1L
  list(
    Bcoef = matrix(rnorm(q * t, sd = 0.7 / sqrt(q)), q, t),
    gamma = rnorm(q, sd = 1 / sqrt(q))
  )
}

size_data <- function(cfg, coefs) {
  XS <- matrix(rnorm(cfg$n * cfg$s), cfg$n, cfg$s)
  D <- cbind(1, XS)
  XT <- D %*% coefs$Bcoef + matrix(rnorm(cfg$n * cfg$t), cfg$n, cfg$t)
  Y <- drop(D %*% coefs$gamma + rnorm(cfg$n))
  PR <- residual_basis(D)
  Bhat <- qr.coef(qr(D), XT)
  G <- PR %*% XT
  rY <- drop(Y - D %*% qr.coef(qr(D), Y))
  list(D = D, PR = PR, Bhat = Bhat, G = G, rY = rY, X = cbind(XS, XT), Y = Y)
}

statistic_context <- function(PR, rY) {
  PRt <- t(PR)
  list(
    rd = drop(PR %*% rY),
    H2 = crossprod(PRt, PRt * as.vector(rY))
  )
}

stat_values <- function(G, context, lambda = 1) {
  rd <- context$rd
  H2 <- context$H2
  raw_vec <- drop(crossprod(G, rd))
  raw <- sum(raw_vec^2)

  poly_vec <- vapply(seq_len(ncol(G)), function(j) {
    gj <- G[, j, drop = FALSE]
    drop(t(gj) %*% H2 %*% gj)
  }, numeric(1))
  poly <- raw + 0.25 * sum(poly_vec^2)

  d <- nrow(G)
  t_dim <- ncol(G)
  if (t_dim <= d) {
    ridge <- spd_quadratic_form(
      crossprod(G) + diag(lambda, t_dim),
      raw_vec
    )
  } else {
    H <- tcrossprod(G)
    ridge <- sum(rd^2) - lambda * spd_quadratic_form(
      H + diag(lambda, d),
      rd
    )
  }

  c("raw score" = raw, "poly2 score" = poly, "ridge gain" = ridge)
}

oracle_stats <- function(D, PR, Bcoef, context, M, lambda = 1) {
  n <- nrow(D)
  t_dim <- ncol(Bcoef)
  out <- matrix(
    NA_real_,
    M,
    3L,
    dimnames = list(NULL, c("raw score", "poly2 score", "ridge gain"))
  )

  for (m in seq_len(M)) {
    Xtilde <- D %*% Bcoef + matrix(rnorm(n * t_dim), n, t_dim)
    G <- PR %*% Xtilde
    out[m, ] <- stat_values(G, context, lambda)
  }
  out
}

mvn_stats <- function(
    Gobs,
    D,
    PR,
    Bhat,
    context,
    M,
    lambda = 1) {
  d <- nrow(Gobs)
  t_dim <- ncol(Gobs)
  out <- matrix(
    NA_real_,
    M,
    3L,
    dimnames = list(NULL, c("raw score", "poly2 score", "ridge gain"))
  )

  if (d >= t_dim) {
    C <- chol(crossprod(Gobs))
    for (m in seq_len(M)) {
      Gdraw <- stiefel_matrix(d, t_dim) %*% C
      Xtilde <- D %*% Bhat + t(PR) %*% Gdraw
      G <- PR %*% Xtilde
      out[m, ] <- stat_values(G, context, lambda)
    }
    return(out)
  }

  for (m in seq_len(M)) {
    Gdraw <- haar_orthogonal(d) %*% Gobs
    Xtilde <- D %*% Bhat + t(PR) %*% Gdraw
    G <- PR %*% Xtilde
    out[m, ] <- stat_values(G, context, lambda)
  }
  out
}

plugin_stats <- function(
    Gobs,
    D,
    PR,
    Bhat,
    context,
    denom,
    M,
    lambda = 1) {
  n <- nrow(D)
  Sigmahat <- crossprod(Gobs) / denom
  L <- psd_factor(Sigmahat)
  rank_sigma <- ncol(L)
  out <- matrix(
    NA_real_,
    M,
    3L,
    dimnames = list(NULL, c("raw score", "poly2 score", "ridge gain"))
  )

  for (m in seq_len(M)) {
    noise <- if (rank_sigma == 0L) {
      matrix(0, n, nrow(L))
    } else {
      matrix(rnorm(n * rank_sigma), n, rank_sigma) %*% t(L)
    }
    Xtilde <- D %*% Bhat + noise
    G <- PR %*% Xtilde
    out[m, ] <- stat_values(G, context, lambda)
  }
  out
}

# One seeded replicate returns p-values and observed statistics. The first
# replicate also retains the generated data and all randomized statistics.
size_replication <- function(b, cfg, coefs, setting_index, seed_base, M,
                              lambda_ridge, keep_example = FALSE) {
  stat_names <- c("raw score", "poly2 score", "ridge gain")
  method_names <- c("Oracle CRT", "Plug-in MLE CRT", "Plug-in CRT (d.f. adjustment)", "MVN-CRT")

  replication_seed <- seed_base + 100000L * setting_index + b
  set.seed(replication_seed)

  dat <- size_data(cfg, coefs)
  context <- statistic_context(dat$PR, dat$rY)
  T0 <- stat_values(dat$G, context, lambda = lambda_ridge)

  randomized_stats <- list(
    "Oracle CRT" = oracle_stats(
      dat$D,
      dat$PR,
      coefs$Bcoef,
      context,
      M = M,
      lambda = lambda_ridge
    ),
    "Plug-in MLE CRT" = plugin_stats(
      dat$G,
      dat$D,
      dat$PR,
      dat$Bhat,
      context,
      denom = cfg$n,
      M = M,
      lambda = lambda_ridge
    ),
    "Plug-in CRT (d.f. adjustment)" = plugin_stats(
      dat$G,
      dat$D,
      dat$PR,
      dat$Bhat,
      context,
      denom = cfg$d,
      M = M,
      lambda = lambda_ridge
    ),
    "MVN-CRT" = mvn_stats(
      dat$G,
      dat$D,
      dat$PR,
      dat$Bhat,
      context,
      M = M,
      lambda = lambda_ridge
    )
  )

  pvals <- matrix(
    NA_real_,
    length(stat_names),
    length(method_names),
    dimnames = list(stat_names, method_names)
  )
  for (method in method_names) {
    for (statistic in stat_names) {
      pvals[statistic, method] <- mc_pvalue(
        T0[[statistic]],
        randomized_stats[[method]][, statistic]
      )
    }
  }


  list(pvalue = pvals, statistic = T0,
       diagnostics = if (keep_example) list(data = dat, context = context,
                          randomized_statistics = randomized_stats) else NULL)
}
