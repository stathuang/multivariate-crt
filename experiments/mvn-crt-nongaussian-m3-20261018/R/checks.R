# Deterministic implementation checks plus separately labelled sampling diagnostics.
bank <- readRDS("seeds.rds")
stopifnot(identical(bank$config, config))
set_stream(bank$check_seed)
cells <- cells_for(config)
checks <- list()
add_check <- function(name, error, tolerance = 1e-10, detail = "") {
  passed <- is.finite(error) && error <= tolerance
  checks[[length(checks) + 1L]] <<- data.frame(name = name, error = error,
     tolerance = tolerance, passed = passed, detail = detail)
}
relative_error <- function(a, b) max(abs(a - b)) / max(1, abs(b))
examples <- list()
for (cid in seq_len(nrow(cells))) {
  n <- cells$n[cid]
  distr <- cells$distribution[cid]
  dat <- generate_data(n, distr, config)
  ctx <- prepare_statistic(dat$Z, dat$U, dat$Y, config$lambda)
  tag <- paste0(distr, "-n", n, ":")
  add_check(paste0(tag, "rank and residual dimension"),
    as.numeric(ctx$qD$rank != 13L || ctx$d != n - 13L), 0)
  add_check(paste0(tag, "row scales"),
    relative_error(dat$e, sweep(dat$G, 1, dat$scale, "*")))
  if (distr != "Gaussian") {
    nu <- if (distr == "t5") 5 else 3
    add_check(paste0(tag, "t standardization"),
       relative_error(dat$scale^2 * dat$V, rep(nu - 2, n)))
    add_check(paste0(tag, "one separate scale per row"),
       as.numeric(length(dat$V) != n || length(unique(dat$V)) != n), 0)
  }
  add_check(paste0(tag, "response formula"), relative_error(dat$Y,
    rowSums(dat$Z) / sqrt(12) + (dat$Z[, 1]^2 - 1) / sqrt(2) + dat$epsilon))
  add_check(paste0(tag, "target conditional mean"), relative_error(dat$U - dat$e,
    (dat$Z[, 1:8] + dat$Z[, 2:9]) / sqrt(2)))
  add_check(paste0(tag, "no quadratic column in background fit"),
     as.numeric(!identical(unname(ctx$D), unname(cbind(1, dat$Z)))), 0)
  PR <- residual_basis(ctx$D)
  Gobs <- PR %*% dat$U
  rd <- as.vector(PR %*% dat$Y)
  add_check(paste0(tag, "observed residual-coordinate statistic"),
    relative_error(ctx$T0, reference_ridge(Gobs, rd)))
  direct <- drop(t(ctx$ry) %*% ctx$r %*% solve(crossprod(ctx$r) + diag(8)) %*%
                 t(ctx$r) %*% ctx$ry)
  add_check(paste0(tag, "observed original matrix formula"), relative_error(ctx$T0, direct))
  C <- chol(crossprod(Gobs))
  H <- C %*% solve(crossprod(Gobs) + diag(8)) %*% t(C)
  eigH <- eigen(H, symmetric = TRUE)
  add_check(paste0(tag, "spectral ridge weights"),
      relative_error(eigH$values, ctx$weights))
  Bhat <- qr.coef(ctx$qD, dat$U)
  copy_errors <- matrix(0, 3, 5)
  for (m in 1:3) {
    Q <- stiefel_matrix(ctx$d, 8L)
    Gdraw <- Q %*% C
    Udraw <- ctx$D %*% Bhat + t(PR) %*% Gdraw
    rdraw <- qr.resid(ctx$qD, Udraw)
    Tdirect <- drop(t(ctx$ry) %*% rdraw %*% solve(crossprod(rdraw) + diag(8)) %*%
                    t(rdraw) %*% ctx$ry)
    Tref <- reference_ridge(PR %*% Udraw, rd)
    sphere_coordinates <- as.vector(crossprod(eigH$vectors, crossprod(Q, rd)))
    Tspectral <- sum(ctx$weights * sphere_coordinates^2)
    copy_errors[m, ] <- c(
      relative_error(crossprod(ctx$D, Udraw), crossprod(ctx$D, dat$U)),
      relative_error(crossprod(Udraw), crossprod(dat$U)),
      relative_error(crossprod(rdraw), ctx$Gram),
      relative_error(Tdirect, Tref), relative_error(Tdirect, Tspectral))
  }
  for (j in 1:5) add_check(paste0(tag, c("D'U preserved", "U'U preserved",
       "residual Gram preserved", "copy score vs reference", "copy score vs spectral formula")[j]),
       max(copy_errors[, j]))
  examples[[cid]] <- list(cell = cells[cid, ], data = dat, T0 = ctx$T0)
}
# Check that the locally retained reference applies the required QR sign correction.
before <- .Random.seed
z <- matrix(rnorm(37 * 8), 37, 8)
qz <- qr(z)
expected <- sweep(qr.Q(qz), 2L, sign(diag(qr.R(qz))), "*")
set_stream(before)
actual <- stiefel_matrix(37, 8)
add_check("QR sign correction retained", relative_error(actual, expected))
add_check("Stiefel orthonormality", relative_error(crossprod(actual), diag(8)))
add_check("positive corrected QR diagonal", as.numeric(any(diag(crossprod(actual, z)) <= 0)), 0)
add_check("p-value: all ties", abs(rank_pvalue(1, rep(1, 199)) - 1), 0)
add_check("p-value: minimum", abs(rank_pvalue(1, rep(0, 199)) - 0.005), 0)
at_boundary <- rank_pvalue(1, c(rep(2, 9), rep(0, 190)))
above_boundary <- rank_pvalue(1, c(rep(2, 10), rep(0, 189)))
add_check("p-value: 0.05 rejected", as.numeric(at_boundary != 0.05 || !(at_boundary <= 0.05)), 0)
add_check("p-value: 0.055 retained", as.numeric(above_boundary != 0.055 || above_boundary <= 0.05), 0)
add_check("p-value: greater-than-or-equal tie convention", abs(rank_pvalue(1, c(1, 0, 2)) - 0.75), 0)

# The following Monte Carlo diagnostics are descriptive, not tuning rules.
dat <- examples[[1]]$data
ctx <- prepare_statistic(dat$Z, dat$U, dat$Y)
Ndiag <- 50000L
quick <- draw_ridge_scores(ctx, Ndiag)
expected_mean <- ctx$response_ss * sum(ctx$weights) / ctx$d
expected_var <- 2 * ctx$response_ss^2 / (ctx$d * (ctx$d + 2)) *
  (sum(ctx$weights^2) - sum(ctx$weights)^2 / ctx$d)
PR <- residual_basis(ctx$D)
Gobs <- PR %*% dat$U
rd <- as.vector(PR %*% dat$Y)
C <- chol(crossprod(Gobs))
explicit <- replicate(1500L, reference_ridge(stiefel_matrix(ctx$d, 8) %*% C, rd))
diagnostics <- data.frame(
  fast_draws = Ndiag, explicit_draws = length(explicit),
  expected_mean = expected_mean, fast_mean = mean(quick), explicit_mean = mean(explicit),
  fast_mean_standardized_error = (mean(quick) - expected_mean) / sqrt(expected_var / Ndiag),
  explicit_mean_standardized_error = (mean(explicit) - expected_mean) /
    sqrt(expected_var / length(explicit)),
  expected_variance = expected_var, fast_variance = var(quick), explicit_variance = var(explicit),
  ks_distance = unname(ks.test(quick, explicit)$statistic))
write.csv(diagnostics, "checks/sampling_diagnostics.csv", row.names = FALSE)
saveRDS(list(fast_scores = quick, explicit_scores = explicit, example = examples[[1]]),
        "checks/sampling_diagnostics.rds")
saveRDS(examples, "checks/examples.rds")

check_table <- do.call(rbind, checks)
write.csv(check_table, "checks/invariants.csv", row.names = FALSE)
cat("Implementation checks:", nrow(check_table), "; passed:", sum(check_table$passed), "\n")
print(diagnostics)
if (!all(check_table$passed)) {
  print(check_table[!check_table$passed, ])
  stop("Implementation checks failed; formal simulation not allowed")
}
saveRDS(list(passed_at = singapore_time(), checks = nrow(check_table),
  code_md5 = tools::md5sum(c("config.R", "R/core.R", "R/reference.R", "R/checks.R", "run.R"))),
  "checks/PASS.rds")
