# Mathematical equivalence checks, original-output parity, and failure handling.
root <- "experiments/appendix-b-lw-20261001"
source(file.path(root, "R/run_tools.R"))
source(file.path(root, "R/lw_run_tools.R"))
loaded <- load_lw_sources(root, globalenv(), "E01.R")
reference <- new.env()
source(file.path(root, "reference/dcrt-author-functions.R"), local = reference)
check <- function(condition, label) {
  if (!isTRUE(condition)) stop(label)
  cat("PASS:", label, "\n")
}
near <- function(a, b) isTRUE(all.equal(a, b, tolerance = 1e-10, check.attributes = FALSE))
settings <- data.frame(n = c(60L, 80L, 80L), s = c(45L, 50L, 60L), t = c(10L, 10L, 40L), d = c(14L, 29L, 19L))
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
for (i in 1:3) {
  cfg <- settings[i, ]
  coefs <- size_coefficients(cfg$s, cfg$t, i)
  old <- size_replication(1L, cfg, coefs, i, 20260429L, 19L, 1, TRUE)
  old_rng <- .Random.seed
  new <- size_replication_lw(1L, cfg, coefs, i, 20260429L, 19L, 1, TRUE)
  check(identical(new$pvalue[, 1:4], old$pvalue), paste("E01 original p-values", i))
  check(identical(new$diagnostics$rng_after_original, old_rng), paste("E01 original RNG", i))
  check(identical(new$diagnostics$randomized_statistics[1:4], old$diagnostics$randomized_statistics), paste("E01 candidate statistics", i))
  X <- old$diagnostics$data$X
  fit <- lw_covariance(X)
  expected_unadjusted <- reference$linshrink_cov(X, k = 0, normalize = FALSE)
  check(near(fit$unadjusted, expected_unadjusted), paste("author unadjusted shrinkage covariance", i))
  invisible(capture.output(expected <- reference$linshrink_cov(X, k = 0, normalize = TRUE)))
  check(near(fit$covariance, expected), paste("author covariance incl. centering", i))
  Xc <- sweep(X, 2, colMeans(X), "-")
  invisible(capture.output(expected <- reference$linshrink_cov(Xc, k = 1, normalize = TRUE)))
  check(near(fit$covariance, expected), paste("D.10 k=1 on centered input", i))
  check(near(1 / diag(solve(fit$covariance)), fit$residual_variance), paste("post-adjustment conditional variances", i))
  conditional <- lw_conditional(X, cfg$s)
  target <- seq.int(cfg$s + 1L, ncol(X))
  check(near(conditional$covariance, solve(solve(fit$covariance)[target, target])), paste("joint Schur complement", i))
  check(near(crossprod(conditional$root), conditional$covariance), paste("draw root orientation", i))
  # Independent scalar regression expressions and exact zero residual mean.
  for (j in c(1L, ncol(X))) {
    beta <- solve(fit$unadjusted[-j, -j], fit$unadjusted[-j, j])
    check(near(mean((Xc[, j] - Xc[, -j] %*% beta)^2), fit$residual_variance[j]), paste("nodewise residual", i, j))
  }
  coefs <- power_coefficients(cfg$s, cfg$t, i)
  old <- power_replication(1L, cfg, coefs, .5, i, 3L, 20260429L, 19L, 1, TRUE)
  old_rng <- .Random.seed
  new <- power_replication_lw(1L, cfg, coefs, .5, i, 3L, 20260429L, 19L, 1, TRUE)
  check(identical(new$pvalue[1:4], old$pvalue), paste("E02 original p-values", i))
  check(identical(new$diagnostics$rng_after_original, old_rng), paste("E02 original RNG", i))
  check(identical(new$diagnostics$data, old$diagnostics$data), paste("same paired data", i))
  check(new$diagnostics$lw_status == "ok" && is.finite(new$pvalue[5]), paste("LW valid output", i))
}
set.seed(42)
X <- matrix(rnorm(180), 30, 6)
check(near(lw_covariance(X)$covariance, lw_covariance(sweep(X, 2, 1:6, "+"))$covariance), "intercept/translation invariance")
zero <- try(lw_covariance(matrix(0, 10, 3)), silent = TRUE)
check(inherits(zero, "try-error"), "undefined zero-variance input fails explicitly")
saved_copy_stats <- lw_copy_stats
lw_copy_stats <- function(...) stop("Deliberate validation failure")
failed <- append_lw(old, settings[3, ], 19L, 1, TRUE, FALSE)
check(identical(failed$pvalue[1:4], old$pvalue), "LW failure preserves other method outputs")
check(is.na(failed$pvalue[5]) && failed$diagnostics$lw_status == "failed", "LW failure remains missing")
rows <- data.frame(method = names(failed$pvalue), status = "ok", error = "", warnings = "")
rows <- mark_lw_failures(rows, list(status = "ok", value = failed, warnings = ""))
check(identical(rows$status, c(rep("ok", 4), "failed")), "method-specific failure denominator")
lw_copy_stats <- saved_copy_stats
cat("All checks passed.\n")
