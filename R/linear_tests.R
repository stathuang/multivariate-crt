# Statistical input checks stay with the test, not with file management.
check_test_inputs <- function(X, Y, target) {
  .ms_matrix(X)
  .ms_target(target, ncol(X))
  if (!is.numeric(Y) || length(Y) != nrow(X) || any(!is.finite(Y)))
    stop("Y must be a finite numeric response of length nrow(X)")
}

# LM-SST sums squared target t statistics; LM-SSR uses negative deviance.
# Use the same glm arithmetic as the accepted reconstruction.
linear_statistic <- function(X, Y, target, statistic) {
  fit <- stats::glm(Y ~ X, family = "gaussian")
  if (statistic == "LM-SSR") return(-fit$deviance)
  coefficients <- summary(fit)$coefficients
  if (fit$rank < ncol(X) + 1L || fit$df.residual <= 0L)
    stop("SST requires an identifiable fit with positive residual degrees of freedom")
  sum(coefficients[1L + target, 3L]^2)
}

# Return the p-value, observed statistic, and the small randomized-statistic vector.
gcrt_linear_test <- function(X, Y, target, copies, statistic, p_rule = "discrete") {
  check_test_inputs(X, Y, target)
  statistic <- match.arg(statistic, c("LM-SST", "LM-SSR"))
  if (!is.list(copies) || !is.list(copies$X_copies) ||
      !identical(copies$reduced, TRUE) ||
      !identical(as.integer(copies$target), as.integer(target)))
    stop("Expected reduced sampler output in target order")

  candidates <- c(list(X[, target, drop = FALSE]), copies$X_copies)
  values <- numeric(length(candidates))
  for (m in seq_along(candidates)) {
    candidate <- candidates[[m]]
    .ms_matrix(candidate, "copy")
    if (!identical(dim(candidate), c(nrow(X), length(target))))
      stop("Copy dimensions differ from the target matrix")
    design <- X
    design[, target] <- candidate
    values[m] <- linear_statistic(design, Y, target, statistic)
  }
  if (any(!is.finite(values))) stop("Nonfinite candidate statistic")
  pvalue <- rank_pvalue(values[1L], values[-1L], p_rule)
  list(pvalue = pvalue, statistic = values[1L],
       diagnostics = list(randomized_statistics = values[-1L], p_rule = p_rule))
}

# Partial F-test: retain both fits, their residuals, and the original anova call.
partial_f_test <- function(X, Y, target) {
  check_test_inputs(X, Y, target)
  nuisance <- X[, setdiff(seq_len(ncol(X)), target), drop = FALSE]
  null_formula <- if (ncol(nuisance)) Y ~ nuisance else Y ~ 1
  fit_null <- stats::glm(null_formula, family = "gaussian")
  fit_full <- stats::glm(Y ~ X, family = "gaussian")
  residual_null <- fit_null$residuals
  residual_full <- fit_full$residuals
  comparison <- stats::anova(fit_null, fit_full, test = "F")
  pvalue <- comparison[2L, 6L]
  if (length(pvalue) != 1L || !is.finite(pvalue) || pvalue < 0 || pvalue > 1)
    stop("Baseline objective is uncertified or p-value invalid")
  list(pvalue = as.numeric(pvalue), statistic = comparison[2L, 5L],
       diagnostics = list(rss_null = sum(residual_null^2),
                          rss_full = sum(residual_full^2),
                          df_numerator = comparison[2L, 3L],
                          df_denominator = fit_full$df.residual))
}
