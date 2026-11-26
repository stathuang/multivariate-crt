check_binary_inputs <- function(X, Y, target) {
  check_test_inputs(X, Y, target)
  if (!all(Y %in% c(0, 1))) stop("Binomial response must use numeric 0/1 coding")
}

# Fit the background logistic Lasso once per test, retaining CV/refit order.
logistic_background_fit <- function(nuisance, Y) {
  if (!ncol(nuisance)) return(rep(mean(as.numeric(Y)), nrow(nuisance)))
  if (ncol(nuisance) == 1L) nuisance <- cbind(nuisance, zero_padding = 0)
  cv <- glmnet::cv.glmnet(nuisance, Y, family = "binomial", nfolds = 5L)
  fit <- glmnet::glmnet(nuisance, Y, family = "binomial")
  as.numeric(stats::predict(fit, nuisance, s = cv$lambda.1se, type = "response"))
}

# GLM-Dev: logistic regression on all covariates, using negative deviance.
# GLM-L1-D: logistic regression on candidate targets plus background probabilities.
# GLM-L1-R-SST: linear regression of Y minus those probabilities on candidate targets.
gcrt_logistic_test <- function(X, Y, target, copies, statistic, p_rule = "discrete") {
  check_binary_inputs(X, Y, target)
  statistic <- match.arg(statistic, c("GLM-Dev", "GLM-L1-D", "GLM-L1-R-SST"))
  candidates <- target_candidates(X, target, copies)
  fitted_background <- NULL
  residual_response <- NULL
  if (statistic != "GLM-Dev") {
    nuisance <- X[, setdiff(seq_len(ncol(X)), target), drop = FALSE]
    fitted_background <- logistic_background_fit(nuisance, Y)
    residual_response <- as.numeric(Y) - as.numeric(fitted_background)
  }

  values <- numeric(length(candidates))
  for (m in seq_along(candidates)) {
    if (statistic == "GLM-L1-R-SST") {
      values[m] <- linear_statistic(candidates[[m]], residual_response,
                                   seq_along(target), "LM-SST")
    } else {
      if (statistic == "GLM-Dev") {
        design <- X
        design[, target] <- candidates[[m]]
      } else {
        design <- cbind(candidates[[m]], fitted_background)
      }
      fit <- stats::glm(Y ~ design, family = "binomial")
      values[m] <- -fit$deviance
    }
  }
  if (any(!is.finite(values))) stop("nonfinite candidate statistic")
  list(pvalue = rank_pvalue(values[1L], values[-1L], p_rule),
       statistic = values[1L],
       diagnostics = list(randomized_statistics = values[-1L],
                          fitted_background = fitted_background,
                          residual_response = residual_response, p_rule = p_rule))
}

# Keep the verified unpenalized likelihood supremum and its objective certificates.
# The chi-square tail remains the original reference used by the comparator.
binary_lrt_test <- function(X, Y, target) {
  check_binary_inputs(X, Y, target)
  diagnostic <- binary_lrt_supremum(X, Y, target)
  diagnostic$objective_status <- "numerically_certified"
  diagnostic$inference_status <- "working_chisquare_reference_calibration_not_proved"
  diagnostic$numerical_status_policy <- "binary_face_exact_dual_v3"
  if (length(diagnostic$pvalue) != 1L || !is.finite(diagnostic$pvalue) ||
      diagnostic$pvalue < 0 || diagnostic$pvalue > 1)
    stop(structure(list(message = "Baseline objective is uncertified or p-value invalid",
                        call = NULL, diagnostics = diagnostic),
                   class = c("BASELINE_UNCERTIFIED", "error", "condition")))
  list(pvalue = as.numeric(diagnostic$pvalue), statistic = diagnostic$statistic,
       diagnostics = diagnostic)
}

# Preserve the SIHR filter and typed zero-retained error from the shared adapter.
cgm_logistic_test <- function(X, Y, target) {
  check_binary_inputs(X, Y, target)
  fit <- cgm_lf_diagnostic(X, Y, diag(ncol(X))[, target, drop = FALSE],
                          model = "logistic", intercept = TRUE, intercept.loading = FALSE)
  coordinate_z <- abs(fit$est.debias.vec / fit$se.vec)
  individual <- 2 * stats::pnorm(coordinate_z, lower.tail = FALSE)
  list(pvalue = .baseline_bonferroni(individual), statistic = max(coordinate_z),
       diagnostics = list(coordinate_z = coordinate_z, coordinate_pvalues = individual,
                          combination = "Bonferroni"))
}
