# Existing comparator methods. See docs/BASELINES.md for scope and provenance.
# Randomness belongs to the caller; these functions never reset the RNG.
.baseline_require <- function(package) {
  if (!requireNamespace(package, quietly = TRUE))
    stop("Required baseline package is unavailable: ", package, call. = FALSE)
}

.baseline_bonferroni <- function(pvalues) {
  if (!length(pvalues) || any(!is.finite(pvalues)) ||
      any(pvalues < 0 | pvalues > 1)) stop("Invalid individual baseline p-values")
  min(1, length(pvalues) * min(pvalues))
}

# Restricted reconstruction of Liu's d_CRT: CDR='No', k=0, MC_free=TRUE,
# Sigma_X=NULL, Gen_X=NULL. Other upstream branches are deliberately excluded.
# Algorithmic source attribution and unresolved redistribution status are recorded
# in vendor/PROVENANCE.md; this function is a derivative, not a licensed release.
.baseline_dcrt <- function(X, y, target, response_model) {
  .baseline_require("glmnet")
  if (response_model == "RF") .baseline_require("randomForest")
  p <- ncol(X)
  if (p < 3L) stop("Historical dCRT requires at least two nuisance predictors")
  individual <- numeric(length(target))
  for (index in seq_along(target)) {
    j <- target[index]
    nuisance <- X[, -j, drop = FALSE]
    # Keep the historical fit sequence, including a separate glmnet refit.
    conditional_cv <- glmnet::cv.glmnet(nuisance, X[, j], alpha = 1,
                                       intercept = TRUE, dfmax = as.integer(p / 2))
    conditional_fit <- glmnet::glmnet(nuisance, X[, j], alpha = 1,
                                     lambda = conditional_cv$lambda.min,
                                     intercept = TRUE, dfmax = as.integer(p / 2))
    x_residual <- X[, j] - predict(conditional_fit, nuisance)
    conditional_variance <- mean(x_residual^2)
    if (response_model == "RF") {
      response_fit <- stock_randomForest(x = nuisance, y = y, ntree = 1000L)
      y_residual <- y - response_fit$predicted
    } else {
      response_cv <- glmnet::cv.glmnet(nuisance, y, alpha = 1,
                                     family = response_model, lambda = NULL,
                                     dfmax = as.integer(p / 2))
      response_fit <- glmnet::glmnet(nuisance, y, alpha = 1,
                                    lambda = response_cv$lambda.min,
                                    family = response_model,
                                    dfmax = as.integer(p / 2))
      fitted_link <- predict(response_fit, nuisance)
      fitted_response <- if (response_model == "binomial")
        1 / (1 + exp(-fitted_link)) else fitted_link
      y_residual <- y - fitted_response
    }
    importance <- abs(mean(x_residual * y_residual))
    variance <- mean(y_residual^2) * conditional_variance
    if (!is.finite(variance) || variance <= 0)
      stop("dCRT residual variance is nonpositive or nonfinite")
    individual[index] <- 2 * stats::pnorm(-sqrt(nrow(X)) * importance / sqrt(variance))
  }
  .baseline_bonferroni(individual)
}

# Return one group p-value; failures remain errors for the runner to record.
baseline_test_detail <- function(X, y, target, method, response_family = "gaussian") {
  allowed <- c("F-test", "Chi-square", "LRT", "De-Lasso", "CGM", "dCRT (L1)", "dCRT (RF)")
  if (length(method) != 1L || !method %in% allowed) stop("Unknown baseline method")
  if (length(response_family) != 1L ||
      !response_family %in% c("gaussian", "binomial", "multinomial"))
    stop("Unsupported response family")
  if (!is.matrix(X) || !is.numeric(X) || any(!is.finite(X)) ||
      nrow(X) != length(y) || nrow(X) < 2L || anyNA(y)) stop("Invalid baseline data")
  if (!is.numeric(target) || !length(target) || anyNA(target) ||
      any(target != as.integer(target)) || anyDuplicated(target) ||
      any(target < 1L | target > ncol(X))) stop("Invalid target indices")
  if (is.numeric(y) && any(!is.finite(y))) stop("Nonfinite response")
  if (response_family == "gaussian" && !is.numeric(y)) stop("Gaussian response must be numeric")
  # Legacy binomial calls use numeric 0/1; do not silently change coding.
  if (response_family == "binomial" && (!is.numeric(y) || !all(y %in% c(0, 1))))
    stop("Binomial response must use numeric 0/1 coding")
  diagnostic <- list(objective_status="not_applicable", inference_status="original_method_reference")
  nuisance <- X[, setdiff(seq_len(ncol(X)), target), drop = FALSE]
  reduced_formula <- if (ncol(nuisance)) y ~ nuisance else y ~ 1
  if (method %in% c("F-test", "Chi-square")) {
    required_family <- if (method == "F-test") "gaussian" else "binomial"
    if (response_family != required_family) stop(method, " requires ", required_family)
    if (method == "Chi-square") {
      diagnostic <- binary_lrt_supremum(X, y, target)
      diagnostic$objective_status <- "numerically_certified"
      diagnostic$inference_status <- "working_chisquare_reference_calibration_not_proved"
      diagnostic$numerical_status_policy <- "binary_face_exact_dual_v3"
      value <- diagnostic$pvalue
    } else {
      reduced <- stats::glm(reduced_formula, family = response_family)
      full <- stats::glm(y ~ X, family = response_family)
      value <- stats::anova(reduced, full, test = "F")[2, 6]
    }
  } else if (method == "LRT") {
    if (response_family != "multinomial") stop("Historical LRT requires multinomial response")
    diagnostic <- bcr_supremum_lrt(X,y,target)
    # VGAM's historical reference uses observed classes when factor levels are absent.
    # Keep the helper's declared-level tail as a separate diagnostic.
    diagnostic$declared_reference_df <- diagnostic$df
    diagnostic$declared_reference_pvalue <- diagnostic$pvalue
    diagnostic$declared_classes <- nlevels(as.factor(y))
    diagnostic$observed_classes <- nlevels(droplevels(as.factor(y)))
    diagnostic$df <- length(target) * (diagnostic$observed_classes - 1L)
    if (is.finite(diagnostic$pvalue)) diagnostic$pvalue <- stats::pchisq(
      max(0, diagnostic$statistic), diagnostic$df, lower.tail=FALSE)
    diagnostic$reference_policy <- "historical_observed_class_chisquare"
    diagnostic$numerical_status_policy <- "multinomial_face_entropy_dual_v5"
    value <- diagnostic$pvalue
  } else if (method == "De-Lasso") {
    if (response_family == "multinomial") stop("De-Lasso does not support multinomial")
    .baseline_require("hdi")
    fit <- hdi::lasso.proj(X, y, family = response_family, multiplecorr.method = "holm")
    value <- .baseline_bonferroni(fit$pval[target])
  } else if (method == "CGM") {
    if (response_family == "multinomial") stop("CGM does not support multinomial")
    .baseline_require("SIHR")
    fit <- cgm_lf_diagnostic(X, y, diag(ncol(X))[, target, drop = FALSE],
                    model = if (response_family == "gaussian") "linear" else "logistic",
                    intercept = TRUE, intercept.loading = FALSE)
    value <- .baseline_bonferroni(2 * stats::pnorm(abs(fit$est.debias.vec / fit$se.vec), lower.tail = FALSE))
  } else {
    # Historical BCR comparator: scalar numeric class labels, Gaussian Lasso/RF.
    model <- if (method == "dCRT (RF)") "RF" else
      if (response_family == "binomial") "binomial" else "gaussian"
    value <- .baseline_dcrt(X, as.numeric(y), target, model)
  }
  if (length(value) != 1L || !is.finite(value) || value < 0 || value > 1)
    stop(structure(list(message="Baseline objective is uncertified or p-value invalid",call=NULL,
                        diagnostics=diagnostic),class=c("BASELINE_UNCERTIFIED","error","condition")))
  list(pvalue=as.numeric(value), diagnostics=diagnostic)
}

# Scalar compatibility interface. Detailed orchestration calls the function above.
baseline_test <- function(X, y, target, method, response_family="gaussian") {
  baseline_test_detail(X, y, target, method, response_family)$pvalue
}
