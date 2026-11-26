# Gaussian, binary-logistic and numerical-response RF dCRT. No seed is set here.
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
# in ../../../provenance/dcrt-source.md (relative to this file).
dcrt_test <- function(X, y, target, response_model, keep_residuals = FALSE) {
  .baseline_require("glmnet")
  if (response_model == "RF") .baseline_require("randomForest")
  check_test_inputs(X, y, target)
  response_model <- match.arg(response_model, c("gaussian", "binomial", "RF"))
  if (response_model == "binomial" && !all(y %in% c(0, 1)))
    stop("Binomial response must use numeric 0/1 coding")
  p <- ncol(X)
  if (p < 3L) stop("Historical dCRT requires at least two nuisance predictors")
  individual <- numeric(length(target))
  z_statistics <- numeric(length(target))
  # E12's oracle comparison reuses the response fit through this explicit
  # optional return. Ordinary experiments do not retain the residual arrays.
  residuals <- if (keep_residuals) vector("list", length(target)) else NULL
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
      # Retain regression on numerical y, including binary 0/1 responses.
      response_fit <- randomForest::randomForest(x = nuisance, y = y, ntree = 100)
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
    z_statistics[index] <- sqrt(nrow(X)) * importance / sqrt(variance)
    individual[index] <- 2 * stats::pnorm(-z_statistics[index])
    if (keep_residuals) residuals[[index]] <- list(
      target = j, covariate = as.numeric(x_residual), response = as.numeric(y_residual),
      conditional_variance = conditional_variance
    )
  }
  names(individual) <- names(z_statistics) <- paste0("X", target)
  diagnostics <- list(coordinate_z = z_statistics, coordinate_pvalues = individual,
                      combination = "Bonferroni")
  if (keep_residuals) diagnostics$residuals <- residuals
  list(pvalue = .baseline_bonferroni(individual),
       statistic = max(z_statistics),
       diagnostics = diagnostics)
}
