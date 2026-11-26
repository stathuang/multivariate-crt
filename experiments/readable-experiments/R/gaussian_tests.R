# Reduced target matrices, with the observed target first in the same order.
target_candidates <- function(X, target, copies) {
  if (!is.list(copies) || !is.list(copies$X_copies) ||
      !identical(copies$reduced, TRUE) ||
      !identical(as.integer(copies$target), as.integer(target)))
    stop("Expected reduced sampler output in target order")
  candidates <- c(list(X[, target, drop = FALSE]), copies$X_copies)
  for (candidate in candidates) {
    .ms_matrix(candidate, "copy")
    if (!identical(dim(candidate), c(nrow(X), length(target))))
      stop("Copy dimensions differ from the target matrix")
  }
  candidates
}

# The nuisance model uses only background covariates. Keep the CV and refit
# calls separate, in their accepted order, and use the one-standard-error rule.
lasso_background_fit <- function(nuisance, Y) {
  if (!ncol(nuisance)) return(rep(mean(as.numeric(Y)), nrow(nuisance)))
  if (ncol(nuisance) == 1L) nuisance <- cbind(nuisance, zero_padding = 0)
  cv <- glmnet::cv.glmnet(nuisance, Y, family = "gaussian", nfolds = 5L)
  fit <- glmnet::glmnet(nuisance, Y, family = "gaussian")
  as.numeric(stats::predict(fit, nuisance, s = cv$lambda.1se, type = "response"))
}

gcrt_lasso_residual_test <- function(X, Y, target, copies, statistic,
                                    p_rule = "discrete") {
  check_test_inputs(X, Y, target)
  statistic <- match.arg(statistic, c("LM-L1-R-SST", "LM-L1-R-SSR"))
  candidates <- target_candidates(X, target, copies)
  nuisance <- X[, setdiff(seq_len(ncol(X)), target), drop = FALSE]
  fitted_background <- lasso_background_fit(nuisance, Y)
  residual_response <- as.numeric(Y) - as.numeric(fitted_background)

  # Fit the same residual response to each target candidate.
  linear_method <- if (statistic == "LM-L1-R-SST") "LM-SST" else "LM-SSR"
  values <- numeric(length(candidates))
  for (m in seq_along(candidates)) {
    values[m] <- linear_statistic(candidates[[m]], residual_response,
                                 seq_along(target), linear_method)
  }
  if (any(!is.finite(values))) stop("nonfinite candidate statistic")
  list(pvalue = rank_pvalue(values[1L], values[-1L], p_rule),
       statistic = values[1L],
       diagnostics = list(randomized_statistics = values[-1L],
                          fitted_background = fitted_background,
                          residual_response = residual_response, p_rule = p_rule))
}

# RF: predict Y from the candidate target columns and all conditioning covariates.
# RF-D: predict Y from the candidate target columns and the background predictions.
# RF-RR: predict Y minus the background predictions from the candidate target columns.
gcrt_rf_test <- function(X, Y, target, copies, statistic, p_rule = "discrete",
                         response_family = "gaussian") {
  if (identical(response_family, "multinomial")) {
    .ms_matrix(X)
    .ms_target(target, ncol(X))
    if (length(Y) != nrow(X) || anyNA(Y)) stop("Invalid categorical response")
  } else check_test_inputs(X, Y, target)
  response_family <- match.arg(response_family, c("gaussian", "binomial", "multinomial"))
  if (response_family == "binomial" && !all(Y %in% c(0, 1)))
    stop("Binomial response must use numeric 0/1 coding")
  statistic <- match.arg(statistic, c("RF", "RF-D", "RF-RR"))
  candidates <- target_candidates(X, target, copies)
  nuisance <- X[, setdiff(seq_len(ncol(X)), target), drop = FALSE]
  fitted_background <- NULL
  if (statistic != "RF") {
    # Binary RF-D uses predicted classes; RF-RR regresses numerical 0/1 responses.
    background_response <- if (statistic == "RF-D" && response_family %in% c("binomial", "multinomial"))
      as.factor(Y) else Y
    if (ncol(nuisance)) {
      fit_background <- randomForest::randomForest(
        nuisance, background_response, ntree = 100L, importance = TRUE
      )
      fitted_background <- fit_background$predicted  # Out-of-bag predictions, shared across candidates.
    } else {
      fitted_background <- if (is.factor(background_response))
        factor(rep(names(which.max(table(Y))), length(Y)), levels = levels(as.factor(Y))) else
        rep(mean(as.numeric(Y)), nrow(X))
    }
  }
  outcome <- if (statistic == "RF-RR")
    as.numeric(Y) - as.numeric(fitted_background) else Y

  values <- numeric(length(candidates))
  for (m in seq_along(candidates)) {
    if (statistic == "RF") {
      design <- X
      design[, target] <- candidates[[m]]
      tested <- target
    } else if (statistic == "RF-D") {
      # cbind keeps the accepted factor-to-code conversion for binary predictions.
      design <- cbind(candidates[[m]], fitted_background)
      tested <- seq_along(target)
    } else {
      design <- candidates[[m]]
      tested <- seq_along(target)
    }
    classification <- response_family %in% c("binomial", "multinomial") && statistic != "RF-RR"
    forest_response <- outcome
    if (classification) {
      forest_response <- droplevels(as.factor(outcome))
      if (nlevels(forest_response) < 2L)
        stop("RF classification requires at least two observed response classes")
    }
    fit <- randomForest::randomForest(design, forest_response, ntree = 100L, importance = TRUE)
    importance_column <- if (classification) ncol(fit$importance) - 1L else 1L
    # Each statistic sums the positive parts of the target variables' permutation importances.
    values[m] <- sum(pmax(fit$importance[tested, importance_column], 0))
  }
  if (any(!is.finite(values))) stop("nonfinite candidate statistic")
  list(pvalue = rank_pvalue(values[1L], values[-1L], p_rule),
       statistic = values[1L],
       diagnostics = list(randomized_statistics = values[-1L],
                          fitted_background = fitted_background, p_rule = p_rule))
}
