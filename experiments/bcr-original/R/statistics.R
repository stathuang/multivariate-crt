# Named G-CRT statistics, preserving the historical fit settings.
.ms_require <- function(package) {
  if (!requireNamespace(package, quietly = TRUE))
    stop("Required package is unavailable: ", package)
}

.ms_glm_stat <- function(X, y, target, family, sst = FALSE) {
  fit <- stats::glm(y ~ X, family = family)
  if (!sst) return(-fit$deviance)
  coefficients <- summary(fit)$coefficients
  if (fit$rank < ncol(X) + 1L || fit$df.residual <= 0L)
    stop("SST requires an identifiable fit with positive residual degrees of freedom")
  sum(coefficients[1L + target, 3L]^2)
}

.ms_rf_stat <- function(X, y, target, classification = FALSE) {
  .ms_require("randomForest")
  if (classification) y <- as.factor(y)
  fit <- randomForest::randomForest(X, y, ntree = 100L, importance = TRUE)
  column <- if (classification) ncol(fit$importance) - 1L else 1L
  sum(pmax(fit$importance[target, column], 0))
}

.ms_lasso_null <- function(X, y, family) {
  # glmnet needs >=2 columns. The intercept-only and one-column boundaries
  # have explicit equivalents; no artificial target information is introduced.
  if (!ncol(X)) return(rep(mean(as.numeric(y)), nrow(X)))
  .ms_require("glmnet")
  if (ncol(X) == 1L) X <- cbind(X, zero_padding = 0)
  cv <- glmnet::cv.glmnet(X, y, family = family, nfolds = 5L)
  fit <- glmnet::glmnet(X, y, family = family)
  as.numeric(stats::predict(fit, X, s = cv$lambda.1se, type = "response"))
}

.ms_rf_null <- function(X, y, classification = FALSE) {
  if (!ncol(X)) {
    if (classification) return(factor(rep(names(which.max(table(y))), length(y)),
                                     levels = levels(as.factor(y))))
    return(rep(mean(as.numeric(y)), nrow(X)))
  }
  .ms_require("randomForest")
  if (classification) y <- as.factor(y)
  randomForest::randomForest(X, y, ntree = 100L, importance = TRUE)$predicted
}

# Returns the M+1 candidate statistics for testing/auditing. The observed
# matrix is first; every candidate is evaluated by the same function.
gcrt_statistics <- function(X, y, target, copies, statistic,
                            response_family = "gaussian") {
  .ms_matrix(X)
  target <- .ms_target(target, ncol(X))
  response_family <- match.arg(response_family, c("gaussian", "binomial", "multinomial"))
  if (length(y) != nrow(X) || anyNA(y) ||
      (!is.factor(y) && (!is.numeric(y) || any(!is.finite(y)))))
    stop("y must be a finite numeric vector or factor of length nrow(X)")
  if (response_family == "gaussian" && is.factor(y)) stop("gaussian y must be numeric")
  labels <- c("LM-SST", "LM-SSR", "LM-L1-R-SST", "LM-L1-R-SSR",
              "RF", "RF-D", "RF-RR", "GLM-Dev", "GLM-L1-D", "GLM-L1-R-SST")
  statistic <- match.arg(statistic, labels)
  if (response_family == "multinomial" && statistic != "RF")
    stop("Historical multinomial G-CRT supports RF only")
  if (response_family != "binomial" && startsWith(statistic, "GLM"))
    stop("GLM statistics require binomial response_family")
  if (response_family != "gaussian" && startsWith(statistic, "LM"))
    stop("LM statistics require gaussian response_family")
  if (response_family == "binomial" &&
      (is.factor(y) || !all(y %in% c(0, 1)))) stop("binomial y must be numeric 0/1")
  full_copies <- FALSE
  if (is.list(copies) && !is.null(copies$X_copies)) {
    full_copies <- identical(copies$reduced, FALSE)
    if (identical(copies$reduced, TRUE) && !identical(as.integer(copies$target), target))
      stop("reduced sampler target order differs from requested target")
    copies <- copies$X_copies
  }
  if (!is.list(copies)) stop("copies must be a list of target matrices or sampler output")
  candidates <- c(list(X[, target, drop = FALSE]), lapply(copies, function(Z) {
    .ms_matrix(Z, "copy")
    if (nrow(Z) != nrow(X)) stop("copy row count differs from X")
    if (full_copies || (ncol(Z) == ncol(X) && ncol(Z) != length(target))) Z <- Z[, target, drop = FALSE]
    if (ncol(Z) != length(target)) stop("copies must have target columns in target order")
    Z
  }))
  nuisance <- X[, setdiff(seq_len(ncol(X)), target), drop = FALSE]
  procedure <- if (statistic %in% c("LM-L1-R-SST", "LM-L1-R-SSR", "GLM-L1-R-SST", "RF-RR"))
    "residual" else if (statistic %in% c("RF-D", "GLM-L1-D")) "distill" else "direct"
  fitted <- NULL
  if (procedure != "direct") {
    if (startsWith(statistic, "RF")) {
      # RF-RR always regresses numerical Y, including binary 0/1 responses.
      fitted <- .ms_rf_null(nuisance, y, statistic == "RF-D" && response_family != "gaussian")
    } else fitted <- .ms_lasso_null(nuisance, y, response_family)
  }
  outcome <- if (procedure == "residual") as.numeric(y) - as.numeric(fitted) else y
  evaluate <- function(Z) {
    if (procedure == "direct") {
      design <- X; design[, target] <- Z; tested <- target
    } else {
      design <- if (procedure == "distill") cbind(Z, fitted) else Z
      tested <- seq_along(target)
    }
    if (startsWith(statistic, "RF"))
      return(.ms_rf_stat(design, outcome, tested,
                        procedure != "residual" && response_family != "gaussian"))
    family <- if (statistic %in% c("GLM-Dev", "GLM-L1-D")) "binomial" else "gaussian"
    .ms_glm_stat(design, outcome, tested, family, grepl("SST$", statistic))
  }
  values <- vapply(candidates, evaluate, numeric(1))
  if (any(!is.finite(values))) stop("nonfinite candidate statistic")
  values
}

gcrt_test <- function(X, y, target, copies, statistic,
                      response_family = "gaussian", p_rule = "discrete") {
  values <- gcrt_statistics(X, y, target, copies, statistic, response_family)
  rank_pvalue(values[1L], values[-1L], p_rule)
}
