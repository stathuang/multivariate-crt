# Exact cv.glmnet / separate glmnet refit sequence from
# readable-experiments/R/gaussian_tests.R:18-24, with metadata retained.
background_lasso <- function(nuisance, Y, seed) {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
  cv <- glmnet::cv.glmnet(nuisance, Y, family = "gaussian", nfolds = 5L)
  fit <- glmnet::glmnet(nuisance, Y, family = "gaussian")
  fitted <- as.numeric(stats::predict(fit, nuisance, s = cv$lambda.1se,
                                      type = "response"))
  list(residual = as.numeric(Y) - fitted, fitted = fitted,
       lambda_1se = cv$lambda.1se, lambda_min = cv$lambda.min)
}
