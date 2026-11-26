# Evaluate one design in signal/method order. The entry script makes data,
# optional De-Lasso cache, and copies; this function never reseeds the generator.
evaluate_gaussian <- function(design, sampling, cache, config, replicate,
                              keep_example = FALSE,
                              graph_name = paste0("K=", config$bandwidth)) {
  rows <- list()
  example <- NULL
  for (theta_index in seq_along(config$theta)) {
    theta <- config$theta[theta_index]
    if (design$status == "ok") {
      X <- design$value$X
      response <- capture_step(switch(
        config$setting,
        l_l_w =, h_l_w = linear_response(X, config$target, theta,
                                        design$value$beta, design$value$noise_uniform),
        l_l_m =, h_l_m = nonlinear_response(X, theta, design$value$noise_uniform),
        stop("Unsupported Gaussian setting")
      ))
    } else {
      response <- design
    }

    tests <- list()
    for (statistic in config$statistics) {
      if (response$status == "failed") {
        test <- response
      } else if (sampling$status == "failed") {
        test <- sampling
      } else {
        Y <- response$value
        test <- capture_step(switch(
          statistic,
          "LM-SST" =, "LM-SSR" =
            gcrt_linear_test(X, Y, config$target, sampling$value, statistic, config$p_rule),
          "LM-L1-R-SST" =, "LM-L1-R-SSR" =
            gcrt_lasso_residual_test(X, Y, config$target, sampling$value, statistic, config$p_rule),
          "RF" =, "RF-D" =, "RF-RR" =
            gcrt_rf_test(X, Y, config$target, sampling$value, statistic, config$p_rule),
          stop("Unsupported Gaussian statistic")
        ))
        test$warnings <- join_warnings(sampling$warnings, test$warnings)
        test$elapsed_seconds <- test$elapsed_seconds + sampling$elapsed_seconds
      }
      tests[[statistic]] <- test
    }

    for (method in config$baselines) {
      if (response$status == "failed") {
        test <- response
      } else if (method == "De-Lasso" && cache$enabled) {
        test <- cache$result
        if (test$status == "ok") {
          seed <- de_lasso_seed(config, replicate, "response", theta_index)
          test <- capture_step(with_seed(seed, de_lasso_test(
            X, response$value, config$target, cache$result$value,
            config$de_lasso_suppress_unused_group_testing
          )))
          test$warnings <- join_warnings(cache$result$warnings, test$warnings)
          test$elapsed_seconds <- test$elapsed_seconds + cache$result$elapsed_seconds
        }
      } else {
        Y <- response$value
        test <- capture_step(switch(
          method,
          "F-test" = partial_f_test(X, Y, config$target),
          "De-Lasso" = de_lasso_test(X, Y, config$target),
          "CGM" = cgm_linear_test(X, Y, config$target),
          "dCRT (L1)" = dcrt_test(X, Y, config$target, "gaussian"),
          "dCRT (RF)" = dcrt_test(X, Y, config$target, "RF"),
          stop("Unsupported Gaussian comparator")
        ))
      }
      tests[[method]] <- test
    }

    for (method in names(tests)) {
      test <- tests[[method]]
      test$warnings <- join_warnings(design$warnings, response$warnings, test$warnings)
      rows[[length(rows) + 1L]] <- test_row(test, config, replicate, theta, method,
                                          graph = graph_name)
    }
    if (keep_example && is.null(example))
      example <- list(theta = theta, Y = response$value, tests = tests)
  }
  results <- do.call(rbind, rows)
  if (cache$enabled) {
    is_de_lasso <- results$method == "De-Lasso"
    results$de_lasso_policy <- config$de_lasso_policy
    results$de_lasso_seed_policy <- config$de_lasso_seed_policy
    results$de_lasso_nodewise_seed <- ifelse(is_de_lasso, cache$nodewise_seed, NA_integer_)
    results$de_lasso_response_seed <- vapply(seq_len(nrow(results)), function(i) {
      if (is_de_lasso[i]) de_lasso_seed(config, replicate, "response",
                                       match(results$theta[i], config$theta)) else NA_integer_
    }, integer(1))
  }
  list(results = results, example = example)
}
