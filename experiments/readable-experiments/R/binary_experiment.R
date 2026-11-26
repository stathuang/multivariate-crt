# Keep the accepted response/statistic/comparator order without reseeding here.
evaluate_binary <- function(design, sampling, config, replicate, keep_example = FALSE) {
  rows <- list()
  example <- NULL
  for (theta in config$theta) {
    if (design$status == "ok") {
      X <- design$value$X
      response <- capture_step(switch(
        config$setting,
        l_gl_w =, h_gl_w = logistic_response(X, config$target, theta,
                                            design$value$beta, design$value$noise_uniform),
        l_gl_m =, h_gl_m = nonlinear_binary_response(X, theta, design$value$noise_uniform),
        stop("Unsupported binary setting")
      ))
    } else response <- design

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
          "GLM-Dev" =, "GLM-L1-D" =, "GLM-L1-R-SST" =
            gcrt_logistic_test(X, Y, config$target, sampling$value, statistic, config$p_rule),
          "RF" =, "RF-D" =, "RF-RR" =
            gcrt_rf_test(X, Y, config$target, sampling$value, statistic,
                         config$p_rule, response_family = "binomial"),
          stop("Unsupported binary statistic")
        ))
        test$warnings <- join_warnings(sampling$warnings, test$warnings)
        test$elapsed_seconds <- test$elapsed_seconds + sampling$elapsed_seconds
      }
      tests[[statistic]] <- test
    }

    for (method in config$baselines) {
      if (response$status == "failed") {
        test <- response
      } else {
        Y <- response$value
        test <- capture_step(switch(
          method,
          "Chi-square" = binary_lrt_test(X, Y, config$target),
          "CGM" = cgm_logistic_test(X, Y, config$target),
          "dCRT (L1)" = dcrt_test(X, Y, config$target, "binomial"),
          "dCRT (RF)" = dcrt_test(X, Y, config$target, "RF"),
          stop("Unsupported binary comparator")
        ))
      }
      tests[[method]] <- test
    }

    for (method in names(tests)) {
      test <- tests[[method]]
      test$warnings <- join_warnings(design$warnings, response$warnings, test$warnings)
      row <- test_row(test, config, replicate, theta, method)
      # Store likelihood certificates and typed failures alongside scalar results.
      diagnostic <- if (method == "Chi-square") test$value$diagnostics else
        if (method %in% config$baselines)
          list(objective_status = "not_applicable", inference_status = "original_method_reference") else NULL
      row <- add_result_diagnostics(row, diagnostic, test$condition)
      rows[[length(rows) + 1L]] <- row
    }
    if (keep_example && is.null(example))
      example <- list(theta = theta, Y = response$value, tests = tests)
  }
  list(results = bind_diagnostic_rows(rows), example = example)
}
