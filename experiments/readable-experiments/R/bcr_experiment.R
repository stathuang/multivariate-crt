# One dataset: RF first, then LRT and the two dCRT comparisons in the original order.
evaluate_bcr <- function(generated, sampling, bcr, config, replicate,
                         keep_example = FALSE) {
  tests <- list()
  for (method in config$methods) {
    if (generated$status == "failed") test <- generated
    else if (method %in% config$statistics && sampling$status == "failed") test <- sampling
    else {
      X <- generated$value$X
      Y <- generated$value$Y
      test <- capture_step(switch(method,
        "RF" = gcrt_rf_test(X, Y, bcr$target, sampling$value, "RF", config$p_rule,
                            response_family = "multinomial"),
        "LRT" = bcr_lrt_test(X, Y, bcr$target),
        # Keep declared factor codes for these regression comparisons. Only the
        # classification statistic removes unused levels inside each RF fit.
        "dCRT (L1)" = dcrt_test(X, as.numeric(Y), bcr$target, "gaussian"),
        "dCRT (RF)" = dcrt_test(X, as.numeric(Y), bcr$target, "RF"),
        stop("Unknown breast-cancer method")
      ))
      # Retain the accepted application's narrow nonconvergence warning rule.
      if (method == "LRT" && test$status == "ok" &&
          identical(config$numerical_status_policy, "bcr_lrt_nonconvergence_is_error_v1") &&
          grepl("convergence not obtained|failed to converge|did not converge|iteration limit.*reached|maximum number of iterations.*reached",
                test$warnings, ignore.case = TRUE)) {
        test$status <- "failed"
        test$error <- "BCR LRT fit did not converge; numeric return retained only as raw_pvalue"
      }
      if (method %in% config$statistics) test$warnings <- join_warnings(sampling$warnings, test$warnings)
    }
    test$warnings <- join_warnings(generated$warnings, test$warnings)
    tests[[method]] <- test
  }
  rows <- list()
  for (method in names(tests)) {
    test <- tests[[method]]
    row <- test_row(test, config, replicate, config$theta, method, graph = "fitted")
    row$status[row$status == "failed"] <- "error"
    row$shared_elapsed_seconds <- generated$elapsed_seconds + sampling$elapsed_seconds
    row$dataset <- "bcr"
    row$design <- config$design
    row$panel <- bcr$panel
    row$raw_pvalue <- if (method == "LRT" && !is.null(test$value)) test$value$pvalue else NA_real_
    row$numerical_status_policy <- config$numerical_status_policy
    diagnostic <- if (method == "LRT") test$value$diagnostics else NULL
    row <- add_result_diagnostics(row, diagnostic, test$condition)
    rows[[method]] <- row
  }
  list(results = bind_diagnostic_rows(rows),
       example = if (keep_example) list(X = generated$value$X, Y = generated$value$Y, tests = tests) else NULL)
}
