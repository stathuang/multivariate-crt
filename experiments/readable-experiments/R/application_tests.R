# One observed or generated application dataset, with one shared set of copies.
# All methods remain in their configured order, including unplotted methods.
evaluate_stock <- function(generated, sampling, stock, config, replicate,
                           keep_example = FALSE) {
  target <- stock$targets[[config$group]]
  tests <- list()
  for (method in config$methods) {
    if (generated$status == "failed") test <- generated
    else if (method %in% config$statistics && sampling$status == "failed") test <- sampling
    else {
      X <- generated$value$X
      Y <- generated$value$Y
      test <- capture_step(switch(method,
        "LM-L1-R-SSR" = gcrt_lasso_residual_test(X, Y, target, sampling$value, method, config$p_rule),
        "RF-RR" =, "RF-D" = gcrt_rf_test(X, Y, target, sampling$value, method, config$p_rule),
        "De-Lasso" = de_lasso_test(X, Y, target),
        "dCRT (L1)" = dcrt_test(X, Y, target, "gaussian"),
        "dCRT (RF)" = dcrt_test(X, Y, target, "RF"),
        stop("Unknown stock method")
      ))
      if (method %in% config$statistics) test$warnings <- join_warnings(sampling$warnings, test$warnings)
    }
    test$warnings <- join_warnings(generated$warnings, test$warnings)
    tests[[method]] <- test
  }
  rows <- list()
  for (method in names(tests)) {
    test <- tests[[method]]
    row <- test_row(test, config, replicate, config$theta, method, graph = "fitted")
    # Keep the application result schema's historical error status.
    row$status[row$status == "failed"] <- "error"
    row$shared_elapsed_seconds <- generated$elapsed_seconds + sampling$elapsed_seconds
    row$dataset <- "stock"
    row$design <- config$design
    row$panel <- "stocks"
    rows[[method]] <- row
  }
  list(results = do.call(rbind, unname(rows)),
       example = if (keep_example) list(X = generated$value$X, Y = generated$value$Y, tests = tests) else NULL)
}
