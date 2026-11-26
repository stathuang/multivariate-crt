# One R process per Python worker. Binary task files stay in the new directory.
args <- commandArgs(trailingOnly = TRUE)
source(file.path(args[1L], "background.R"))
suppressPackageStartupMessages(library(glmnet))
input <- file("stdin", "r")
repeat {
  line <- readLines(input, n = 1L, warn = FALSE)
  if (!length(line)) break
  warnings <- character()
  result <- tryCatch(withCallingHandlers({
    task <- jsonlite::fromJSON(line)
    con <- file(task$input, "rb")
    values <- readBin(con, "double", n = task$n * (task$p + 1L),
                      size = 8L, endian = "little")
    close(con)
    stopifnot(length(values) == task$n * (task$p + 1L))
    data <- matrix(values, nrow = task$n)
    fit <- background_lasso(data[, seq_len(task$p), drop = FALSE],
                            data[, task$p + 1L], task$seed)
    stopifnot(all(is.finite(fit$residual)))
    con <- file(task$output, "wb")
    writeBin(fit$residual, con, size = 8L, endian = "little")
    close(con)
    list(ok = TRUE, lambda_1se = fit$lambda_1se, lambda_min = fit$lambda_min)
  }, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) list(ok = FALSE, error = conditionMessage(e)))
  result$warnings <- warnings
  cat(jsonlite::toJSON(result, auto_unbox = TRUE, digits = 17), "\n")
  flush(stdout())
}
close(input)
