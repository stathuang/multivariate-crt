# Capture warnings and failures without substituting successful-looking values.
capture_step <- function(expr) {
  warnings <- character()
  error <- ""
  condition <- NULL
  start <- proc.time()[["elapsed"]]
  value <- tryCatch(
    withCallingHandlers(force(expr), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }),
    error = function(e) {
      error <<- conditionMessage(e)
      condition <<- e
      NULL
    }
  )
  list(
    value = value, condition = condition, status = if (nzchar(error)) "failed" else "ok",
    error = error, warnings = paste(unique(warnings), collapse = " | "),
    elapsed_seconds = unname(proc.time()[["elapsed"]] - start)
  )
}

join_warnings <- function(...) {
  messages <- c(...)
  paste(messages[nzchar(messages)], collapse = " | ")
}

test_row <- function(test, config, replicate, theta, method,
                     graph = paste0("K=", config$bandwidth)) {
  answer <- test$value
  if (test$status == "ok" &&
    (length(answer$pvalue) != 1L || !is.finite(answer$pvalue) ||
      answer$pvalue < 0 || answer$pvalue > 1)) {
    test$status <- "failed"
    test$error <- "Method returned an invalid p-value"
  }
  data.frame(
    scenario = config$id, replicate = as.integer(replicate), theta = theta,
    graph = graph, method = method,
    pvalue = if (test$status == "ok") answer$pvalue else NA_real_,
    statistic = if (test$status == "ok") answer$statistic else NA_real_,
    status = test$status, error = test$error, warnings = test$warnings,
    elapsed_seconds = test$elapsed_seconds, stringsAsFactors = FALSE
  )
}

# Derive labels from the parameters actually used by the calculation.
prepare_gaussian_config <- function(config) {
  config$methods <- c(config$statistics, config$baselines)
  if (config$family != "supergraph") config$graph_bandwidths <- config$bandwidth
  if (!length(config$graph_bandwidths) || anyNA(config$graph_bandwidths) ||
    any(config$graph_bandwidths != floor(config$graph_bandwidths)) ||
    any(config$graph_bandwidths < config$bandwidth | config$graph_bandwidths >= config$p) ||
    anyDuplicated(config$graph_bandwidths))
    stop("Graph bandwidths must be distinct integers between the true bandwidth and p - 1")
  config$graph_ids <- paste0("K=", config$graph_bandwidths)
  if (length(config$repetitions) != 1L || !is.finite(config$repetitions) ||
    config$repetitions < 1 || config$repetitions != floor(config$repetitions))
    stop("repetitions must be a positive integer")
  seeds <- config$seed_base + c(1, config$repetitions)
  if (length(seeds) != 2L || any(!is.finite(seeds)) ||
    any(seeds < 0 | seeds > .Machine$integer.max | seeds != floor(seeds)))
    stop("Replicate seeds must be nonnegative R integers")
  config
}


write_report_csv <- function(x, path) {
  # Export 17 significant digits. Use RDS for exact binary replay;
  # decimal parsing in R can differ in the final floating-point bit.
  y <- x
  for (name in names(y)) {
    if (is.double(y[[name]])) {
      y[[name]] <- vapply(y[[name]], function(value) {
        if (is.na(value)) NA_character_ else sprintf("%.17g", value)
      }, character(1))
    }
  }
  write.csv(y, path, row.names = FALSE, na = "")
}

mark_lw_failures <- function(rows, value) {
  if (value$status != "ok") return(rows)
  diagnostic <- value$value$diagnostics
  selected <- rows$method == "LW plug-in CRT (D.10)"
  rows$status[selected] <- diagnostic$lw_status
  rows$error[selected] <- diagnostic$lw_error
  rows$warnings[selected] <- join_warnings(value$warnings, diagnostic$lw_warnings)
  rows
}
