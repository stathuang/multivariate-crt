# Summarize rejection indicators; failed and missing attempts stay visible.
report_cell_keys <- c("scenario", "theta", "graph", "method")
report_key <- function(x, cols = report_cell_keys) {
  do.call(paste, c(x[cols], sep = "\034"))
}

summarize_results <- function(results, expected, alpha = 0.05, metric = "rejection") {
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1)
    stop("alpha must be a number strictly between zero and one")
  metric <- match.arg(metric, c("rejection", "FDP", "power"))
  required <- c(report_cell_keys, "replicate", "status", if (metric == "rejection") "pvalue" else metric)
  if (!all(required %in% names(results))) stop("Missing result columns")
  if (!all(c(report_cell_keys, "expected") %in% names(expected)))
    stop("Expected-cell schema is incomplete")
  if (anyNA(results[c(report_cell_keys, "replicate", "status")]) || anyNA(expected))
    stop("Missing identifiers/expected counts")
  if (any(expected$expected < 1 | expected$expected != as.integer(expected$expected)))
    stop("Expected counts must be positive integers")

  result_keys <- report_key(results)
  expected_keys <- report_key(expected)
  if (anyDuplicated(expected_keys)) stop("Duplicate expected cells")
  if (anyDuplicated(report_key(results, c(report_cell_keys, "replicate"))))
    stop("Duplicate replicate/method records")
  if (any(!result_keys %in% expected_keys)) stop("Unexpected result cells")
  if (any(results$replicate < 1 | results$replicate != as.integer(results$replicate)))
    stop("Invalid replicate identifiers")

  rows <- vector("list", nrow(expected))
  for (i in seq_len(nrow(expected))) {
    cell <- results[result_keys == expected_keys[i], , drop = FALSE]
    planned <- expected$expected[i]
    if (any(cell$replicate > planned)) stop("Replicate exceeds expected count")
    raw <- if (metric == "rejection") cell$pvalue else cell[[metric]]
    valid <- cell$status == "ok" & is.finite(raw) & raw >= 0 & raw <= 1
    values <- if (metric == "rejection") as.numeric(raw[valid] < alpha) else raw[valid]
    warnings <- if ("warnings" %in% names(cell)) {
      ifelse(is.na(cell$warnings), "", cell$warnings)
    } else {
      rep("", nrow(cell))
    }

    row <- expected[i, , drop = FALSE]
    row$alpha <- alpha
    row$metric <- metric
    row$attempts <- nrow(cell)
    row$valid <- length(values)
    row$failed <- sum(!valid)
    row$missing <- planned - nrow(cell)
    row$warning_records <- sum(nzchar(warnings))
    row$severe_convergence_warnings <- sum(grepl(
      "convergence not obtained|did not converge|failed to converge|non.?convergen",
      warnings, ignore.case = TRUE
    ))
    row$rejections <- if (metric == "rejection") sum(values) else NA_integer_
    row$estimate <- if (length(values)) mean(values) else NA_real_
    row$se <- if (length(values) > 1L) sd(values) / sqrt(length(values)) else NA_real_
    row$complete <- row$valid == planned && row$attempts == planned
    rows[[i]] <- row
  }
  do.call(rbind, rows)
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
