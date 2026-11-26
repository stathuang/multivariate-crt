# Separate computational success, objective accuracy and inferential reference.
add_result_diagnostics <- function(row, diagnostic = NULL, condition = NULL) {
  if (is.null(diagnostic) && inherits(condition, "condition")) diagnostic <- condition$diagnostics
  defaults <- list(raw_pvalue = NA_real_, error_code = "", objective_status = "not_applicable",
    inference_status = "original_method_reference", numerical_status_policy = "original_method",
    likelihood_statistic = NA_real_, statistic_lower = NA_real_, statistic_upper = NA_real_,
    df = NA_real_, declared_reference_df = NA_real_, declared_reference_pvalue = NA_real_,
    declared_classes = NA_integer_, observed_classes = NA_integer_, reference_policy = "original_method", full_objective_gap = NA_real_, reduced_objective_gap = NA_real_,
    n = NA_integer_, n_retained = NA_integer_, prob_filter = NA_real_,
    selected_lambda = NA_real_, initial_beta_supplied = NA)
  for (name in names(defaults)) if (!name %in% names(row)) row[[name]] <- defaults[[name]]
  if (!is.null(diagnostic)) {
    direct <- intersect(names(defaults), names(diagnostic))
    for (name in direct) row[[name]] <- diagnostic[[name]]
    if (!is.null(diagnostic$pvalue)) row$raw_pvalue <- diagnostic$pvalue
    if (!is.null(diagnostic$statistic)) row$likelihood_statistic <- diagnostic$statistic
    if (!is.null(diagnostic$full$gap)) row$full_objective_gap <- diagnostic$full$gap
    if (!is.null(diagnostic$reduced$gap)) row$reduced_objective_gap <- diagnostic$reduced$gap

  }
  if (inherits(condition, "condition")) {
    row$error_code <- class(condition)[1L]
    if (all(row$method %in% c("Chi-square", "LRT"))) {
      row$objective_status <- "uncertified"
      row$inference_status <- "working_chisquare_reference_calibration_not_proved"

    }
    if (inherits(condition, "CGM_NO_RETAINED_OBSERVATIONS")) {
      row$objective_status <- "undefined_domain"
      row$numerical_status_policy <- "cgm_retained_count_guard_v1"
      for (name in c("n", "n_retained", "prob_filter", "selected_lambda", "initial_beta_supplied"))
        row[[name]] <- condition[[name]]

    }

  }
  attr(row, "diagnostics") <- list(diagnostic)
  row

}
bind_diagnostic_rows <- function(rows) {
  answer <- do.call(rbind, rows)
  attr(answer, "diagnostics") <- unlist(lapply(rows, function(x) attr(x, "diagnostics")), recursive = FALSE)
  answer

}
