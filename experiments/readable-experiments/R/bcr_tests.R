# Full and background-only multinomial likelihoods use the accepted solver.
# Preserve the historical chi-square degrees of freedom based on observed classes.
bcr_lrt_test <- function(X, Y, target) {
  diagnostic <- bcr_supremum_lrt(X, Y, target)
  diagnostic$declared_reference_df <- diagnostic$df
  diagnostic$declared_reference_pvalue <- diagnostic$pvalue
  diagnostic$declared_classes <- nlevels(as.factor(Y))
  diagnostic$observed_classes <- nlevels(droplevels(as.factor(Y)))
  diagnostic$df <- length(target) * (diagnostic$observed_classes - 1L)
  if (is.finite(diagnostic$pvalue)) diagnostic$pvalue <- stats::pchisq(
    max(0, diagnostic$statistic), diagnostic$df, lower.tail = FALSE)
  diagnostic$reference_policy <- "historical_observed_class_chisquare"
  diagnostic$numerical_status_policy <- "multinomial_face_entropy_dual_v5"
  if (length(diagnostic$pvalue) != 1L || !is.finite(diagnostic$pvalue) ||
      diagnostic$pvalue < 0 || diagnostic$pvalue > 1)
    stop(structure(list(message = "Baseline objective is uncertified or p-value invalid",
                        call = NULL, diagnostics = diagnostic),
                   class = c("BASELINE_UNCERTIFIED", "error", "condition")))
  list(pvalue = as.numeric(diagnostic$pvalue), statistic = diagnostic$statistic,
       diagnostics = diagnostic)
}
