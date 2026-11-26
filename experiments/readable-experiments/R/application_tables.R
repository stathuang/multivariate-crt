# Preserve the established table format, selecting cells by row and method names.
format_application_pvalue <- function(p) {
  ifelse(!is.finite(p), "--", ifelse(p < .0005, "$<0.001$", sprintf("%.3f", p)))
}

write_application_table <- function(rows, row_ids, methods, path, row_labels = row_ids,
                                    method_labels = methods, observed = FALSE,
                                    stock = TRUE, trial_run = FALSE) {
  key <- paste(rows$row_id, rows$method, sep = "\034")
  if (anyDuplicated(key)) stop("Table contains duplicate named cells")
  endrow <- paste0(" ", intToUtf8(92), intToUtf8(92))
  body <- character(length(row_ids))
  for (i in seq_along(row_ids)) {
    index <- match(paste(row_ids[i], methods, sep = "\034"), key)
    cell <- rows[index, , drop = FALSE]
    if (observed) {
      value <- format_application_pvalue(cell$pvalue)
    } else {
      value <- ifelse(!is.finite(cell$estimate), "--", paste0(
        sprintf("%.3f", cell$estimate), " (",
        ifelse(is.finite(cell$se), sprintf("%.3f", cell$se), "--"), ")"
      ))
    }
    incomplete <- is.na(cell$complete) | !cell$complete
    value[incomplete] <- paste0(value[incomplete], "$^{\\dagger}$")
    body[i] <- paste0(paste(c(row_labels[i], value), collapse = " & "), endrow)
  }
  title <- if (stock) "$\\mc{T}$" else ""
  columns <- paste0("l", if (stock) "|" else "", paste(rep("c", length(methods)), collapse = ""))
  header <- paste0(paste(c(title, method_labels), collapse = " & "), endrow)
  lines <- c("% New local reproduction; see run.rds and counts.csv for provenance.",
    "% Dagger: incomplete cell. -- means unavailable, never zero.",
    paste0("\\begin{tabular}{", columns, "}"), "\\hline",
    if (!stock) paste0(" & \\multicolumn{", length(methods), "}{c}{Method}", endrow),
    header, "\\hline", body, "\\hline",
    if (trial_run) paste0("\\multicolumn{", length(methods) + 1L,
      "}{l}{\\footnotesize TRIAL: code verification output}", endrow),
    "\\end{tabular}")
  writeLines(lines, path)
}
