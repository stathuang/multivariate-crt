# Run bookkeeping for this independent extension only.
load_lw_sources <- function(experiment_dir, envir, entry) {
  files <- paste0("R/", c("run_tools.R", "oracle_common.R", "oracle_size.R", "oracle_power.R",
                         "oracle_reporting.R", "lw.R", "lw_run_tools.R"))
  bytes <- setNames(lapply(file.path(experiment_dir, files), function(path)
    readBin(path, "raw", n = file.info(path)$size)), files)
  for (content in bytes) eval(parse(text = rawToChar(content)), envir = envir)
  bytes[[entry]] <- readBin(file.path(experiment_dir, entry), "raw", n = file.info(file.path(experiment_dir, entry))$size)
  list(bytes = bytes, sha256 = vapply(bytes, digest::digest, "", algo = "sha256", serialize = FALSE))
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

write_lw_size_table <- function(summary, settings, methods, path, trial_run = FALSE) {
  x <- summary[summary$statistic=='ridge gain',]
  lines <- c('\\begin{tabular}{rrrrccccc}', '\\hline',
    '$n$ & $s$ & $t$ & $d$ & Oracle CRT & Plug-in MLE & Plug-in (d.f.) & MVN-CRT & LW (D.10) \\\\', '\\hline')
  for (i in seq_len(nrow(settings))) {
    cells <- x[x$setting==settings$setting[i],]
    cells <- cells[match(methods,cells$method),]
    value <- ifelse(is.finite(cells$empirical_size), sprintf('%.3f (%.3f)',cells$empirical_size,cells$se),'--')
    incomplete <- is.na(cells$complete) | !cells$complete
    value[incomplete] <- paste0(value[incomplete],'$^{\\dagger}$')
    lines <- c(lines,paste0(paste(c(unlist(settings[i,c('n','s','t','d')]),value),collapse=' & '),' \\\\'))
  }
  lines <- c(lines,'\\hline',if(trial_run)'\\multicolumn{9}{l}{\\footnotesize TRIAL: code verification output} \\\\', '\\end{tabular}')
  writeLines(lines,path)
}

plot_lw_power <- function(summary_rows, settings, theta_grid, method_names, alpha, trial_run = FALSE) {
  figure <- plot_oracle_power(summary_rows, settings, theta_grid, method_names, alpha, trial_run)
  # Preserve the four existing ggplot hue colors and shapes; add one new key.
  colors <- setNames(c(scales::hue_pal()(4), "#A6761D"), method_names)
  shapes <- setNames(c(16, 17, 15, 3, 8), method_names)
  figure + ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_shape_manual(values = shapes) +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 2, byrow = TRUE),
                   shape = ggplot2::guide_legend(nrow = 2, byrow = TRUE))
}
