# E01/E02 retain their original p <= alpha rule and binomial plug-in MC SE.
# These differ from the sample-SD convention of E03–E25; do not silently unify them.
summarize_oracle <- function(rows, repetitions, M, alpha, settings, methods, power = FALSE) {
  keys <- if (power) c('setting','panel','n','s','t','p','d','theta','method') else
    c('setting','n','s','t','p','d','statistic','method')
  cells <- unique(rows[keys])
  result <- vector('list',nrow(cells))
  for (i in seq_len(nrow(cells))) {
    selected <- rep(TRUE,nrow(rows))
    for (name in keys) selected <- selected & rows[[name]] == cells[[name]][i]
    x <- rows[selected,,drop=FALSE]
    if (anyDuplicated(x$replication)) stop('Duplicate replicate in a benchmark cell')
    valid <- x$status=='ok' & is.finite(x$p_value)
    count <- sum(valid)
    estimate <- if (count) mean(x$p_value[valid] <= alpha) else NA_real_
    row <- cells[i,,drop=FALSE]
    row[[if(power)'rejection_rate' else 'empirical_size']] <- estimate
    row$se <- if(count) sqrt(estimate*(1-estimate)/count) else NA_real_
    row$replications <- repetitions
    row$randomizations <- M
    row$randomization_mode <- 'explicit full randomized covariate matrices'
    row$valid <- count
    row$failed <- sum(!valid)
    row$missing <- repetitions-nrow(x)
    row$complete <- count==repetitions && nrow(x)==repetitions
    result[[i]] <- row
  }
  out <- do.call(rbind,result)
  within <- if(power)out$theta else match(out$statistic,c('raw score','poly2 score','ridge gain'))
  out <- out[order(match(out$setting,settings$setting),within,match(out$method,methods)),]
  rownames(out) <- NULL
  out
}

write_oracle_size_table <- function(summary, settings, methods, path, trial_run = FALSE) {
  x <- summary[summary$statistic=='ridge gain',]
  lines <- c('\\begin{tabular}{rrrrcccc}', '\\hline',
    '$n$ & $s$ & $t$ & $d$ & Oracle CRT & Plug-in MLE & Plug-in (d.f.) & MVN-CRT \\\\', '\\hline')
  for (i in seq_len(nrow(settings))) {
    cells <- x[x$setting==settings$setting[i],]
    cells <- cells[match(methods,cells$method),]
    value <- ifelse(is.finite(cells$empirical_size), sprintf('%.3f (%.3f)',cells$empirical_size,cells$se),'--')
    incomplete <- is.na(cells$complete) | !cells$complete
    value[incomplete] <- paste0(value[incomplete],'$^{\\dagger}$')
    lines <- c(lines,paste0(paste(c(unlist(settings[i,c('n','s','t','d')]),value),collapse=' & '),' \\\\'))
  }
  lines <- c(lines,'\\hline',if(trial_run)'\\multicolumn{8}{l}{\\footnotesize TRIAL: code verification output} \\\\', '\\end{tabular}')
  writeLines(lines,path)
}

plot_oracle_power <- function(summary_rows, settings, theta_grid, method_names, alpha, trial_run = FALSE) {
  plot_data <- summary_rows
  plot_data$method <- factor(plot_data$method, levels = method_names)
  plot_data$panel <- factor(plot_data$panel, levels = settings$panel)
  
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = theta,
      y = rejection_rate,
      color = method,
      shape = method
    )
  ) +
    ggplot2::geom_hline(
      yintercept = alpha,
      linetype = "dashed",
      color = "grey45",
      linewidth = 0.3
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = pmax(0, rejection_rate - se),
        ymax = pmin(1, rejection_rate + se)
      ),
      width = 0.025,
      linewidth = 0.35
    ) +
    ggplot2::geom_line(linewidth = 0.55) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::facet_wrap(~panel, nrow = 1) +
    ggplot2::scale_x_continuous(breaks = theta_grid) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.25)
    ) +
    ggplot2::labs(
      x = expression(theta),
      y = "Rejection probability",
      color = NULL,
      shape = NULL
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(nrow = 1, byrow = TRUE),
      shape = ggplot2::guide_legend(nrow = 1, byrow = TRUE)
    ) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 8),
      legend.key.width = grid::unit(0.55, "lines"),
      legend.spacing.x = grid::unit(0.15, "lines"),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey92", color = "grey55")
    )
  if (trial_run) p <- p + ggplot2::labs(caption = "TRIAL: code verification output")
  if (any(!summary_rows$complete)) p <- p + ggplot2::labs(subtitle = "INCOMPLETE: see counts.csv")
  p
}
