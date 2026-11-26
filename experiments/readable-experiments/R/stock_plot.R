# E22: retain the historical four curves, order, widths and absence of error bars.
plot_stock_noise <- function(summary, alpha = .05, trial_run = FALSE) {
  methods <- c("LM-L1-R-SSR", "RF-RR", "dCRT (L1)", "dCRT (RF)")
  plotted <- summary[summary$method %in% methods, , drop = FALSE]
  plotted$Test <- factor(plotted$method, levels = methods)
  caption <- if (trial_run) paste0("TRIAL: ", paste(unique(summary$expected), collapse = ", "),
    " repetitions per noise level; code preview only.") else ""
  if (any(!summary$complete)) caption <- paste(caption, "INCOMPLETE: see counts.csv.", sep = "\n")
  level_label <- paste0("Level = ", sub("^0", "", format(alpha, trim = TRUE)))
  ggplot2::ggplot(plotted, ggplot2::aes(x = theta, y = estimate, color = Test,
    shape = Test, linetype = Test, linewidth = Test)) +
    ggplot2::geom_line(na.rm = TRUE) + ggplot2::geom_point(size = 3, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = scales::hue_pal()(4), breaks = methods, drop = FALSE) +
    ggplot2::scale_shape_manual(values = 1:4, breaks = methods, drop = FALSE) +
    ggplot2::scale_linetype_manual(values = c("solid", "dashed", "dotted", "twodash"), breaks = methods, drop = FALSE) +
    ggplot2::scale_linewidth_manual(values = c(.6, .6, 1.2, .6), breaks = methods, drop = FALSE) +
    ggplot2::theme_bw() + ggplot2::theme(text = ggplot2::element_text(size = 18),
      axis.text = ggplot2::element_text(size = 18), legend.position = "right",
      legend.key.height = grid::unit(1.2, "cm"), legend.key.width = grid::unit(1.5, "cm"),
      plot.caption = ggplot2::element_text(size = 10)) +
    ggplot2::labs(x = expression(sigma), y = "Rejection rate", color = "Test", shape = "Test",
      linetype = "Test", linewidth = "Test", caption = caption) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::scale_y_continuous(breaks = seq(0, 1, .25), labels = sprintf("%.2f", seq(0, 1, .25))) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = alpha, size = level_label), linetype = "dotted") +
    ggplot2::scale_size_manual(values = setNames(.5, level_label)) +
    ggplot2::guides(size = ggplot2::guide_legend(title = "\n "))
}
