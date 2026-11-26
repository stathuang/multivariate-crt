# Keep the accepted 8 x 6 inch panel, method order and visual encoding.
# Results/counts retain all methods; the plot shows only those supplied in methods.
plot_curve <- function(summary, methods, trial_run = FALSE) {
  alpha <- unique(summary$alpha)
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1)
    stop("Summary must contain one valid alpha")
  level_label <- paste0("Level = ", sub("^0", "", format(alpha, trim = TRUE)))
  labels <- methods
  labels[labels == "dCRT (L1)"] <- "dCRT"
  n_methods <- length(methods)
  plotted <- summary[summary$method %in% methods, , drop = FALSE]
  plotted$Test <- factor(plotted$method, levels = methods)
  caption <- character()
  if (trial_run) caption <- paste0("TRIAL: ", paste(unique(summary$expected), collapse = ", "),
                                    " repetitions per signal level; code preview only.")
  if (any(!summary$complete)) {
    caption <- c(caption, "INCOMPLETE: see counts.csv for failed or missing results.")
  }
  upper <- min(1, max(c(alpha, plotted$estimate +
    ifelse(is.finite(plotted$se), plotted$se, 0)), na.rm = TRUE) * 1.2)

  ggplot2::ggplot(plotted, ggplot2::aes(
    x = theta, y = estimate, color = Test, shape = Test,
    linetype = Test, linewidth = Test
  )) +
    ggplot2::geom_line(na.rm = TRUE) +
    ggplot2::geom_point(size = 3, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = scales::hue_pal()(n_methods),
                               breaks = methods, labels = labels, drop = FALSE) +
    ggplot2::scale_shape_manual(values = seq_len(n_methods),
                               breaks = methods, labels = labels, drop = FALSE) +
    ggplot2::scale_linetype_manual(values = c("dotdash", "dashed", "solid", "twodash", "longdash")[seq_len(n_methods)],
                                  breaks = methods, labels = labels, drop = FALSE) +
    ggplot2::scale_linewidth_manual(values = rep(0.6, n_methods),
                                   breaks = methods, labels = labels, drop = FALSE) +
    ggplot2::theme_bw() +
    ggplot2::theme(text = ggplot2::element_text(size = 20),
                   axis.text = ggplot2::element_text(size = 20),
                   legend.position = "right",
                   legend.key.height = grid::unit(1.2, "cm"),
                   legend.key.width = grid::unit(1.5, "cm"),
                   plot.caption = ggplot2::element_text(size = 10)) +
    ggplot2::labs(x = expression(theta), y = "Power", color = "Test", shape = "Test",
                  linetype = "Test", linewidth = "Test",
                  caption = paste(caption, collapse = "\n")) +
    ggplot2::coord_cartesian(ylim = c(0, upper)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = estimate - se, ymax = estimate + se),
                          linetype = "solid", width = 0.1, linewidth = 0.4,
                          show.legend = FALSE, na.rm = TRUE) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = alpha, size = level_label),
                         linetype = "dotted") +
    ggplot2::scale_size_manual(values = setNames(0.5, level_label)) +
    ggplot2::guides(size = ggplot2::guide_legend(title = "\n "))
}
