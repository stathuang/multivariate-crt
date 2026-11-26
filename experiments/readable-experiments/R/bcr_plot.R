# E24 uses exactly the LRT rows computed for E23; no second simulation is run.
plot_bcr_histogram <- function(results, summary, panel, trial_run = FALSE) {
  x <- results[results$method == "LRT" & results$status == "ok" & is.finite(results$pvalue), , drop = FALSE]
  caption <- paste(c(paste("Panel:", panel), if (trial_run) "TRIAL: code verification output",
                     if (any(!summary$complete)) "INCOMPLETE: see counts.csv."), collapse = "\n")
  ggplot2::ggplot(x, ggplot2::aes(x = pvalue)) +
    ggplot2::geom_histogram(bins = 20, fill = "grey", color = "white") +
    ggplot2::theme_minimal() +
    ggplot2::theme(text = ggplot2::element_text(size = 20),
      axis.title = ggplot2::element_text(size = 18),
      plot.title = ggplot2::element_text(size = 20, hjust = .5),
      plot.caption = ggplot2::element_text(size = 10)) +
    ggplot2::labs(title = "Histogram", x = "Simulated p-value", y = "Frequency", caption = caption)
}
