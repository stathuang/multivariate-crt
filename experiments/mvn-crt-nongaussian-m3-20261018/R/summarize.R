paths <- list.files("results/formal", pattern = "^cell-", full.names = TRUE)
stopifnot(length(paths) == 12L)
all_rows <- list()
summary_rows <- list()
for (i in seq_along(paths)) {
  result <- readRDS(file.path(paths[i], "results.rds"))
  tab <- result$replications
  all_rows[[i]] <- tab
  ok <- tab$status == "ok"
  N <- sum(ok)
  k <- sum(tab$reject[ok])
  ci <- wilson_interval(k, N)
  summary_rows[[i]] <- data.frame(cell_id = tab$cell_id[1], distribution = tab$distribution[1],
    n = tab$n[1], attempted = nrow(tab), completed = N, failed = sum(!ok), rejections = k,
    rejection_rate = k / N, mc_se = sqrt((k / N) * (1 - k / N) / N),
    wilson_lower = ci[1], wilson_upper = ci[2],
    M = config$M, nominal_level = config$alpha, elapsed_seconds = result$elapsed_seconds)
}
summary <- do.call(rbind, summary_rows)
all_rows <- do.call(rbind, all_rows)
write.csv(summary, "results/summary.csv", row.names = FALSE)
write.csv(all_rows, "results/replications.csv", row.names = FALSE)
saveRDS(summary, "results/summary.rds")
table_lines <- c("# Rejection rates under the null", "",
 "Each cell uses 1,000 independent datasets, M = 199, and nominal level 0.05.",
 "Intervals are pointwise 95% Wilson binomial intervals; MC SE is estimated from the observed rejection rate.", "",
 "| Error distribution | n | Rejections / completed | Failed | Rejection rate | MC SE | 95% Wilson interval |",
 "|---|---:|---:|---:|---:|---:|---:|")
for (i in seq_len(nrow(summary))) {
  r <- summary[i, ]
  table_lines <- c(table_lines, sprintf("| %s | %d | %d / %d | %d | %.4f | %.4f | [%.4f, %.4f] |",
      r$distribution, r$n, r$rejections, r$completed, r$failed, r$rejection_rate,
      r$mc_se, r$wilson_lower, r$wilson_upper))
}
writeLines(table_lines, "results/SUMMARY.md")

draw_figure <- function() {
  par(mar = c(5.0, 4.5, 3.2, 1.2), mgp = c(2.6, 0.8, 0), las = 1,
      family = "sans", cex.axis = 0.95, cex.lab = 1.05)
  plot(NA, xlim = c(42, 480), ylim = c(0, max(0.1, summary$wilson_upper + 0.005)),
       log = "x", xaxt = "n", xlab = "Sample size n", ylab = "Null rejection rate", bty = "l")
  axis(1, at = config$sample_sizes, labels = config$sample_sizes)
  abline(h = seq(0, 0.15, 0.025), col = "#E8E8E8", lwd = 0.8)
  abline(h = config$alpha, lty = 2, col = "#666666", lwd = 1.4)
  colors <- c("#333333", "#0072B2", "#D55E00")
  shapes <- c(16, 17, 15)
  offsets <- c(-0.035, 0, 0.035)
  for (j in seq_along(config$distributions)) {
    r <- summary[summary$distribution == config$distributions[j], ]
    x <- r$n * exp(offsets[j])
    arrows(x, r$wilson_lower, x, r$wilson_upper, code = 3, angle = 90,
           length = 0.045, col = colors[j], lwd = 1.25)
    lines(x, r$rejection_rate, type = "b", col = colors[j], pch = shapes[j], lwd = 1.5, cex = 1.0)
  }
  legend("topright", legend = c("Gaussian", "Standardized multivariate t5",
       "Standardized multivariate t3", "Nominal level 0.05"),
       col = c(colors, "#666666"), pch = c(shapes, NA), lty = c(1, 1, 1, 2),
       lwd = c(1.5, 1.5, 1.5, 1.4), bty = "n", cex = 0.85)
  title(main = "Gaussian MVN-CRT under non-Gaussian target errors", cex.main = 1.05)
  mtext("1,000 repetitions per cell; 199 copies; pointwise 95% Wilson intervals", side = 3, line = 0.3, cex = 0.8)
  mtext("Points are slightly offset horizontally for readability.", side = 1, line = 4, cex = 0.7)
}
pdf("figures/rejection-rates.pdf", width = 9, height = 6, useDingbats = FALSE)
draw_figure()
dev.off()
# PDF needs no X11/Cairo. Poppler provides the PNG preview without those R libraries.
renderer <- Sys.which("pdftoppm")
if (nzchar(renderer)) {
  dir.create("tmp/cache", recursive = TRUE, showWarnings = FALSE)
  Sys.setenv(XDG_CACHE_HOME = normalizePath("tmp/cache"))
  status <- system2(renderer, c("-f", "1", "-singlefile", "-r", "200", "-png",
                               "figures/rejection-rates.pdf", "figures/rejection-rates"))
  if (status != 0L) stop("PDF generated, but Poppler PNG rendering failed")
} else {
  message("PDF generated. Install Poppler/pdftoppm to produce the PNG preview.")
}
print(summary, row.names = FALSE)
