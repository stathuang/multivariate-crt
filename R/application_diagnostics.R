# Stock graph illustration and the three historical graph goodness-of-fit tests.
# Sampling uses R/core.R; all fitting here is deterministic and does not set RNG.
.stock_diagnostic_graph <- function(graph, p) {
  if (!is.matrix(graph) || !identical(dim(graph), c(p,p)) || any(!is.finite(graph)) ||
      any(!graph %in% c(0,1)) || any(diag(graph) != 0) || any(graph != t(graph)))
    stop("Expected a symmetric binary graph with zero diagonal")
  invisible(graph)
}

stock_local_F <- function(X, graph) {
  .ms_matrix(X); .stock_diagnostic_graph(graph, ncol(X))
  n <- nrow(X); p <- ncol(X)
  if (max(rowSums(graph)) + 3 > n) stop("Insufficient residual degrees of freedom for local F tests")
  centered <- scale(X, center = TRUE, scale = FALSE)
  Fvalue <- matrix(0, p, p); pvalue <- matrix(1, p, p)
  # The input is already centered, so these fits do not add an intercept.
  regression_residual <- function(response, predictors) {
    if (!ncol(predictors)) return(response)
    fit <- stats::lm.fit(predictors, response)
    if (fit$rank != ncol(predictors)) stop("Rank-deficient regression design")
    fit$residuals
  }
  for (i in seq_len(p)) {
    neighbors <- setdiff(which(graph[i, ] != 0), i)
    missing <- setdiff(seq_len(p), c(i, neighbors))
    if (!length(missing)) next

    # Null model: regress node i on its graph neighbors.
    residual_null <- regression_residual(
      centered[, i], centered[, neighbors, drop = FALSE]
    )
    df <- n - 2L - length(neighbors)
    for (j in missing) {
      # Larger model: add the nonneighbor j to the same regression.
      residual_larger <- regression_residual(
        centered[, i], centered[, c(neighbors, j), drop = FALSE]
      )
      gain <- sum((residual_null - residual_larger)^2)
      # Compute RSS directly: subtracting two nearly equal sums of squares
      # can lose precision when the larger model fits almost perfectly.
      rss_larger <- sum(residual_larger^2)
      if (rss_larger <= 0) stop("Degenerate augmented residual variance")
      Fvalue[i, j] <- gain / (rss_larger / df)
      pvalue[i, j] <- stats::pf(Fvalue[i, j], 1L, df, lower.tail = FALSE)
    }
  }
  list(Fvalue = Fvalue, pvalue = pvalue)
}
stock_F_sum <- function(X, graph) sum(stock_local_F(X, graph)$Fvalue)
stock_VV <- function(X, graph) {
  .ms_matrix(X); .stock_diagnostic_graph(graph, ncol(X))
  # A complete graph has no missing-edge hypotheses to test. Return before
  # local regressions, whose degrees-of-freedom check is irrelevant here.
  if (all(graph[upper.tri(graph)] == 1)) return(1)
  pv <- stock_local_F(X, graph)$pvalue
  bonf_pvalues <- rep(1, ncol(X))
  for (i in seq_len(ncol(X))) {
    m_i <- sum(graph[i, -i] == 0)
    # A node with no nonneighbors contributes 1, not an artificial zero p-value.
    if (m_i > 0) bonf_pvalues[i] <- m_i * min(pv[i, ])
  }
  # Apply the outer Bonferroni factor across nodes and cap at 1.
  min(1, ncol(X) * min(bonf_pvalues))
}
stock_DP <- function(X, graph) {
  .ms_matrix(X); .stock_diagnostic_graph(graph, ncol(X))
  n <- nrow(X); p <- ncol(X)
  if (2 * max(rowSums(graph)) + 3 > n) stop("Legacy pairwise residual test degree restriction violated")
  centered <- scale(X, center = TRUE, scale = FALSE)
  neighbors <- lapply(seq_len(p), function(i) setdiff(which(graph[i, ] != 0), i))
  pairs <- which(upper.tri(graph) & graph == 0, arr.ind = TRUE)
  if (!nrow(pairs)) return(1)
  pv <- numeric(nrow(pairs))
  for (k in seq_len(nrow(pairs))) {
    i <- pairs[k, 1L]; j <- pairs[k, 2L]; nuisance <- union(neighbors[[i]], neighbors[[j]])
    residual_i <- centered[, i]
    residual_j <- centered[, j]
    if (length(nuisance)) {
      predictors <- centered[, nuisance, drop = FALSE]
      # Regress each endpoint on the union of their graph neighborhoods.
      # Centering has already removed the intercept; fit the two models explicitly.
      fit_i <- stats::lm.fit(predictors, centered[, i])
      fit_j <- stats::lm.fit(predictors, centered[, j])
      if (fit_i$rank != length(nuisance) || fit_j$rank != length(nuisance))
        stop("Rank-deficient pairwise neighborhood")
      residual_i <- fit_i$residuals
      residual_j <- fit_j$residuals
    }
    # With an empty neighborhood, the centered endpoints are already residuals.
    r <- stats::cor(residual_i, residual_j); df <- n - length(nuisance) - 2L
    if (!is.finite(r) || abs(r) > 1 || df <= 0) stop("Invalid partial-correlation test")
    t <- sqrt(df) * r / sqrt(1 - r^2)
    pv[k] <- 2 * stats::pt(abs(t), df, lower.tail = FALSE)
  }
  min(1, min(pv) * length(pv))
}

stock_graph_gof <- function(stock, M = 1000L, L = 3L) {
  # Caller sets seed2023. Preserve historical VV, DP, sampling, F_sum order.
  time <- proc.time()[["elapsed"]]
  vv <- stock_VV(stock$X, stock$graph)
  dp <- stock_DP(stock$X, stock$graph)
  analytical_seconds <- proc.time()[["elapsed"]] - time
  time <- proc.time()[["elapsed"]]
  samples <- mscrt_sample(stock$X, stock$graph, target = seq_len(ncol(stock$X)),
                          M = M, L = L, reduced = FALSE)
  sampling_seconds <- proc.time()[["elapsed"]] - time
  time <- proc.time()[["elapsed"]]
  observed <- stock_F_sum(stock$X, stock$graph)
  randomized <- vapply(samples$X_copies, stock_F_sum, 0.0, graph = stock$graph)
  p <- rank_pvalue(observed, randomized, rule = "discrete")
  statistic_seconds <- proc.time()[["elapsed"]] - time
  pvalues <- c(VV = vv, DP = dp, F_sum = p)
  .stock_diagnostic_pvalues(pvalues)
  list(pvalues = pvalues, observed = observed, randomized = randomized,
       M = M, L = L, analytical_seconds = analytical_seconds, sampling_seconds = sampling_seconds,
       statistic_seconds = statistic_seconds,
       scope = "Fixed fitted graph on saved transformed stock input; graph and transformations selected on these data")
}

plot_stock_graph <- function(stock, path) {
  sectors <- factor(stock$sectors)
  colors <- c("tomato1", "orange", "gold1", "darkolivegreen1", "lightblue", "slateblue2",
              "mediumorchid2", "green", "lightpink1", "royalblue4", "cyan")
  g <- igraph::graph_from_adjacency_matrix(stock$graph, mode = "undirected", diag = FALSE)
  grDevices::pdf(path, width = 8, height = 6)
  on.exit(grDevices::dev.off())
  plot(g, vertex.shape = "circle", vertex.size = 10, vertex.label = stock$stocks,
       vertex.color = colors[as.integer(sectors)], vertex.label.color = "black",
       vertex.label.cex = 0.35, edge.color = "black", edge.width = 0.8)
  graphics::legend("bottomright", legend = levels(sectors), col = colors, pch = 19, pt.cex = 1,
                   cex = 0.5, bty = "o", ncol = 2, lwd = 0)
  invisible(path)
}

.stock_diagnostic_pvalues <- function(pvalues) {
  if (!identical(names(pvalues), c("VV", "DP", "F_sum")) ||
      length(pvalues) != 3L || any(!is.finite(pvalues)) || any(pvalues < 0 | pvalues > 1))
    stop("Expected finite VV, DP and F_sum p-values in [0,1]")
  invisible(pvalues)
}

write_stock_gof_table <- function(result, path) {
  .stock_diagnostic_pvalues(result$pvalues)
  if (result$M != 1000L || result$L != 3L) stop("Historical presentation requires M1000 L3; pilot results stay separate")
  pv <- sprintf("%.3f", result$pvalues[c("VV", "DP", "F_sum")])
  writeLines(c("\\begin{tabular}{l|ccc}", "\\hline",
    " & \\multicolumn{3}{c}{Method} \\\\",
    " & \\nameVV & \\nameDP & \\multicolumn{1}{c}{F$_{\\Sigma}$} \\\\",
    "\\hline", "\\hline",
    paste0("p-value & ", paste(pv, collapse = " & "), " \\\\"),
    "\\hline", "\\end{tabular}"), path)
}
