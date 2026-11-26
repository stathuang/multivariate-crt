# E11: retain one Gaussian dataset and compare graphs with additional edges.
# The same permutation is applied to every graph. Creating these graphs uses
# no random numbers, so the following sampling and fitting order is preserved.
supergraph_design <- function(config) {
  design <- gaussian_design(config$n, config$p, config$bandwidth,
                            config$precision_offdiag)
  distance <- abs(outer(seq_len(config$p), seq_len(config$p), `-`))
  graphs <- vector("list", length(config$graph_bandwidths))
  for (i in seq_along(graphs)) {
    graph <- 1L * (distance > 0 & distance <= config$graph_bandwidths[i])
    graphs[[i]] <- graph[design$permutation, design$permutation, drop = FALSE]
  }
  design$graphs <- setNames(graphs, config$graph_ids)
  design
}

# E14/E15: Student-t innovations and a finite autoregression across columns.
# The first column is unscaled t; subsequent innovations use noise_scale.
# Lag 1 always refers to the immediately preceding column, including the
# initial columns where fewer than six predecessors are available.
t3_design <- function(n, p, bandwidth, lag_coefficients, df, noise_scale) {
  if (n < 2L || p < 3L || bandwidth < 0L || bandwidth >= p ||
      !length(lag_coefficients) || any(!is.finite(lag_coefficients)) ||
      !is.finite(df) || df <= 0 || !is.finite(noise_scale) || noise_scale <= 0)
    stop("Invalid Student-t design parameters")
  X <- matrix(0, n, p)
  X[, 1L] <- stats::rt(n, df)
  for (j in 2:p) {
    lag <- seq_len(min(j - 1L, length(lag_coefficients)))
    previous <- X[, j - lag, drop = FALSE]
    X[, j] <- drop(previous %*% lag_coefficients[lag]) +
      noise_scale * stats::rt(n, df)
  }

  # Unlike the Gaussian design, draw all innovations before the permutation.
  permutation <- sample.int(p)
  X <- X[, permutation, drop = FALSE]
  distance <- abs(outer(seq_len(p), seq_len(p), `-`))
  graph <- 1L * (distance > 0 & distance <= bandwidth)
  graph <- graph[permutation, permutation, drop = FALSE]
  beta <- numeric(p)
  active <- seq_len(min(p, 20L))
  beta[active] <- stats::runif(length(active), 1, 2) / sqrt(length(active))
  noise_uniform <- stats::runif(n)
  list(X = X, graph = graph, beta = beta, noise_uniform = noise_uniform,
       permutation = permutation)
}
