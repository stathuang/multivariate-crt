# E11: finish the signal grid for one graph before sampling the next graph.
# All graphs share X, beta and response noise; each has its own hub and copies.
evaluate_supergraphs <- function(design, config, replicate, keep_example = FALSE) {
  rows <- list()
  examples <- list()
  for (graph_name in config$graph_ids) {
    if (design$status == "ok") {
      graph <- design$value$graphs[[graph_name]]
      sampling <- capture_step(mscrt_sample(
        design$value$X, graph, config$target, M = config$M, L = config$L,
        reduced = TRUE
      ))
    } else {
      sampling <- design
    }
    result <- evaluate_gaussian(design, sampling, list(enabled = FALSE),
                                 config, replicate, keep_example, graph_name)
    rows[[graph_name]] <- result$results
    if (keep_example)
      examples[[graph_name]] <- list(sampling = sampling, example = result$example)
  }
  results <- do.call(rbind, unname(rows))
  list(results = results, example = examples)
}

# One panel per statistic, with the established graph labels and curve style.
plot_supergraphs <- function(summary, statistic, graph_bandwidths, true_bandwidth,
                             trial_run = FALSE) {
  selected <- summary[summary$method == statistic, , drop = FALSE]
  graph_ids <- paste0("K=", graph_bandwidths)
  labels <- setNames(rep(NA_character_, length(graph_ids)), graph_ids)
  labels[graph_bandwidths == true_bandwidth] <- "Faithful Graph"
  # Number larger graphs by bandwidth; the sampling order is left unchanged.
  oversized_bandwidths <- sort(graph_bandwidths[graph_bandwidths > true_bandwidth])
  for (i in seq_along(oversized_bandwidths)) {
    graph_id <- paste0("K=", oversized_bandwidths[i])
    labels[graph_id] <- paste("Oversized Graph", i)
  }
  selected$method <- unname(labels[selected$graph])
  if (anyNA(selected$method)) stop("Summary contains a graph absent from the plotting configuration")
  plot_curve(selected, unname(labels), trial_run) +
    ggplot2::labs(color = "Graph", shape = "Graph", linetype = "Graph", linewidth = "Graph")
}
