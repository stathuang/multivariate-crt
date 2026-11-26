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
