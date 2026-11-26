# E11: one Gaussian design, the true graph and two graphs with additional edges.
# Open Experiments.Rproj in RStudio, then run these sections from top to bottom.

# 1. Parameters ---------------------------------------------------------------
trial_run <- isTRUE(getOption("readable.trial", FALSE))  # Set TRUE for a two-repeat trial.
experiment_id <- "E11"
config <- list(
  id = "supergraph_h_l_w",
  family = "supergraph",
  setting = "h_l_w",
  n = 80L,
  p = 120L,
  target = 1:8,
  theta = c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5),
  repetitions = 400L,
  M = 100L,
  L = 3L,
  p_rule = "discrete",
  response_family = "gaussian",
  statistics = c("LM-L1-R-SST", "LM-L1-R-SSR"),
  baselines = character(),
  seed_base = 1000000L,
  bandwidth = 6L,                    # True Gaussian design graph.
  precision_offdiag = 0.2,
  graph_bandwidths = c(6L, 12L, 18L), # Sampling graphs, in execution order.
  coefficient_policy = "redraw_each_replicate"
)
if (trial_run) config$repetitions <- 2L  # Keep n, p, M, L and all signal levels.
alpha <- 0.05

# 2. Load shared functions ----------------------------------------------------
experiment_dir <- if (file.exists("Experiments.Rproj")) getwd() else
  file.path(getwd(), "experiments", "readable-experiments")
experiment_dir <- normalizePath(experiment_dir, mustWork = TRUE)
packages <- c("MASS", "glmnet", "randomForest", "ggplot2", "scales", "digest", "hdi", "SIHR")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))
source(file.path(experiment_dir, "R", "run_tools.R"), local = TRUE)
loaded_sources <- load_experiment_sources(experiment_dir, environment(), extensions = TRUE)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

# 3. One replicate: inspect each graph and its own randomized copies ----------
config <- prepare_gaussian_config(config)
first_context <- capture_run_context(config, loaded_sources)
replicate <- 1L
set.seed(config$seed_base + replicate)
design <- capture_step(supergraph_design(config))
# Sampling and all signal levels are completed for K=6 before K=12, then K=18.
# Moving all sampling calls ahead of fitting would change the random-number stream.
first <- evaluate_supergraphs(design, config, replicate, keep_example = TRUE)
first$context <- first_context
first$rng_after <- .Random.seed
X <- design$value$X
graphs <- design$value$graphs
beta <- design$value$beta
examples <- first$example
Y <- examples[[1L]]$example$Y
copies <- examples[[1L]]$sampling$value  # First graph; examples holds all three.
tests <- examples[[1L]]$example$tests
print(first$results[, c("graph", "theta", "method", "pvalue", "statistic", "status")])

# 4. Repeat the experiment ----------------------------------------------------
config <- prepare_gaussian_config(config)
assert_first_matches(first, config, loaded_sources)
output_dir <- new_output_directory(experiment_dir, experiment_id, trial_run)
run <- start_run_record(config, first, loaded_sources, output_dir, trial_run, alpha)
replicate_results <- vector("list", config$repetitions)
replicate_results[[1L]] <- first$results
saveRDS(first, file.path(output_dir, "replicate-001.rds"))
saveRDS(list(design = design, example = examples), file.path(output_dir, "example.rds"))
for (replicate in seq_len(config$repetitions)[-1L]) {
  set.seed(config$seed_base + replicate)
  design_next <- capture_step(supergraph_design(config))
  result <- evaluate_supergraphs(design_next, config, replicate)
  result$rng_after <- .Random.seed
  replicate_results[[replicate]] <- result$results
  saveRDS(result, file.path(output_dir, sprintf("replicate-%03d.rds", replicate)))
  message("Completed replicate ", replicate, " / ", config$repetitions)
}
results <- do.call(rbind, replicate_results)
run <- save_run_results(run, results, output_dir)

# 5. Summarize by graph and draw one panel for each statistic ------------------
expected <- expand.grid(scenario = run$config$id, theta = run$config$theta,
                        graph = run$config$graph_ids, method = run$config$methods,
                        stringsAsFactors = FALSE)
expected$expected <- run$config$repetitions
summary <- summarize_results(results, expected, alpha = alpha)
write_report_csv(summary, file.path(output_dir, "counts.csv"))
for (statistic in run$config$statistics) {
  plot <- plot_supergraphs(summary, statistic, run$config$graph_bandwidths,
                          true_bandwidth = run$config$bandwidth, trial_run = trial_run)
  # Preserve the established SSR filename abbreviation.
  label <- if (statistic == "LM-L1-R-SSR") "LM-R-SSR" else statistic
  filename <- paste0("CRT-", run$config$setting, "-", label, "-OversizedGraphs")
  ggplot2::ggsave(file.path(output_dir, paste0(filename, ".pdf")), plot, width = 8, height = 6)
  ggplot2::ggsave(file.path(output_dir, paste0(filename, ".png")), plot, width = 8, height = 6, dpi = 150)
}
run$state <- "completed"
saveRDS(run, file.path(output_dir, "run.rds"))
print(summary[, c("graph", "theta", "method", "valid", "failed", "estimate", "se")])
message("Saved ", experiment_id, " outputs in: ", output_dir)
