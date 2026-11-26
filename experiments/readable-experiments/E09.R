# E09: low-dimensional nonlinear binary response.
# Open Experiments.Rproj in RStudio, then run these sections from top to bottom.

# 1. Parameters ---------------------------------------------------------------
trial_run <- isTRUE(getOption("readable.trial", FALSE))  # Set TRUE for a two-repeat trial.
experiment_id <- "E09"
config <- list(
  id = "main_l_gl_m",
  family = "main",
  setting = "l_gl_m",
  n = 50L,
  p = 20L,
  # Test columns 1:8 conditional on the rest; theta scales transformed columns 1:10.
  # Columns 9:10 are conditioning variables; changing target does not move the signal block.
  target = c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L),
  theta = c(0L, 1L, 2L, 3L, 5L, 7L),
  repetitions = 800L,
  M = 100L,
  L = 3L,
  p_rule = "discrete",
  response_family = "binomial",
  statistics = c("RF", "RF-D", "RF-RR"),
  baselines = c("Chi-square", "dCRT (L1)", "dCRT (RF)"),
  seed_base = 700000L,
  bandwidth = 6L,
  precision_offdiag = 0.2,
  coefficient_policy = "redraw_each_replicate",
  # Keep the verified likelihood-supremum calculation and its numerical certificates.
  binary_likelihood_policy = "binary_face_exact_dual_v3"
)
config$methods <- c(config$statistics, config$baselines)
if (trial_run) config$repetitions <- 2L  # Keep n, p, M, L and all signal levels.
alpha <- 0.05
plot_methods <- c("RF", "RF-D", "RF-RR", "Chi-square", "dCRT (RF)")

# 2. Load shared functions ----------------------------------------------------
# Works from this RStudio project or from the repository root.
experiment_dir <- if (file.exists("Experiments.Rproj")) getwd() else
  file.path(getwd(), "experiments", "readable-experiments")
experiment_dir <- normalizePath(experiment_dir, mustWork = TRUE)
packages <- c("MASS", "glmnet", "randomForest", "ggplot2", "scales", "digest", "hdi", "SIHR",
              "lpSolve", "gmp", "Rmpfr")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))
source(file.path(experiment_dir, "R", "run_tools.R"), local = TRUE)
loaded_sources <- load_experiment_sources(experiment_dir, environment(), binary = TRUE)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

# 3. One replicate: inspect data, copies, statistics and p-values ---------------
config <- prepare_gaussian_config(config)
first_context <- capture_run_context(config, loaded_sources)
replicate <- 1L
set.seed(config$seed_base + replicate)
design <- capture_step(gaussian_design(config$n, config$p, config$bandwidth,
                                 config$precision_offdiag))
if (design$status == "ok") {
  X <- design$value$X
  graph <- design$value$graph
  beta <- design$value$beta  # Unused in this response; its draw preserves the RNG sequence.
  sampling <- capture_step(mscrt_sample(X, graph, config$target,
                                       M = config$M, L = config$L, reduced = TRUE))
} else {
  sampling <- design
}
first <- evaluate_binary(design, sampling, config, replicate, keep_example = TRUE)
first$context <- first_context
first$rng_after <- .Random.seed
Y <- first$example$Y  # First signal level; the full grid is in first$results.
copies <- sampling$value
tests <- first$example$tests
print(first$results[, c("theta", "method", "pvalue", "statistic", "status")])
# For example: tests[["RF"]]$value$diagnostics$randomized_statistics
# or: tests[["dCRT (L1)"]]$value$diagnostics$coordinate_pvalues

# 4. Repeat the experiment ----------------------------------------------------
config <- prepare_gaussian_config(config)
assert_first_matches(first, config, loaded_sources)
output_dir <- new_output_directory(experiment_dir, experiment_id, trial_run)
run <- start_run_record(config, first, loaded_sources, output_dir, trial_run, alpha)
replicate_results <- vector("list", config$repetitions)
replicate_results[[1L]] <- first$results
saveRDS(first, file.path(output_dir, "replicate-001.rds"))
saveRDS(list(design = design, sampling = sampling, example = first$example),
        file.path(output_dir, "example.rds"))

for (replicate in seq_len(config$repetitions)[-1L]) {
  set.seed(config$seed_base + replicate)
  design_next <- capture_step(gaussian_design(config$n, config$p, config$bandwidth,
                                        config$precision_offdiag))
  if (design_next$status == "ok") {
    sampling_next <- capture_step(mscrt_sample(
      design_next$value$X, design_next$value$graph, config$target,
      M = config$M, L = config$L, reduced = TRUE
    ))
  } else {
    sampling_next <- design_next
  }
  result <- evaluate_binary(design_next, sampling_next, config, replicate)
  result$rng_after <- .Random.seed
  replicate_results[[replicate]] <- result$results
  saveRDS(result, file.path(output_dir, sprintf("replicate-%03d.rds", replicate)))
  message("Completed replicate ", replicate, " / ", config$repetitions)
}
results <- bind_diagnostic_rows(replicate_results)
run <- save_run_results(run, results, output_dir)

# 5. Summarize rejection indicators and draw the historical panel ---
expected <- expand.grid(scenario = run$config$id, theta = run$config$theta,
                        graph = run$config$graph_ids, method = run$config$methods,
                        stringsAsFactors = FALSE)
expected$expected <- run$config$repetitions
summary <- summarize_results(results, expected, alpha = alpha)
write_report_csv(summary, file.path(output_dir, "counts.csv"))
plot <- plot_curve(summary, plot_methods, trial_run)
ggplot2::ggsave(file.path(output_dir, paste0("CRT-", run$config$setting, ".pdf")), plot, width = 8, height = 6)
ggplot2::ggsave(file.path(output_dir, paste0("CRT-", run$config$setting, ".png")), plot,
                width = 8, height = 6, dpi = 150)
run$summary_alpha <- alpha
run$plot_methods <- plot_methods
run$state <- "completed"
saveRDS(run, file.path(output_dir, "run.rds"))
print(summary[, c("theta", "method", "valid", "failed", "estimate", "se")])
message("Saved ", experiment_id, " outputs in: ", output_dir)
