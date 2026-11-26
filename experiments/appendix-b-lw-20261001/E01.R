# E01: original comparisons plus the D.10 Ledoit-Wolf plug-in method.
# 1. Parameters ---------------------------------------------------------------
trial_run <- isTRUE(getOption("readable.trial", FALSE))
experiment_id <- "E01"
config <- list(id = "oracle_plugin_benchmark",
  family = "oracle_benchmark", response_family = "gaussian",
  master_seed = 20260429L, repetitions = 5000L, M = 999L, lambda_ridge = 1,
  workers = getOption("lw.workers", 1L))
if (trial_run) config$repetitions <- 2L  # Retain every setting, signal level and copy.
config$repetitions <- getOption("lw.repetitions", config$repetitions)
alpha <- .05
settings <- data.frame(setting = c("LD-tight", "LD-moderate", "p>n target"),
  n = c(60L, 80L, 80L), s = c(45L, 50L, 60L), t = c(10L, 10L, 40L))
settings$p <- settings$s + settings$t
settings$d <- settings$n - settings$s - 1L
settings$panel <- sprintf("n=%d, s=%d, t=%d, d=%d", settings$n, settings$s, settings$t, settings$d)
config$settings <- settings
config$methods <- c("Oracle CRT", "Plug-in MLE CRT", "Plug-in CRT (d.f. adjustment)", "MVN-CRT", "LW plug-in CRT (D.10)")
# The first actual replicate seed, used by the common configuration guard.
config$seed_base <- config$master_seed + 100000L
if (.Platform$OS.type == "windows" && config$workers > 1L) {
  warning("Parallel execution is unavailable with mclapply on Windows; using one worker.")
  config$workers <- 1L
}
if (any(!is.finite(c(config$repetitions, config$M, config$workers))) ||
    any(c(config$repetitions, config$M, config$workers) < 1) ||
    any(c(config$repetitions, config$M, config$workers) != floor(c(config$repetitions, config$M, config$workers))))
  stop("Repetitions, randomizations and workers must be positive integers")

# 2. Load shared mathematical and run tools -----------------------------------
experiment_dir <- if (file.exists("Experiments.Rproj")) getwd() else
  file.path(getwd(), "experiments", "appendix-b-lw-20261001")
experiment_dir <- normalizePath(experiment_dir, mustWork = TRUE)
source(file.path(experiment_dir, "R/run_tools.R"), local = TRUE)
source(file.path(experiment_dir, "R/lw_run_tools.R"), local = TRUE)
loaded_sources <- load_lw_sources(experiment_dir, environment(), paste0(experiment_id, ".R"))
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

# 3. Inspect the first dataset and all five randomizations ---------------------
first_context <- capture_run_context(config, loaded_sources)
coefs <- size_coefficients(settings$s[1], settings$t[1], 1L)
first <- capture_step(size_replication_lw(1L, settings[1,], coefs,
  setting_index = 1L,
  seed_base = config$master_seed, M = config$M, lambda_ridge = config$lambda_ridge, keep_example = TRUE))
first$context <- first_context
first$rng_after <- .Random.seed
X <- first$value$diagnostics$data$X
Y <- first$value$diagnostics$data$Y
tests <- first$value
print(tests$pvalue)
print(first$elapsed_seconds)

# 4. Repeat over settings, retaining the first result -----------
config$settings <- settings
assert_first_matches(first, config, loaded_sources)
output_dir <- new_output_directory(experiment_dir, experiment_id, trial_run)
run <- start_run_record(config, first, loaded_sources, output_dir, trial_run, alpha)
seed_plan <- expand.grid(replicate = seq_len(config$repetitions),
  setting_index = seq_len(nrow(settings)))
seed_plan$seed <- config$master_seed + 100000L * seed_plan$setting_index + seed_plan$replicate
run$seeds <- seed_plan
saveRDS(run, file.path(output_dir, "run.rds"))
saveRDS(first, file.path(output_dir, "example.rds"))
all_rows <- list()
for (i in seq_len(nrow(settings))) {
  cfg <- settings[i,]
  coefs <- size_coefficients(cfg$s, cfg$t, i)
  replicate_ids <- seq_len(config$repetitions)
  if (i == 1L) replicate_ids <- replicate_ids[-1L]
  calculate <- function(b) {
    value <- capture_step(size_replication_lw(b, cfg, coefs,
      setting_index = i,
      seed_base = config$master_seed, M = config$M, lambda_ridge = config$lambda_ridge))
    value$rng_after <- .Random.seed
    value
  }
  # Keep the original optional process parallelism and explicit per-replicate seeds.
  values <- if (config$workers > 1L) parallel::mclapply(replicate_ids, calculate,
    mc.cores = config$workers, mc.preschedule = TRUE, mc.set.seed = FALSE) else
    lapply(replicate_ids, calculate)
  names(values) <- as.character(replicate_ids)
  if (i == 1L) values <- c(list("1" = first), values)
  for (b in seq_len(config$repetitions)) {
    value <- values[[as.character(b)]]
    if (inherits(value, "try-error")) stop("A parallel worker failed: ", value)
    saveRDS(value, file.path(output_dir, sprintf("setting-%02d-replicate-%05d.rds", i, b)))
    rows <- expand.grid(replication = b, statistic = c("raw score", "poly2 score", "ridge gain"),
                        method = config$methods, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    rows$setting <- cfg$setting
    rows$n <- cfg$n; rows$s <- cfg$s; rows$t <- cfg$t; rows$p <- cfg$p; rows$d <- cfg$d
    rows$p_value <- if (value$status == "ok") as.vector(value$value$pvalue) else NA_real_
    rows$observed_statistic <- if (value$status == "ok") unname(value$value$statistic[rows$statistic]) else NA_real_
    rows$status <- value$status
    rows$error <- value$error
    rows$warnings <- value$warnings
    rows <- mark_lw_failures(rows, value)
    all_rows[[length(all_rows) + 1L]] <- rows
  }
  message(cfg$setting, ": completed ", config$repetitions, " replicates")
}
results <- do.call(rbind, all_rows)
run <- save_run_results(run, results, output_dir)

# 5. Original summary convention and output style -----------------------------
summary <- summarize_oracle(results, config$repetitions, config$M, alpha, settings,
                             config$methods, power = FALSE)
write_report_csv(summary, file.path(output_dir, "counts.csv"))
write_report_csv(results, file.path(output_dir, "oracle_plugin_benchmark_pvalues.csv"))
write_report_csv(summary, file.path(output_dir, "oracle_plugin_benchmark_summary.csv"))
write_lw_size_table(summary, settings, config$methods,
  file.path(output_dir, "oracle_plugin_benchmark_size.tex"), trial_run)
run$state <- "completed"
saveRDS(run, file.path(output_dir, "run.rds"))
message("Saved ", experiment_id, " outputs in: ", output_dir)
