# Capture warnings and failures without substituting successful-looking values.
capture_step <- function(expr) {
  warnings <- character()
  error <- ""
  condition <- NULL
  start <- proc.time()[["elapsed"]]
  value <- tryCatch(
    withCallingHandlers(force(expr), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }),
    error = function(e) {
      error <<- conditionMessage(e)
      condition <<- e
      NULL
    }
  )
  list(value = value, condition = condition, status = if (nzchar(error)) "failed" else "ok",
       error = error, warnings = paste(unique(warnings), collapse = " | "),
       elapsed_seconds = unname(proc.time()[["elapsed"]] - start))
}

join_warnings <- function(...) {
  messages <- c(...)
  paste(messages[nzchar(messages)], collapse = " | ")
}

test_row <- function(test, config, replicate, theta, method,
                     graph = paste0("K=", config$bandwidth)) {
  answer <- test$value
  if (test$status == "ok" &&
      (length(answer$pvalue) != 1L || !is.finite(answer$pvalue) ||
       answer$pvalue < 0 || answer$pvalue > 1)) {
    test$status <- "failed"
    test$error <- "Method returned an invalid p-value"
  }
  data.frame(
    scenario = config$id, replicate = as.integer(replicate), theta = theta,
    graph = graph, method = method,
    pvalue = if (test$status == "ok") answer$pvalue else NA_real_,
    statistic = if (test$status == "ok") answer$statistic else NA_real_,
    status = test$status, error = test$error, warnings = test$warnings,
    elapsed_seconds = test$elapsed_seconds, stringsAsFactors = FALSE
  )
}

# All output stays under the experiment directory; each run has its own folder.
new_output_directory <- function(experiment_dir, experiment_id, trial_run) {
  category <- if (trial_run) "trial" else "full"
  parent <- file.path(experiment_dir, "outputs", experiment_id, category)
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  stamp <- format(Sys.time(), "%Y%m%d-%H%M%S", tz = "Asia/Singapore")
  path <- tempfile(paste0(stamp, "-"), tmpdir = parent)
  if (!dir.create(path)) stop("Cannot create a fresh output directory")
  path
}

# Derive labels from the parameters actually used by the calculation.
prepare_gaussian_config <- function(config) {
  config$methods <- c(config$statistics, config$baselines)
  if (config$family != "supergraph") config$graph_bandwidths <- config$bandwidth
  if (!length(config$graph_bandwidths) || anyNA(config$graph_bandwidths) ||
      any(config$graph_bandwidths != floor(config$graph_bandwidths)) ||
      any(config$graph_bandwidths < config$bandwidth | config$graph_bandwidths >= config$p) ||
      anyDuplicated(config$graph_bandwidths))
    stop("Graph bandwidths must be distinct integers between the true bandwidth and p - 1")
  config$graph_ids <- paste0("K=", config$graph_bandwidths)
  if (length(config$repetitions) != 1L || !is.finite(config$repetitions) ||
      config$repetitions < 1 || config$repetitions != floor(config$repetitions))
    stop("repetitions must be a positive integer")
  seeds <- config$seed_base + c(1, config$repetitions)
  if (length(seeds) != 2L || any(!is.finite(seeds)) ||
      any(seeds < 0 | seeds > .Machine$integer.max | seeds != floor(seeds)))
    stop("Replicate seeds must be nonnegative R integers")
  config
}

# Evaluate the captured bytes, so recorded module hashes describe loaded code.
load_experiment_sources <- function(experiment_dir, envir, binary = FALSE,
                                    extensions = FALSE, oracle = FALSE, group = FALSE,
                                    stock = FALSE, bcr = FALSE, oracle_benchmark = FALSE) {
  files <- paste0("R/", c("sampler.R", "data.R", "linear_tests.R", "gaussian_tests.R",
                          "dcrt.R", "hdi_cache.R", "cgm_diagnostics.R", "comparators.R",
                          "run_tools.R", "experiment.R", "summary.R", "plot.R"))
  if (binary) files <- c(files, paste0("R/", c(
    "binary_data.R", "binary_likelihood.R", "result_diagnostics.R",
    "binary_tests.R", "binary_experiment.R"
  )))
  if (extensions) files <- c(files, paste0("R/", c(
    "extended_designs.R", "supergraph_experiment.R"
  )))
  if (oracle) files <- c(files, "R/oracle_dcrt.R")
  if (group) stop("Group selection is outside the submitted manuscript.")
  if (stock) files <- c(files, paste0("R/", c(
    "stock_data.R", "application_diagnostics.R", "application_tests.R",
    "application_run_tools.R", "application_tables.R", "stock_plot.R"
  )))
  if (bcr) files <- c(files, paste0("R/", c(
    "bcr_data.R", "multinomial_likelihood.R", "bcr_tests.R", "bcr_experiment.R",
    "result_diagnostics.R", "application_run_tools.R", "application_tables.R", "bcr_plot.R"
  )))
  if (oracle_benchmark) files <- c(files, paste0("R/", c(
    "oracle_common.R", "oracle_size.R", "oracle_power.R", "oracle_reporting.R"
  )))
  bytes <- setNames(lapply(file.path(experiment_dir, files), function(path)
    readBin(path, "raw", n = file.info(path)$size)), files)
  hashes <- vapply(bytes, digest::digest, "", algo = "sha256", serialize = FALSE)
  for (content in bytes) eval(parse(text = rawToChar(content)), envir = envir)
  list(bytes = bytes, sha256 = hashes)
}

computation_config <- function(config) {
  config$repetitions <- NULL
  config
}

capture_run_context <- function(config, loaded_sources) {
  packages <- c("MASS", "glmnet", "randomForest", "ggplot2", "scales", "digest", "hdi", "SIHR")
  if (config$response_family == "binomial") packages <- c(packages, "lpSolve", "gmp", "Rmpfr")
  if (config$family == "group_selection") packages <- c(packages, "JuliaCall", "knockoffsr")
  if (config$response_family == "multinomial")
    packages <- c(packages, "lpSolve", "gmp", "Rmpfr", "CVXR", "rags2ridges")
  first_seed <- if (!is.null(config$observed_seed)) config$observed_seed else config$seed_base + 1L
  list(config = computation_config(config), first_seed = first_seed,
       rng_kind = RNGkind(), R = R.version.string,
       packages = vapply(packages, function(p) as.character(utils::packageVersion(p)), ""),
       source_sha256 = loaded_sources$sha256)
}

assert_first_matches <- function(first, config, loaded_sources) {
  current <- capture_run_context(config, loaded_sources)
  if (is.null(first$context) || !identical(first$context, current))
    stop("Calculation configuration or loaded code changed; rerun section 3 before section 4.")
  invisible(TRUE)
}

start_run_record <- function(config, first, loaded_sources, output_dir, trial_run, alpha) {
  record <- c(list(state = "running", started_at = Sys.time(),
                  purpose = if (trial_run) "readability trial" else paste("full", config$id, "run"),
                  config = config, seeds = if (!is.null(config$observed_seed)) config$observed_seed else
                    config$seed_base + seq_len(config$repetitions),
                  summary_alpha = alpha), first$context[setdiff(names(first$context), "config")])
  # Save the loaded bytes, even if files on disk changed since section 2.
  for (file in names(loaded_sources$bytes)) {
    path <- file.path(output_dir, "source", file)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    writeBin(loaded_sources$bytes[[file]], path)
  }
  saveRDS(record, file.path(output_dir, "run.rds"))
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
  record
}

save_run_results <- function(record, results, output_dir) {
  record$state <- "computed"
  record$results <- results
  saveRDS(record, file.path(output_dir, "run.rds"))
  write_report_csv(results, file.path(output_dir, "results.csv"))
  record
}
