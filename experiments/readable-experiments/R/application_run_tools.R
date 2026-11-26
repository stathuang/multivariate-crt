# Application identifiers and seed checks; input hashing remains outside tests.
prepare_application_config <- function(config) {
  config$methods <- c(config$statistics, config$baselines)
  config$graph_ids <- "fitted"
  if (length(config$repetitions) != 1L || !is.finite(config$repetitions) ||
      config$repetitions < 1 || config$repetitions != floor(config$repetitions))
    stop("repetitions must be a positive integer")
  if (config$design == "observed" && config$repetitions != 1L)
    stop("Observed analyses use exactly one repetition")
  seeds <- if (config$design == "observed") config$observed_seed else
    config$seed_base + c(1, config$repetitions)
  if (!length(seeds) || any(!is.finite(seeds) | seeds < 0 |
                            seeds > .Machine$integer.max | seeds != floor(seeds)))
    stop("Invalid application seeds")
  config
}

read_prepared_application <- function(path) {
  if (!file.exists(path)) stop("Prepared data missing. Run the corresponding prepare script first: ", path)
  list(data = readRDS(path), sha256 = digest::digest(file = path, algo = "sha256"))
}

# Observed applications have a single p-value per cell, with no simulation SE.
summarize_observed <- function(results) {
  if (anyDuplicated(report_key(results))) stop("Duplicate observed method results")
  out <- results[c(report_cell_keys, "pvalue", "status")]
  out$complete <- out$status == "ok" & is.finite(out$pvalue) & out$pvalue >= 0 & out$pvalue <= 1
  out
}
