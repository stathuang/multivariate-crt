#!/usr/bin/env Rscript
# Resumable experiment runner. Each repeat has explicit provenance and seed.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
source(file.path(dirname(sub("^--file=", "", file_arg[1])), "common.R"))
require_study_packages()
root <- study_root_from_script()
args <- commandArgs(trailingOnly = TRUE)
opts <- list(campaign = "corrected-v1", workers = "4", purpose = "production", backend = "psock")
i <- 1L
while (i <= length(args)) {
  if (args[i] %in% c("--list", "--retry-errors")) {
    opts[[substring(args[i], 3)]] <- TRUE; i <- i + 1L
  } else {
    if (!startsWith(args[i], "--") || i == length(args)) stop("Expected --option value.")
    opts[[substring(args[i], 3)]] <- args[i + 1L]; i <- i + 2L
  }
}
env <- load_study(root)
configs <- study_configs(env)
if (isTRUE(opts$list)) {
  for (config in configs) cat(config$id, config$family, config$repetitions, "repeats\n")
  quit(status = 0)
}
if (is.null(opts$scenario) || !opts$scenario %in% names(configs)) stop("Choose --scenario from --list.")
if (!grepl("^[A-Za-z0-9_-]+$", opts$campaign)) stop("Campaign must be a simple directory name.")
if (!opts$purpose %in% c("production", "timing")) stop("Purpose must be production or timing.")
if (!opts$backend %in% c("psock","serial")) stop("Backend must be psock or serial; native fork is unsupported.")
config <- configs[[opts$scenario]]
workers <- as.integer(opts$workers)
if (is.na(workers) || workers < 1 || workers > 14) stop("Workers must be between 1 and 14.")
if (config$family == "group_selection" && workers != 1L) {
  stop("Julia-backed group selection currently requires --workers 1 in a fresh R process.")
}
if (opts$backend == "serial" && workers != 1L) stop("Serial backend requires --workers 1.")
if (config$family == "group_selection") opts$backend <- "serial"
requested <- if (is.null(opts$limit)) as.integer(config$repetitions) else as.integer(opts$limit)
if (is.na(requested) || requested < 1 || requested > config$repetitions) stop("Invalid repeat limit.")
if (opts$purpose == "production" && requested != config$repetitions) {
  stop("Partial repeats require --purpose timing; never label a pilot as production.")
}
data_dir <- normalizePath(file.path(root, "data"), mustWork = TRUE)
run_dir <- file.path(root, "runs", opts$campaign, config$id)
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
lock <- file.path(run_dir, ".running")
if (!dir.create(lock, showWarnings = FALSE)) stop("Run lock exists. Inspect the recorded process before recovery: ", lock)
writeLines(c(paste("pid", Sys.getpid()), paste("started", Sys.time())), file.path(lock, "owner.txt"))

exit_status <- 0L
tryCatch({
  source_names <- study_source_names(config$family)
  julia_provenance <- NULL
  if (config$family == "group_selection") {
    source_names <- c(source_names, "scripts/runtime_setup.R")
    source(file.path(root, "scripts/runtime_setup.R"))
    configure_julia_runtime(root)
    julia_provenance <- julia_runtime_provenance(root)
  }
  vendor_files <- list.files(file.path(root, "vendor"), recursive = TRUE, full.names = FALSE)
  vendor_files <- vendor_files[grepl("[.]R$", vendor_files)]
  source_names <- c(source_names, file.path("vendor", vendor_files))
  source_hashes <- setNames(lapply(file.path(root, source_names), sha256_file), source_names)
  data_files <- if (config$family == "application") list.files(data_dir, recursive = TRUE, full.names = FALSE) else character()
  data_files <- data_files[!grepl("(^|/)[.]", data_files)]
  input_hashes <- setNames(lapply(file.path(data_dir, data_files), sha256_file), data_files)
  provenance <- list(config = config, sources = source_hashes, inputs = input_hashes,
                     R = R.version.string, platform = R.version$platform, packages = package_versions(),
                     RNG = RNGkind(), purpose = opts$purpose, julia = julia_provenance, backend = opts$backend)
  fingerprint <- sha256_object(provenance)
  manifest_path <- file.path(run_dir, "manifest.rds")
  if (file.exists(manifest_path)) {
    manifest <- readRDS(manifest_path)
    if (!identical(manifest$fingerprint, fingerprint)) stop("Code/configuration/input/environment differs; use a new campaign.")
    if (!identical(manifest$requested_repetitions, requested)) stop("Repeat limit differs; use a new timing campaign.")
  } else {
    snapshot <- file.path(run_dir, "source")
    for (name in source_names) {
      destination <- file.path(snapshot, name)
      dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
      if (!file.copy(file.path(root, name), destination)) stop("Could not freeze ", name)
    }
    atomic_json(config, file.path(snapshot, "config.json"))
    manifest <- list(schema_version = 2L, fingerprint = fingerprint, scenario = config$id,
       campaign = opts$campaign, purpose = opts$purpose, started = format(Sys.time(), tz = "UTC", usetz = TRUE),
       planned_repetitions = config$repetitions, requested_repetitions = requested,
       workers = workers, provenance = provenance, state = "created")
  }
  snapshot_hashes <- setNames(lapply(file.path(run_dir, "source", source_names), sha256_file), source_names)
  if (!identical(source_hashes, snapshot_hashes)) stop("Frozen source snapshot differs from recorded source hashes.")
  # Worker code is frozen. Inputs are hash-checked again before completion.
  source(file.path(run_dir,"source","scripts","worker.R"))
  env <- load_study(file.path(run_dir, "source"), family = config$family)
  options(mscrt.study_root = root)
  rep_dir <- file.path(run_dir, "replicates")
  dir.create(rep_dir, showWarnings = FALSE)
  rep_path <- function(rep) file.path(rep_dir, sprintf("replicate-%05d.rds", rep))
  valid_saved <- function(rep) {
    if (!file.exists(rep_path(rep))) return(FALSE)
    saved <- readRDS(rep_path(rep))
    if (!identical(saved$fingerprint, fingerprint) || !identical(saved$seed, replicate_seed(config, rep))) stop("Saved repeat provenance mismatch.")
    if (!identical(saved$replicate, rep) || !identical(saved$result_hash, sha256_object(saved$results))) stop("Saved repeat identity/content hash mismatch.")
    if (!identical(saved$diagnostics_hash,sha256_object(saved$diagnostics)) ||
        !identical(saved$generation_hash,sha256_object(saved$generation))) stop("Saved diagnostic/generation hash mismatch.")
    stopifnot(identical(saved$stock_audit_hash,sha256_object(saved$stock_audit)))
    validate_replicate(saved$results, config, rep)
    if (isTRUE(opts$`retry-errors`) && any(saved$results$status != "ok")) return(FALSE)
    TRUE
  }
  completed <- vapply(seq_len(requested), valid_saved, logical(1))
  pending <- which(!completed)
  manifest$state <- "running"; manifest$remaining <- length(pending)
  atomic_rds(manifest, manifest_path); atomic_json(manifest, file.path(run_dir, "manifest.json"))
  capture.output(sessionInfo(), file = file.path(run_dir, "sessionInfo.txt"))
  cat(config$id, "pending", length(pending), "of", requested, "workers", workers, "\n")
  source_root <- normalizePath(file.path(run_dir,"source"),mustWork=TRUE)
  cluster <- NULL
  if (opts$backend == "psock") cluster <- study_start_workers(workers,source_root,data_dir,config$family) else
    study_worker_initialize(source_root,data_dir,config$family)
  save_answer <- function(saved) {
    rep <- saved$replicate
    stopifnot(identical(saved$stock_audit_hash,sha256_object(saved$stock_audit)))
    validate_replicate(saved$results,config,rep)
    if (!identical(saved$result_hash,sha256_object(saved$results)) ||
        !identical(saved$diagnostics_hash,sha256_object(saved$diagnostics)) ||
        !identical(saved$generation_hash,sha256_object(saved$generation))) stop("Worker content hash mismatch.")
    if (file.exists(rep_path(rep))) {
      history <- file.path(run_dir,"retried",paste0("replicate-",rep,"-",format(Sys.time(),"%Y%m%dT%H%M%OS6"),".rds"))
      dir.create(dirname(history),recursive=TRUE,showWarnings=FALSE)
      if (!file.copy(rep_path(rep),history)) stop("Could not preserve failed attempt.")
    }
    atomic_rds(saved,rep_path(rep))
    cat(sprintf("repeat %d saved in %.2fs; failed cells %d\n",rep,saved$elapsed_seconds,sum(saved$results$status!="ok")))
  }
  # Small batches allow periodic progress manifests and a clean interrupt point.
  batches <- split(pending, ceiling(seq_along(pending) / max(1L, workers * 2L)))
  for (batch in batches) {
    tasks <- lapply(batch,function(rep) list(config=config,replicate=rep,
       seed=replicate_seed(config,rep),fingerprint=fingerprint))
    answers <- if (opts$backend == "psock") parallel::parLapply(cluster,tasks,study_compute_replicate) else
      lapply(tasks,study_compute_replicate)
    for (answer in answers) save_answer(answer)
    manifest$remaining <- sum(!vapply(seq_len(requested), valid_saved, logical(1)))
    manifest$updated <- format(Sys.time(), tz = "UTC", usetz = TRUE)
    atomic_rds(manifest, manifest_path); atomic_json(manifest, file.path(run_dir, "manifest.json"))
  }
  if (!is.null(cluster)) {parallel::stopCluster(cluster);cluster <- NULL}
  records <- lapply(seq_len(requested), function(rep) {
    saved <- readRDS(rep_path(rep))
    if (!identical(saved$fingerprint, fingerprint) || !identical(saved$replicate, rep) ||
        !identical(saved$seed, replicate_seed(config, rep)) ||
        !identical(saved$result_hash, sha256_object(saved$results))) stop("Result provenance/content mismatch before aggregation.")
    if (!identical(saved$diagnostics_hash,sha256_object(saved$diagnostics)) ||
        !identical(saved$generation_hash,sha256_object(saved$generation))) stop("Final diagnostic/generation hash mismatch.")
    stopifnot(identical(saved$stock_audit_hash,sha256_object(saved$stock_audit)))
    validate_replicate(saved$results, config, rep); saved
  })
  combined <- bind_result_rows(lapply(records, `[[`, "results"))
  atomic_csv(combined, file.path(run_dir, "results.csv"))
  atomic_rds(combined, file.path(run_dir, "results.rds"))
  final_data_files <- if (config$family == "application") list.files(data_dir, recursive = TRUE, full.names = FALSE) else character()
  final_data_files <- final_data_files[!grepl("(^|/)[.]", final_data_files)]
  if (!identical(data_files, final_data_files) ||
      !identical(input_hashes, setNames(lapply(file.path(data_dir, final_data_files), sha256_file), final_data_files))) {
    stop("Input data changed during run; outputs cannot be certified.")
  }
  seconds <- vapply(records, `[[`, numeric(1), "elapsed_seconds")
  manifest$state <- if (any(combined$status != "ok")) "finished_with_failures" else
    if (opts$purpose == "timing") "timing_completed" else "completed"
  manifest$remaining <- 0L
  manifest$finished <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  manifest$failed_cells <- sum(combined$status != "ok")
  manifest$timing <- list(median_repeat_seconds = median(seconds), max_repeat_seconds = max(seconds),
       total_worker_seconds = sum(seconds), projected_seconds_at_workers = median(seconds) * config$repetitions / workers)
  atomic_rds(manifest, manifest_path); atomic_json(manifest, file.path(run_dir, "manifest.json"))
  cat("Finished", config$id, manifest$state, "median repeat", median(seconds), "seconds\n")
  if (manifest$failed_cells > 0L) exit_status <- 2L
}, error = function(e) {
  message(conditionMessage(e)); exit_status <<- 1L
  atomic_json(list(state = "interrupted_or_failed", message = conditionMessage(e),
                  time = format(Sys.time(), tz = "UTC", usetz = TRUE)), file.path(run_dir, "runner_error.json"))
}, finally = {
  if (exists("cluster") && !is.null(cluster)) parallel::stopCluster(cluster)
  unlink(lock, recursive = TRUE)
})
quit(status = exit_status)
