# Fresh-process execution for native numerical libraries; never fork an R runtime.
# Every worker reads frozen source and distinct tasks. Only the coordinator writes.
study_worker_initialize <- function(source_root, data_dir, family, expected_runtime = NULL) {
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  # Namespace initialization must precede explicit replicate seeding.
  actual_runtime <- list(R = R.version.string, platform = R.version$platform, packages = package_versions())
  if (!is.null(expected_runtime) && !identical(actual_runtime, expected_runtime))
    stop("Worker R/package inventory differs from coordinator provenance.")
  state <- new.env(parent = emptyenv())
  state$env <- load_study(source_root, family)
  state$data_dir <- data_dir
  assign(".study_worker_state", state, envir = .GlobalEnv)
  invisible(TRUE)

}
study_compute_replicate <- function(task) {
  state <- get(".study_worker_state", envir = .GlobalEnv)
  env <- state$env
  config <- task$config
  rep <- task$replicate
  before <- proc.time()[["elapsed"]]
  outer_warnings <- character()
  set.seed(task$seed)
  rows <- tryCatch(withCallingHandlers({
    answer <- if (config$family == "application") env$application_replicate(config, rep, state$data_dir) else
      env$simulation_replicate(config, rep)
    validate_replicate(answer, config, rep)
    answer

  }, warning = function(w) {
    outer_warnings <<- c(outer_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")

  }), error = function(e) {
    answer <- replicate_failure(config, rep, conditionMessage(e), proc.time()[["elapsed"]] - before)
    env$add_result_diagnostics(answer, condition = e)

  })
  if (length(outer_warnings)) rows$warnings <- paste(rows$warnings, paste(unique(outer_warnings), collapse = " | "), sep = " | ")
  stock_audit <- env$stock_audit_export()
  generation <- attr(rows, "generation")
  attr(rows, "generation") <- NULL
  diagnostics <- attr(rows, "diagnostics")
  attr(rows, "diagnostics") <- NULL
  list(fingerprint = task$fingerprint, seed = task$seed, replicate = rep,
    elapsed_seconds = proc.time()[["elapsed"]] - before,
    finished = format(Sys.time(), tz = "UTC", usetz = TRUE),
    result_hash = sha256_object(rows), results = rows, diagnostics = diagnostics,
    diagnostics_hash = sha256_object(diagnostics), generation = generation,
    generation_hash = sha256_object(generation), stock_audit = stock_audit,
    stock_audit_hash = sha256_object(stock_audit), rng_after = .Random.seed)

}
study_start_workers <- function(workers, source_root, data_dir, family) {
  expected_runtime <- list(R = R.version.string, platform = R.version$platform, packages = package_versions())
  cluster <- parallel::makePSOCKcluster(workers, rscript_args = "--vanilla", outfile = "")
  ready <- FALSE
  on.exit(if (!ready) parallel::stopCluster(cluster), add = TRUE)
  parallel::clusterCall(cluster, function(source_root, data_dir, family, expected_runtime) {
    sys.source(file.path(source_root, "scripts/common.R"), envir = .GlobalEnv)
    sys.source(file.path(source_root, "scripts/worker.R"), envir = .GlobalEnv)
    study_worker_initialize(source_root, data_dir, family, expected_runtime)

  }, source_root, data_dir, family, expected_runtime)
  ready <- TRUE
  cluster

}
