# Portable loading, provenance and atomic output helpers.
study_root_from_script <- function() {
  file_arg <- grep("^--file=", commandArgs(), value = TRUE)
  if (!length(file_arg)) stop("Run this entry point with Rscript.")
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."), mustWork = TRUE)
}

require_study_packages <- function() {
  needed <- c("jsonlite", "digest")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install documented dependencies first: ", paste(missing, collapse = ", "))
}

study_module_names <- function(family=NULL) {
  files <- c("core.R", "statistics.R", "binary_likelihood.R", "cgm_diagnostics.R",
             "result_diagnostics.R", "multinomial_likelihood.R", "baselines.R")
  if (is.null(family) || family == "application") files <- c(files,"applications.R","stock_audit.R")
  if (is.null(family) || family != "application") files <- c(files,"baseline_cache_candidate.R","simulation.R")
  files
}

study_source_names <- function(family) {
  c(file.path("R",study_module_names(family)),"scripts/common.R","scripts/run.R","scripts/worker.R")
}

load_study <- function(root, family = NULL) {
  options(mscrt.study_root = root)
  env <- new.env(parent = globalenv())
  files <- study_module_names(family)
  for (name in files) {
    path <- file.path(root, "R", name)
    if (!file.exists(path)) stop("Missing implementation: ", path)
    sys.source(path, envir = env)
  }
  env
}

sha256_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
sha256_object <- function(x) digest::digest(x, algo = "sha256", serialize = TRUE)

atomic_rds <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = ".write-", tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(object, temporary, version = 3)
  if (!file.rename(temporary, path)) stop("Could not atomically replace ", path)
}

atomic_json <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = ".write-", tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  jsonlite::write_json(object, temporary, pretty = TRUE, auto_unbox = TRUE,
                       null = "null", na = "null", digits = 16)
  if (!file.rename(temporary, path)) stop("Could not atomically replace ", path)
}

atomic_csv <- function(object, path) {
  temporary <- tempfile(pattern = ".write-", tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  # Explicit 17 significant digits allow a double to round-trip through CSV.
  for (name in names(object)) if (is.double(object[[name]])) {
    value <- object[[name]]
    object[[name]] <- ifelse(is.na(value), NA_character_, sprintf("%.17g", value))
  }
  write.csv(object, temporary, row.names = FALSE, na = "")
  if (!file.rename(temporary, path)) stop("Could not atomically replace ", path)
}

package_versions <- function() {
  pkgs <- c("jsonlite", "digest", "MASS", "Matrix", "mvtnorm", "glmnet", "randomForest",
            "glasso", "bestNormalize", "VGAM", "hdi", "SIHR", "rags2ridges", "ggplot2", "knockoffsr", "JuliaCall", "lpSolve", "gmp", "Rmpfr", "CVXR", "ECOSolveR", "igraph", "scales")
  setNames(lapply(pkgs, function(p) {
    if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else NA_character_
  }), pkgs)
}

study_configs <- function(env, root = getOption("mscrt.study_root")) {
  sims <- env$simulation_configs(file.path(root, "config", "simulations.json"))
  apps <- env$application_configs(file.path(root, "config", "applications.json"))
  configs <- c(sims, apps)
  ids <- vapply(configs, `[[`, character(1), "id")
  if (anyDuplicated(ids)) stop("Duplicate scenario identifiers.")
  names(configs) <- ids
  configs
}

result_key <- function(x) paste(x$theta, x$graph, x$method, sep = "\r")

expected_rows <- function(config) {
  if (is.null(config$methods) || is.null(config$graph_ids)) {
    stop("Scenario lacks explicit methods/graph_ids: ", config$id)
  }
  expand.grid(theta = unlist(config$theta), graph = unlist(config$graph_ids),
              method = unlist(config$methods), stringsAsFactors = FALSE)
}

validate_replicate <- function(rows, config, replicate) {
  required <- c("scenario", "replicate", "theta", "graph", "method", "pvalue", "status",
                "error", "warnings", "elapsed_seconds")
  if (!is.data.frame(rows) || !all(required %in% names(rows))) stop("Malformed replicate result.")
  if (!nrow(rows) || any(rows$scenario != config$id) || any(rows$replicate != replicate)) {
    stop("Replicate identifier mismatch.")
  }
  if (anyDuplicated(result_key(rows))) stop("Duplicate method/cell result.")
  if (!setequal(result_key(rows), result_key(expected_rows(config)))) stop("Incomplete or unexpected result grid.")
  if (anyNA(rows$status) || any(!rows$status %in% c("ok", "error", "failed", "skipped"))) stop("Unknown result status.")
  ok <- rows$status == "ok"
  if (!identical(config$family, "group_selection")) {
    if (any(!is.finite(rows$pvalue[ok])) || any(rows$pvalue[ok] < 0 | rows$pvalue[ok] > 1)) {
      stop("Successful test must have a finite p-value in [0,1].")
    }
  } else {
    if (!all(c("FDP", "power") %in% names(rows)) ||
        any(!is.finite(rows$FDP[ok])) || any(!is.finite(rows$power[ok])) ||
        any(rows$FDP[ok] < 0 | rows$FDP[ok] > 1) || any(rows$power[ok] < 0 | rows$power[ok] > 1)) {
      stop("Successful selection result must have finite FDP/power in [0,1].")
    }
  }
  if (any(!is.finite(rows$elapsed_seconds)) || any(rows$elapsed_seconds < 0)) stop("Invalid elapsed time.")
  invisible(TRUE)
}

replicate_failure <- function(config, replicate, message, seconds, warnings = "") {
  rows <- expected_rows(config)
  rows$scenario <- config$id; rows$replicate <- replicate
  rows$pvalue <- NA_real_; rows$status <- "error"; rows$error <- message
  rows$warnings <- warnings; rows$elapsed_seconds <- seconds
  if (identical(config$family, "group_selection")) {
    rows$FDP <- NA_real_; rows$power <- NA_real_
  }
  rows
}

replicate_seed <- function(config, replicate) {
  seed <- if (!is.null(config$observed_seed)) config$observed_seed else config$seed_base + replicate
  if (length(seed) != 1 || !is.finite(seed) || seed < 0 || seed > .Machine$integer.max) stop("Invalid seed.")
  as.integer(seed)
}

bind_result_rows <- function(results) {
  cols <- unique(unlist(lapply(results, names)))
  do.call(rbind, lapply(results, function(x) {
    for (name in setdiff(cols, names(x))) x[[name]] <- NA
    x[, cols, drop = FALSE]
  }))
}
