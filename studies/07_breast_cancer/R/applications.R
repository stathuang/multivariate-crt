# Application replicate: data generation, randomisation, and method evaluation.
.bcr_lrt_nonconvergence <- function(warnings) {
  any(grepl("convergence not obtained|failed to converge|did not converge|iteration limit.*reached|maximum number of iterations.*reached",
            warnings, ignore.case = TRUE))
}

application_replicate <- function(config, replicate, data_dir) {
  if (length(replicate) != 1L || !is.finite(replicate) || replicate != as.integer(replicate) ||
      replicate < 1L || replicate > config$repetitions) stop("Invalid application replicate")
  path <- file.path(data_dir, config$data_file)
  if (!file.exists(path)) stop("Prepared application data missing; explicitly call prepare_application_data() first")
  data <- readRDS(path)
  is_stock <- config$dataset == "stock"
  policy <- config$numerical_status_policy
  if (!is_stock && !is.null(policy) && !identical(policy, "bcr_lrt_nonconvergence_is_error_v1"))
    stop("Unknown BCR numerical status policy")
  target <- if (is_stock) data$targets[[config$group]] else data$target
  methods <- c(config$statistics, config$baselines)
  stage_warnings <- character(); start <- proc.time()[["elapsed"]]
  generated <- tryCatch(withCallingHandlers({
    X <- data$X; y <- data$y
    if (config$design == "stock_fixed_start_null") {
      # Historical finite kernel from one fixed observed matrix, no stationarity claim.
      X <- mscrt_sample(X, data$graph, target = seq_len(ncol(X)), M = 1L, L = 20L,
                        reduced = FALSE)$X_copies[[1L]]
      null <- data$nulls[[config$group]]
      y <- null$coef[1L] + as.numeric(X[, null$outside, drop = FALSE] %*% null$coef[-1L]) +
        stats::rnorm(nrow(X), sd = null$sd)
    } else if (config$design == "bcr_fitted_null") {
      null <- data$null
      X <- rags2ridges::rmvnormal(nrow(X), mu = null$mean, sigma = null$covariance)
      prob <- stats::predict(null$fit, newx = X[, null$outside, drop = FALSE],
                             type = "response", s = null$lambda)[, , 1L]
      draw <- apply(prob, 1L, function(v) which(stats::rmultinom(1L, 1L, v)[, 1L] == 1L))
      y <- factor(draw, levels = seq_len(ncol(prob)))
    } else if (config$design == "stock_added_noise") {
      y <- y + stats::rnorm(length(y), sd = config$theta)
    } else if (config$design != "observed") stop("Unknown application design")
    list(X = X, y = y)
  }, warning = function(w) { stage_warnings <<- c(stage_warnings, conditionMessage(w)); invokeRestart("muffleWarning") }),
    error = function(e) e)
  # A failed shared sampler marks G-CRT methods failed but does not erase baselines.
  sample_warnings <- character()
  copies <- if (!inherits(generated, "error")) tryCatch(withCallingHandlers(
    mscrt_sample(generated$X, data$graph, target = target, M = config$M, L = config$L, reduced = TRUE),
    warning = function(w) { sample_warnings <<- c(sample_warnings, conditionMessage(w)); invokeRestart("muffleWarning") }),
    error = function(e) e) else generated
  shared_elapsed <- proc.time()[["elapsed"]] - start
  rows <- lapply(methods, function(method) {
    warns <- stage_warnings; method_warnings <- character(); t0 <- proc.time()[["elapsed"]]
    result <- tryCatch(withCallingHandlers({
      if (inherits(generated, "error")) stop(conditionMessage(generated))
      if (method %in% config$statistics) {
        warns <- c(warns, sample_warnings)
        if (inherits(copies, "error")) stop(conditionMessage(copies))
        gcrt_test(generated$X, generated$y, target, copies, statistic = method,
                  response_family = config$response_family, p_rule = config$p_rule)
      } else if (!is_stock && method %in% c("dCRT (L1)", "dCRT (RF)")) {
        # Historical BCR dCRT tests use numeric class codes and gaussian nuisance fits.
        baseline_test(generated$X, as.numeric(generated$y), target, method, response_family = "gaussian")
      } else baseline_test(generated$X, generated$y, target, method, response_family = config$response_family)
    }, warning = function(w) {
      method_warnings <<- c(method_warnings, conditionMessage(w))
      warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") }),
      error = function(e) e)
    error <- if (inherits(result, "error")) conditionMessage(result) else ""
    if (!nzchar(error) && (length(result) != 1L || !is.finite(result) || result < 0 || result > 1))
      error <- "Method returned an invalid p-value"
    raw_pvalue <- if (is.numeric(result) && length(result) == 1L && is.finite(result)) as.numeric(result) else NA_real_
    if (!is_stock && method == "LRT" && identical(policy, "bcr_lrt_nonconvergence_is_error_v1") &&
        .bcr_lrt_nonconvergence(method_warnings) && !nzchar(error))
      error <- "BCR LRT fit did not converge; numeric return retained only as raw_pvalue"
    row <- data.frame(scenario = config$id, replicate = replicate, theta = config$theta, graph = "fitted",
               method = method, pvalue = if (nzchar(error)) NA_real_ else as.numeric(result),
               status = if (nzchar(error)) "error" else "ok", error = error,
               warnings = paste(unique(warns), collapse = " | "),
               elapsed_seconds = proc.time()[["elapsed"]] - t0,
               shared_elapsed_seconds = shared_elapsed, dataset = config$dataset, design = config$design,
               panel = if (is_stock) "stocks" else data$panel, stringsAsFactors = FALSE)
    if (!is_stock) {
      row$raw_pvalue <- if (method == "LRT") raw_pvalue else NA_real_
      row$numerical_status_policy <- if (is.null(policy)) "legacy_returned_value" else policy
    }
    row
  })
  do.call(rbind, rows)
}
