# Existing stock and breast-cancer applications. No downloads or legacy sourcing.
# Callers own replicate RNG; only explicit preparation owns fixed fitting seeds.
application_configs <- function(path = file.path("config", "applications.json")) {
  x <- jsonlite::fromJSON(path, simplifyVector = FALSE)$scenarios
  x <- lapply(x, function(z) {
    for (nm in c("statistics", "baselines", "theta", "target", "methods", "graph_ids"))
      if (!is.null(z[[nm]])) z[[nm]] <- unlist(z[[nm]], use.names = FALSE)
    z
  })
  stats::setNames(x, vapply(x, `[[`, "", "id"))
}

.application_seed <- function(seed, code) {
  existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (existed) old <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (existed) assign(".Random.seed", old, envir = .GlobalEnv) else
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv))
  set.seed(seed)
  force(code)
}

.application_hash <- function(path) digest::digest(file = path, algo = "sha256")

read_geo_expression <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  lines <- readLines(con, warn = FALSE)
  first <- which(lines == "!series_matrix_table_begin")
  last <- which(lines == "!series_matrix_table_end")
  if (length(first) != 1L || length(last) != 1L || last <= first + 2L)
    stop("Expected exactly one GEO expression table")
  tab <- utils::read.delim(text = paste(lines[seq.int(first + 1L, last - 1L)], collapse = "\n"),
                          check.names = FALSE, stringsAsFactors = FALSE)
  ids <- tab[[1L]]
  if (anyDuplicated(ids) || anyDuplicated(names(tab)[-1L])) stop("Duplicate GEO identifiers")
  X <- t(as.matrix(tab[-1L])); storage.mode(X) <- "double"
  colnames(X) <- ids
  list(X = X, samples = rownames(X), probes = ids)
}

bcr_endpoint <- function(clinical) {
  required <- c("GEO.asscession.number", "time.to.relapse.or.last.follow.up..months.",
                "relapse..1.True.", "ER.Status", "Brain.relapses..1.yes..0.no.")
  if (!all(required %in% names(clinical))) stop("Clinical column names differ from documented GSE2034 schema")
  months <- clinical[[required[2L]]]; relapse <- clinical[[required[3L]]]
  brain <- clinical[[required[5L]]]
  if (anyNA(clinical[required]) || any(!relapse %in% 0:1) || any(!brain %in% 0:1) ||
      any(months < 0) || anyDuplicated(clinical[[required[1L]]])) stop("Invalid clinical values")
  if (any(months <= 48 & relapse == 0))
    stop("Early censored patients require an explicit endpoint decision; cannot silently label them")
  y <- integer(nrow(clinical))
  y[months <= 48 & relapse == 1] <- 1L
  y[months <= 48 & relapse == 1 & brain == 1] <- 2L
  factor(y, levels = 0:2)
}

prepare_bcr <- function(input_dir, panel = c("historical_expression_first60", "wang_erplus60"),
                        selection_seed = 2023L, erplus_probes = NULL) {
  panel <- match.arg(panel)
  clinical <- utils::read.csv(file.path(input_dir, "GSE2034-PI.csv"), fileEncoding = "UTF-8-BOM")
  expression <- read_geo_expression(file.path(input_dir, "GSE2034_series_matrix.txt"))
  listed <- utils::read.csv(file.path(input_dir, "ImpVar.csv"), header = FALSE,
                           fileEncoding = "UTF-8-BOM")[[1L]]
  if (length(listed) != 76L || anyDuplicated(listed)) stop("Expected 76 unique listed probes")
  ordered <- expression$probes[expression$probes %in% listed]
  if (length(ordered) != 76L) stop("A listed probe is absent")
  if (panel == "wang_erplus60") {
    if (is.null(erplus_probes) || length(erplus_probes) != 60L || anyDuplicated(erplus_probes) ||
        !identical(as.character(erplus_probes), as.character(listed[1:60])))
      stop("Corrected panel requires the verified Wang Table 3 ER-positive probe IDs in exact order")
  }
  probes <- if (panel == "historical_expression_first60") ordered[1:60] else erplus_probes
  join <- match(clinical$GEO.asscession.number, expression$samples)
  if (anyNA(join) || anyDuplicated(join)) stop("Clinical/expression join is not one to one")
  # Original orderNorm is fitted on ALL 286 patients before restriction to ER+.
  raw <- expression$X[join, probes, drop = FALSE]
  if (any(!is.finite(raw))) stop("Non-finite selected expression data")
  normalized <- apply(raw, 2L, function(v) bestNormalize::orderNorm(v)$x.t)
  erplus <- clinical$ER.Status == "ER+"
  y <- bcr_endpoint(clinical)[erplus]
  X <- normalized[erplus, , drop = FALSE]
  selection <- .application_seed(selection_seed, {
    cv <- glmnet::cv.glmnet(X, y, family = "multinomial", nfolds = 5L)
    fit <- glmnet::glmnet(X, y, family = "multinomial")
    selected <- function(s) which(rowSums(as.matrix(do.call(cbind, stats::coef(fit, s = s)))^2)[-1L] != 0)
    list(strong = selected(cv$lambda.min), condition = selected(cv$lambda.1se),
         lambda_min = cv$lambda.min, lambda_1se = cv$lambda.1se, foldid = cv$foldid)
  })
  target_panel <- setdiff(selection$strong, selection$condition)
  if (!length(target_panel)) stop("Selected target is empty")
  if (!all(selection$condition %in% selection$strong))
    stop("One-SE selected genes are not nested; hypothesis requires explicit resolution")
  mat_X <- unname(X[, selection$strong, drop = FALSE])
  target <- match(target_panel, selection$strong)
  graph <- matrix(1L, ncol(mat_X), ncol(mat_X)); diag(graph) <- 0L
  outside <- setdiff(seq_len(ncol(mat_X)), target)
  null <- .application_seed(2023L, {
    cv <- glmnet::cv.glmnet(mat_X[, outside, drop = FALSE], y, family = "multinomial", nfolds = 5L)
    fit <- glmnet::glmnet(mat_X[, outside, drop = FALSE], y, family = "multinomial")
    list(fit = fit, lambda = cv$lambda.min, seed = 2023L,
         mean = colMeans(mat_X), covariance = stats::cov(mat_X), outside = outside)
  })
  list(X = mat_X, y = y, graph = graph, target = target, null = null,
       panel = panel, panel_probes = probes, selected_probes = probes[selection$strong],
       condition_probes = probes[selection$condition], target_probes = probes[target_panel],
       selection = selection, class_counts = table(y), selection_seed = selection_seed,
       endpoint = c("0" = "No relapse by 48 months", "1" = "Relapse by 48 months, brain-relapse flag absent",
                    "2" = "Relapse by 48 months, brain-relapse flag present"),
       endpoint_limitation = "Brain flag has no separate event date; not proof brain relapse occurred by month 48",
       scope = "Exploratory outcome-selected hypothesis; fixed preprocessing and panel in null simulation")
}

prepare_stock <- function(input_dir, null_seed = 2023L) {
  e <- new.env(parent = emptyenv()); load(file.path(input_dir, "DataForExp.RData"), envir = e)
  pre <- new.env(parent = emptyenv()); load(file.path(input_dir, "Pre-X.RData"), envir = pre)
  names <- utils::read.csv(file.path(input_dir, "list_stocks.csv"))$Symbol
  X <- unname(e$mat_X); y <- as.numeric(e$Y); graph <- unname(e$G1 * 1L)
  if (!identical(dim(X), c(92L, 103L)) || length(y) != nrow(X) ||
      any(!is.finite(X)) || any(!is.finite(y)) || length(names) != ncol(X)) stop("Unexpected saved stock input")
  if (max(abs(X - pre$X)) != 0) stop("Saved stock X disagrees with preprocessing snapshot")
  within_sector <- outer(as.character(pre$industry), as.character(pre$industry), `==`) * 1L
  penalty <- (1 - within_sector) * 0.15; diag(penalty) <- 0
  fitted <- (glasso::glasso(stats::var(X), rho = penalty)$wi != 0) * 1L; diag(fitted) <- 0L
  if (!identical(dim(fitted), dim(graph)) || any(fitted != graph)) stop("Reconstructed graph differs from saved graph")
  targets <- list(T1 = as.integer(e$setT1), T2 = as.integer(e$setT2))
  nulls <- lapply(targets, function(target) .application_seed(null_seed, {
    outside <- setdiff(seq_len(ncol(X)), target)
    fit <- glmnet::cv.glmnet(X[, outside, drop = FALSE], y)
    fitted_y <- stats::predict(fit, newx = X[, outside, drop = FALSE], s = fit$lambda.min)
    list(coef = as.numeric(stats::coef(fit, s = fit$lambda.min)),
         sd = stats::sd(y - fitted_y), lambda = fit$lambda.min, outside = outside, seed = null_seed)
  }))
  list(X = X, y = y, graph = graph, targets = targets, nulls = nulls, stocks = names,
       target_names = lapply(targets, function(t) names[t]), sectors = as.character(pre$industry),
       start_date = as.character(pre$start_date), end_date = as.character(pre$end_date),
       starting_point = "Saved transformed covariates and saved aligned DIA response; raw alignment not reconstructed",
       null_seed_status = "New explicit preparation seed; historical stock null fitting seed was not recorded")
}

# Explicit import: source_dir has Stock/ and Breast/; all subsequent reads use data_dir.
# External data and derivatives remain excluded from the distributable public folder.
prepare_application_data <- function(source_dir, data_dir, overwrite = FALSE) {
  source_dir <- normalizePath(source_dir, mustWork = TRUE)
  # Require an existing package data directory so validation never creates a source child.
  data_dir <- normalizePath(data_dir, mustWork = TRUE)
  if (basename(data_dir) != "data" || !file.exists(file.path(dirname(data_dir), "R", "applications.R")))
    stop("Preparation destination must be the data directory of this reproducibility package")
  if (identical(source_dir, data_dir) || startsWith(source_dir, paste0(data_dir, "/")) ||
      startsWith(data_dir, paste0(source_dir, "/")))
    stop("Source must be separate from destination")
  inputs <- c("Stock/DataForExp.RData", "Stock/Pre-X.RData", "Stock/list_stocks.csv",
              "Breast/GSE2034-PI.csv", "Breast/GSE2034_series_matrix.txt", "Breast/ImpVar.csv")
  cache <- file.path(data_dir, "external-cache")
  dest <- file.path(cache, "inputs", inputs)
  if (any(!file.exists(file.path(source_dir, inputs)))) stop("Required input missing; see data/README.md")
  if (!overwrite && any(file.exists(dest))) stop("Cache already exists; use an empty destination or explicit overwrite")
  for (d in unique(dirname(dest))) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  if (!all(file.copy(file.path(source_dir, inputs), dest, overwrite = overwrite))) stop("Input copy failed")
  if (!all(vapply(dest, .application_hash, "") == vapply(file.path(source_dir, inputs), .application_hash, "")))
    stop("Import hash mismatch")
  prepared <- file.path(cache, "prepared"); dir.create(prepared, showWarnings = FALSE)
  prep_warnings <- character()
  withCallingHandlers({
    stock <- prepare_stock(file.path(cache, "inputs", "Stock"))
    saveRDS(stock, file.path(prepared, "stock.rds"), version = 3L)
    historical <- prepare_bcr(file.path(cache, "inputs", "Breast"), "historical_expression_first60")
    saveRDS(historical, file.path(prepared, "bcr_historical.rds"), version = 3L)
    verified_probes <- utils::read.csv(file.path(data_dir, "public", "wang2005_erplus60.csv"),
                                      stringsAsFactors = FALSE)$probe_id
    corrected <- prepare_bcr(file.path(cache, "inputs", "Breast"), "wang_erplus60",
                             erplus_probes = verified_probes)
    saveRDS(corrected, file.path(prepared, "bcr_erplus60.rds"), version = 3L)
  }, warning = function(w) { prep_warnings <<- c(prep_warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
  files <- c(dest, file.path(prepared, c("stock.rds", "bcr_historical.rds", "bcr_erplus60.rds")))
  manifest <- data.frame(path = substring(files, nchar(data_dir) + 2L),
                         sha256 = vapply(files, .application_hash, ""),
                         bytes = file.info(files)$size, redistribution = "external-cache: not included in public bundle")
  utils::write.csv(manifest, file.path(data_dir, "INPUT_MANIFEST.csv"), row.names = FALSE)
  utils::write.csv(as.data.frame(table(prep_warnings)), file.path(cache, "preparation_warnings.csv"), row.names = FALSE)
  comparison <- data.frame(probe = union(historical$panel_probes, corrected$panel_probes))
  comparison$historical_panel <- comparison$probe %in% historical$panel_probes
  comparison$wang_erplus60 <- comparison$probe %in% corrected$panel_probes
  comparison$historical_target <- comparison$probe %in% historical$target_probes
  comparison$corrected_target <- comparison$probe %in% corrected$target_probes
  utils::write.csv(comparison, file.path(cache, "bcr_panel_comparison.csv"), row.names = FALSE)
  invisible(manifest)
}

# Deliberately narrow warning policy for new BCR likelihood-ratio runs.
# A boundary fitted-probability warning alone is not classified as nonconvergence.
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
