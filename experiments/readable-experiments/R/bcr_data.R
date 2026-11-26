# Adopted ER-positive 60-probe panel; input files are read only.
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

prepare_bcr <- function(input_dir, panel = "wang_erplus60",
                        selection_seed = 2023L, erplus_probes = NULL, null_seed = 2023L) {
  if (!identical(panel, "wang_erplus60")) stop("This entry uses the adopted ER-positive 60-probe panel")
  clinical <- utils::read.csv(file.path(input_dir, "GSE2034-PI.csv"), fileEncoding = "UTF-8-BOM")
  expression <- read_geo_expression(file.path(input_dir, "GSE2034_series_matrix.txt"))
  listed <- utils::read.csv(file.path(input_dir, "ImpVar.csv"), header = FALSE,
                           fileEncoding = "UTF-8-BOM")[[1L]]
  if (length(listed) != 76L || anyDuplicated(listed)) stop("Expected 76 unique listed probes")
  ordered <- expression$probes[expression$probes %in% listed]
  if (length(ordered) != 76L) stop("A listed probe is absent")
  if (is.null(erplus_probes)) erplus_probes <- listed[1:60]
  probes <- as.character(erplus_probes)
  if (length(probes) != 60L || anyNA(probes) || anyDuplicated(probes) ||
      !identical(probes, as.character(listed[1:60])))
    stop("Corrected panel requires the verified Wang Table 3 ER-positive probe IDs in exact order")
  # Use this validated list throughout; match preserves its order.
  probe_index <- match(probes, expression$probes)
  if (anyNA(probe_index) || anyDuplicated(probe_index)) stop("Panel probes must be unique and present")
  join <- match(clinical$GEO.asscession.number, expression$samples)
  if (anyNA(join) || anyDuplicated(join)) stop("Clinical/expression join is not one to one")
  # Original orderNorm is fitted on ALL 286 patients before restriction to ER+.
  raw <- expression$X[join, probe_index, drop = FALSE]
  if (!identical(colnames(raw), as.character(probes))) stop("Panel probe order mismatch")
  if (any(!is.finite(raw))) stop("Non-finite selected expression data")
  normalized <- apply(raw, 2L, function(v) bestNormalize::orderNorm(v)$x.t)
  erplus <- clinical$ER.Status == "ER+"
  if (anyNA(erplus) || sum(erplus) != 209L) stop("Expected 209 ER-positive patients")
  y <- bcr_endpoint(clinical)[erplus]
  X <- normalized[erplus, , drop = FALSE]
  # Select the larger panel at lambda.min and background at lambda.1se.
  selection <- with_seed(selection_seed, {
    cv <- glmnet::cv.glmnet(X, y, family = "multinomial", nfolds = 5L)
    fit <- glmnet::glmnet(X, y, family = "multinomial")
    selected <- function(s) {
      coefficients_by_class <- stats::coef(fit, s = s)
      coefficient_matrix <- as.matrix(do.call(cbind, coefficients_by_class))
      squared_coefficients <- coefficient_matrix^2
      total_squared <- rowSums(squared_coefficients)
      # The first row is the intercept; the remaining indices refer to probes.
      probe_squared <- total_squared[-1L]
      which(probe_squared != 0)
    }
    list(strong = selected(cv$lambda.min), condition = selected(cv$lambda.1se),
         lambda_min = cv$lambda.min, lambda_1se = cv$lambda.1se, foldid = cv$foldid)
  })
  target_panel <- setdiff(selection$strong, selection$condition)
  if (!length(target_panel)) stop("Selected target is empty")
  if (!all(selection$condition %in% selection$strong))
    stop("One-SE selected genes are not nested; hypothesis requires explicit resolution")
  mat_X <- unname(X[, selection$strong, drop = FALSE])
  target <- match(target_panel, selection$strong)
  graph <- matrix(1L, ncol(mat_X), ncol(mat_X))
  diag(graph) <- 0L
  outside <- setdiff(seq_len(ncol(mat_X)), target)
  # The null response fit uses only the selected background probes.
  null <- with_seed(null_seed, {
    cv <- glmnet::cv.glmnet(mat_X[, outside, drop = FALSE], y, family = "multinomial", nfolds = 5L)
    fit <- glmnet::glmnet(mat_X[, outside, drop = FALSE], y, family = "multinomial")
    list(fit = fit, lambda = cv$lambda.min, seed = null_seed,
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

# One unrestricted fitted-null dataset. Declared classes may have zero observations.
# Local RF fitting drops unused levels; the response here keeps its original codes.
bcr_null_data <- function(bcr) {
  null <- bcr$null
  X <- rags2ridges::rmvnormal(nrow(bcr$X), mu = null$mean, sigma = null$covariance)
  probabilities <- stats::predict(null$fit, newx = X[, null$outside, drop = FALSE],
                                  type = "response", s = null$lambda)[, , 1L]
  draw <- apply(probabilities, 1L, function(v)
    which(stats::rmultinom(1L, 1L, v)[, 1L] == 1L))
  Y <- factor(draw, levels = seq_len(ncol(probabilities)))
  list(X = X, Y = Y, probabilities = probabilities, class_counts = table(Y))
}
