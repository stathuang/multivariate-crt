# Rebuild the saved stock graph and fixed null response fits from supplied inputs.
# The transformed X and aligned response are the reproducibility starting point.
# Stock downloads and raw-price alignment are outside this preparation step.
prepare_stock <- function(input_dir, null_seed = 2023L) {
  e <- new.env(parent = emptyenv())
  load(file.path(input_dir, "DataForExp.RData"), envir = e)
  pre <- new.env(parent = emptyenv())
  load(file.path(input_dir, "Pre-X.RData"), envir = pre)
  names <- utils::read.csv(file.path(input_dir, "list_stocks.csv"))$Symbol
  X <- unname(e$mat_X)
  y <- as.numeric(e$Y)
  graph <- unname(e$G1 * 1L)
  if (!identical(dim(X), c(92L, 103L)) || length(y) != nrow(X) ||
      any(!is.finite(X)) || any(!is.finite(y)) || length(names) != ncol(X)) stop("Unexpected saved stock input")
  if (max(abs(X - pre$X)) != 0) stop("Saved stock X disagrees with preprocessing snapshot")
  within_sector <- outer(as.character(pre$industry), as.character(pre$industry), `==`) * 1L
  penalty <- (1 - within_sector) * 0.15
  diag(penalty) <- 0
  sample_covariance <- stats::var(X)
  glasso_fit <- glasso::glasso(sample_covariance, rho = penalty)
  precision_matrix <- glasso_fit$wi
  fitted_graph <- (precision_matrix != 0) * 1L
  diag(fitted_graph) <- 0L
  if (!identical(dim(fitted_graph), dim(graph)) || any(fitted_graph != graph))
    stop("Reconstructed graph differs from saved graph")
  targets <- list(T1 = as.integer(e$setT1), T2 = as.integer(e$setT2))
  nulls <- lapply(targets, function(target) with_seed(null_seed, {
    outside <- setdiff(seq_len(ncol(X)), target)
    fit <- glmnet::cv.glmnet(X[, outside, drop = FALSE], y)
    fitted_y <- stats::predict(fit, newx = X[, outside, drop = FALSE], s = fit$lambda.min)
    coefficients <- as.numeric(stats::coef(fit, s = fit$lambda.min))
    residual_sd <- stats::sd(y - fitted_y)
    list(coef = coefficients, sd = residual_sd, lambda = fit$lambda.min,
         outside = outside, seed = null_seed)
  }))
  list(X = X, y = y, graph = graph, targets = targets, nulls = nulls, stocks = names,
       target_names = lapply(targets, function(t) names[t]), sectors = as.character(pre$industry),
       start_date = as.character(pre$start_date), end_date = as.character(pre$end_date),
       starting_point = "Saved transformed covariates and saved aligned DIA response; raw alignment not reconstructed",
       null_seed_status = "New explicit preparation seed; historical stock null fitting seed was not recorded")
}

# E21: each repetition starts at the same observed matrix. These 20 sweeps
# generate the covariates; config$L separately controls the test's copies.
stock_null_data <- function(stock, group, generation_sweeps = 20L) {
  X <- mscrt_sample(stock$X, stock$graph, target = seq_len(ncol(stock$X)),
                    M = 1L, L = generation_sweeps, reduced = FALSE)$X_copies[[1L]]
  null <- stock$nulls[[group]]
  Y <- null$coef[1L] + as.numeric(X[, null$outside, drop = FALSE] %*% null$coef[-1L]) +
    stats::rnorm(nrow(X), sd = null$sd)
  list(X = X, Y = Y)
}

# E22: add fresh Gaussian noise to the saved response at the specified scale.
stock_noise_data <- function(stock, sigma) {
  list(X = stock$X, Y = stock$y + stats::rnorm(length(stock$y), sd = sigma))
}
