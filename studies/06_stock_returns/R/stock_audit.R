# Passive records for this isolated stock campaign. No RNG draws or resets here.
.stock_audit <- new.env(parent = emptyenv())
stock_audit_reset <- function(config) {
  rm(list = ls(.stock_audit, all.names = TRUE), envir = .stock_audit)
  .stock_audit$config <- config
  .stock_audit$current_method <- ""
  .stock_audit$ranks <- list(); .stock_audit$methods <- list()
  .stock_audit$samplers <- list(); .stock_audit$forests <- list()
  .stock_audit$coordinate_pvalues <- list()
  .stock_audit$R <- R.version.string; .stock_audit$RNGkind <- RNGkind()
  invisible(NULL)
}
stock_audit_export <- function() as.list(.stock_audit, all.names = TRUE)

stock_record_rank <- function(observed, randomized, rule, greater, equal, offset,
                              pvalue, rng_before, rng_after) {
  entry <- list(method = .stock_audit$current_method, observed = observed,
    randomized = randomized, M = length(randomized), rule = rule,
    greater = greater, equal_including_observed = equal, offset = offset,
    pvalue = pvalue, rng_before = rng_before, rng_after = rng_after)
  .stock_audit$ranks[[length(.stock_audit$ranks) + 1L]] <- entry
  invisible(NULL)
}

stock_randomForest <- function(...) {
  fit <- randomForest::randomForest(...)
  method <- .stock_audit$current_method
  expected <- .stock_audit$config$forest_trees
  if (!is.null(expected) && fit$ntree != expected) stop("Unexpected actual forest size")
  .stock_audit$forests[[method]] <- c(.stock_audit$forests[[method]], fit$ntree)
  fit
}

.stock_original_sampler <- mscrt_sample
mscrt_sample <- function(X, graph, target = seq_len(ncol(X)), M = 100L, L = 3L,
                         reduced = TRUE) {
  before <- .Random.seed
  answer <- .stock_original_sampler(X, graph, target, M, L, reduced)
  after <- .Random.seed
  .stock_audit$samplers[[length(.stock_audit$samplers) + 1L]] <- list(
    X_hash = sha256_object(X), graph_hash = sha256_object(graph), target = target,
    M = M, L = L, reduced = reduced, output_hash = sha256_object(answer),
    rng_before = before, rng_after = after)
  if (identical(.stock_audit$config$design, "observed"))
    .stock_audit$observed_copies <- answer
  answer
}

stock_method_begin <- function(method, X, y, target) {
  .stock_audit$current_method <- method
  .stock_audit$methods[[method]] <- list(rng_before = .Random.seed,
    X_hash = sha256_object(X), y_hash = sha256_object(y), target = target)
  invisible(NULL)
}
stock_method_end <- function(method) {
  .stock_audit$methods[[method]]$rng_after <- .Random.seed
  .stock_audit$current_method <- ""
  invisible(NULL)
}
.stock_original_gcrt <- gcrt_test
gcrt_test <- function(X, y, target, copies, statistic,
                     response_family = "gaussian", p_rule = "discrete") {
  stock_method_begin(statistic, X, y, target)
  on.exit(stock_method_end(statistic))
  .stock_audit$methods[[statistic]]$copies_hash <- sha256_object(copies)
  .stock_original_gcrt(X, y, target, copies, statistic, response_family, p_rule)
}
.stock_original_baseline <- baseline_test_detail
baseline_test_detail <- function(X, y, target, method, response_family = "gaussian") {
  stock_method_begin(method, X, y, target)
  on.exit(stock_method_end(method))
  .stock_original_baseline(X, y, target, method, response_family)
}
.stock_original_bonferroni <- .baseline_bonferroni
.baseline_bonferroni <- function(pvalues) {
  result <- .stock_original_bonferroni(pvalues)
  .stock_audit$coordinate_pvalues[[.stock_audit$current_method]] <- list(
    individual = pvalues, combined = result)
  result
}
.stock_original_application <- application_replicate
application_replicate <- function(config, replicate, data_dir) {
  if (!identical(config$dataset, "stock")) stop("This campaign permits stocks only")
  stock_audit_reset(config)
  .stock_original_application(config, replicate, data_dir)
}
