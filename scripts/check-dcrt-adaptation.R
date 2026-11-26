# Compare the adopted restricted implementations with the historical dCRT source.
# Parse only the two required function definitions; never source the full legacy file.
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
script <- sub('^--file=', '', grep('^--file=', commandArgs(), value=TRUE)[1L])
root <- normalizePath(file.path(dirname(script), '..'))
legacy <- new.env(parent=globalenv())
suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(randomForest))
for (expr in parse(args[1L])) {
  if (is.call(expr) && as.character(expr[[1L]]) %in% c('<-', '=') &&
      is.symbol(expr[[2L]]) && as.character(expr[[2L]]) %in% c('d_CRT', 'Creat_condition_gaussian'))
    eval(expr, legacy)
}
stopifnot(is.function(legacy$d_CRT), is.function(legacy$Creat_condition_gaussian))
readable <- new.env(parent=globalenv())
source(file.path(root, 'experiments/readable-experiments/R/sampler.R'), local=readable)
source(file.path(root, 'experiments/readable-experiments/R/run_tools.R'), local=readable)
source(file.path(root, 'experiments/readable-experiments/R/comparators.R'), local=readable)
source(file.path(root, 'experiments/readable-experiments/R/linear_tests.R'), local=readable)
source(file.path(root, 'experiments/readable-experiments/R/dcrt.R'), local=readable)
bcr <- new.env(parent=globalenv())
source(file.path(root, 'experiments/bcr-original/R/baselines.R'), local=bcr)
stock <- new.env(parent=globalenv())
source(file.path(root, 'experiments/stock-rf1000/engine/R/baselines.R'), local=stock)
# The campaign's actual stock_randomForest delegates to randomForest and records
# tree counts without drawing or resetting random numbers.
stock$stock_randomForest <- randomForest::randomForest
RNGkind('Mersenne-Twister', 'Inversion', 'Rejection')
set.seed(2034)
X <- matrix(rnorm(120 * 8), 120, 8)
continuous <- 0.3 * X[, 1] - 0.2 * X[, 3] + rnorm(120)
binary <- rbinom(120, 1, plogis(continuous))
target <- c(3L, 1L)
rows <- list()
for (implementation in c('readable', 'bcr', 'stock')) {
  for (model in c('gaussian', 'binomial', 'RF')) {
    y <- if (model == 'binomial') binary else continuous
    trees <- if (implementation == 'stock') 1000L else 100L
    set.seed(2023)
    ref <- legacy$d_CRT(X, y, candidate_set=target, CDR='No', k=0,
      model=switch(model, gaussian='Gaussian_lasso', binomial='Binomial_lasso', RF='RF'),
      RF.num.trees=c(trees, 30L), MC_free=TRUE, return.pvalues=TRUE)$p_values
    rng_ref <- .Random.seed
    set.seed(2023)
    got <- if (implementation == 'readable') readable$dcrt_test(X, y, target, model)$pvalue else
      get(implementation)$.baseline_dcrt(X, y, target, model)
    difference <- abs(got - min(1, length(target) * min(ref)))
    same_rng <- identical(.Random.seed, rng_ref)
    stopifnot(difference < 1e-12, same_rng)
    rows[[length(rows) + 1L]] <- data.frame(implementation, model, trees,
      absolute_pvalue_difference=difference, identical_final_rng=same_rng)
  }
}
out <- file.path(root, 'outputs/dcrt-adaptation')
dir.create(out, recursive=TRUE, showWarnings=FALSE)
write.csv(do.call(rbind, rows), file.path(out, 'parity.csv'), row.names=FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, 'session.txt'))
print(do.call(rbind, rows), row.names=FALSE)
