# Figure S12, Table S12, and observed p-values in Section 8.2.
source('R/setup.R')
study <- '07_breast_cancer'
kind <- argument('--experiment', 'observed')
stopifnot(kind %in% c('observed', 'null'))
observed <- kind == 'observed'
out <- output_directory(study, kind)
old <- new.env(parent = globalenv()); current <- new.env(parent = globalenv())
for (p in c(
  'jsonlite', 'digest', 'MASS', 'Matrix', 'glmnet', 'randomForest', 'rags2ridges',
  'VGAM', 'lpSolve', 'gmp', 'Rmpfr', 'CVXR', 'ECOSolveR', 'nnet'
))
  stopifnot(requireNamespace(p, quietly = TRUE))
sys.source(file.path(root, 'R/sampler.R'), old)
for (f in c('statistics.R', 'baselines.R', 'applications.R'))
  sys.source(file.path(root, 'studies', study, 'R', f), old)
for (f in c('sampler.R', 'data.R', 'linear_tests.R', 'gaussian_tests.R', 'multinomial_likelihood.R', 'bcr_tests.R'))
  sys.source(file.path(root, 'R', f), current)

cfg <- jsonlite::fromJSON(file.path('studies', study, paste0(kind, '_settings.json')), simplifyVector = FALSE)
for (k in c('statistics', 'baselines', 'methods', 'graph_ids')) cfg[[k]] <- unlist(cfg[[k]])
data_dir <- file.path(root, 'data/processed')
save_settings(cfg, out)
raw_gcrt <- old$gcrt_test
saved <- NULL
old$gcrt_test <- function(...) {
  a <- list(...)
  saved <<- list(arguments = a, rng = .Random.seed)
  do.call(raw_gcrt, a)
}
ids <- if (observed) 1L else if (quick) c(1L, 74L) else 1:400
results <- list()
RNGkind('Mersenne-Twister', 'Inversion', 'Rejection')
for (i in ids) {
  set.seed(if (observed) cfg$observed_seed else cfg$seed_base + i)
  legacy <- old$application_replicate(cfg, i, data_dir)
  stopifnot(!is.null(saved))
  # Legacy calls above consume exactly their original random stream, including old failures.
  a <- saved$arguments
  fit <- current$bcr_lrt_test(a[[1L]], a[[2L]], a[[3L]])
  repaired <- legacy
  k <- which(repaired$method == 'LRT')
  repaired$pvalue[k] <- fit$pvalue;repaired$status[k] <- fit$diagnostics$objective_status;repaired$error[k] <- ''
  k <- which(repaired$method == 'RF')
  if (repaired$status[k] != 'ok') {
    assign('.Random.seed', saved$rng, envir = .GlobalEnv)
    fixed <- current$gcrt_rf_test(a[[1L]], a[[2L]], a[[3L]], a[[4L]], 'RF', cfg$p_rule, response_family = 'multinomial')
    repaired$pvalue[k] <- fixed$pvalue;repaired$status[k] <- 'ok';repaired$error[k] <- ''
  }
  stopifnot(all(repaired$status %in% c('ok', 'numerically_certified')), all(is.finite(repaired$pvalue)))
  saveRDS(
    list(results = repaired, lrt_diagnostics = fit$diagnostics),
    file.path(out, sprintf('replicate-%05d.rds', i))
  )
  message('Breast cancer: replicate ', i)
  results[[as.character(i)]] <- repaired
}
write.csv(do.call(rbind, results), file.path(out, 'results.csv'), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, 'sessionInfo.txt'))
