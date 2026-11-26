# Figure 2, Figure S11, and Tables S9-S11.
source('R/setup.R')
study <- '06_stock_returns'
part <- argument('--experiment', 'observed-T1')
configs <- jsonlite::fromJSON(file.path('studies', study, 'settings.json'), simplifyVector = FALSE)
ids <- c(
  'observed-T1' = 'stock_three_stat_replay_T1', 'observed-T2' = 'stock_three_stat_replay_T2',
  'null-T1' = 'stock_typeI_T1', 'null-T2' = 'stock_typeI_T2', 'noise-05' = 'stock_noise_05',
  'noise-10' = 'stock_noise_10', 'noise-15' = 'stock_noise_15', 'noise-20' = 'stock_noise_20', 'noise-25' = 'stock_noise_25'
)
out <- output_directory(study, part)
if (part %in% c('graph', 'gof')) {
  load_modules(c('sampler', 'application_diagnostics'))
  stock <- readRDS('data/processed/stock.rds')
  if (part == 'graph') {
    set.seed(1L)
    plot_stock_graph(stock, file.path(out, 'FigureS11.pdf'))
    saveRDS(list(graph = stock$graph, stocks = stock$stocks, sectors = stock$sectors), file.path(out, 'graph.rds'))
  } else {
    set.seed(2023L); M <- 1000L; L <- 3L
    vv <- stock_VV(stock$X, stock$graph); dp <- stock_DP(stock$X, stock$graph)
    message('Generating ', M, ' copies for the stock graph goodness-of-fit test...')
    copies <- mscrt_sample(stock$X, stock$graph, seq_len(ncol(stock$X)), M, L, reduced = FALSE)
    observed <- stock_F_sum(stock$X, stock$graph)
    randomized <- numeric(M)
    for (m in seq_len(M)) {
      randomized[m] <- stock_F_sum(copies$X_copies[[m]], stock$graph)
      if (m %% 100L == 0L) message('Scored ', m, ' / ', M, ' copies')
    }
    pv <- c(VV = vv, DP = dp, F_sum = rank_pvalue(observed, randomized, rule = 'discrete'))
    write_report_csv(data.frame(method = names(pv), pvalue = unname(pv)), file.path(out, 'results.csv'))
    saveRDS(list(pvalues = pv, observed = observed, randomized = randomized, M = M, L = L), file.path(out, 'gof.rds'))
  }
} else {
  if (!part %in% names(ids)) stop('Unknown stock experiment: ', part)
  cfg <- Filter(function(c) c$id == ids[[part]], configs)[[1L]]
  for (k in c('statistics', 'baselines', 'methods', 'graph_ids')) cfg[[k]] <- unlist(cfg[[k]])
  # Initialize numerical packages before assigning replicate seeds.
  for (p in c(
    'jsonlite', 'digest', 'MASS', 'Matrix', 'mvtnorm', 'glmnet', 'randomForest', 'glasso',
    'bestNormalize', 'VGAM', 'hdi', 'SIHR', 'rags2ridges', 'lpSolve', 'gmp', 'Rmpfr', 'CVXR', 'ECOSolveR'
  ))
    stopifnot(requireNamespace(p, quietly = TRUE))
  e <- new.env(parent = globalenv())
  load_modules(c('binary_likelihood', 'cgm_diagnostics', 'result_diagnostics', 'multinomial_likelihood'), e)
  e$sha256_object <- function(x) digest::digest(x, algo = 'sha256', serialize = TRUE)
  for (m in c('core', 'statistics', 'baselines', 'applications', 'stock_audit'))
    sys.source(file.path('studies', study, 'R', paste0(m, '.R')), e)
  save_settings(cfg, out)
  reps <- if (quick) min(2L, cfg$repetitions) else cfg$repetitions
  rows <- list()
  for (i in seq_len(reps)) {
    set.seed(if (is.null(cfg$observed_seed)) cfg$seed_base + i else cfg$observed_seed)
    rows[[i]] <- e$application_replicate(cfg, i, 'data/processed')
    write_report_csv(do.call(rbind, rows), file.path(out, 'results.csv'))
    message(part, ': ', i, ' / ', reps)
  }
}

writeLines(capture.output(sessionInfo()), file.path(out, 'sessionInfo.txt'))
