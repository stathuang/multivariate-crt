# Figures S1-S7: select one panel, keeping all methods and signal levels.
source('R/setup.R')
study <- '02_graphical_crt'
panel <- argument('--figure', 'S1a')
settings <- file.path('studies', study, 'settings', paste0('Figure', panel, '.R'))
if (!file.exists(settings)) stop('Unknown figure panel: ', panel)
source(settings)
if (quick) config$repetitions <- 2L
load_modules(c(
  'sampler', 'data', 'linear_tests', 'gaussian_tests', 'dcrt', 'hdi_cache',
  'cgm_diagnostics', 'comparators', 'experiment', 'binary_data', 'binary_likelihood',
  'result_diagnostics', 'binary_tests', 'binary_experiment', 'extended_designs', 'supergraph_experiment'
))
config <- prepare_gaussian_config(config)
out <- output_directory(study, paste0('Figure', panel))
save_settings(config, out)
rows <- list()
for (replicate in seq_len(config$repetitions)) {
  set.seed(config$seed_base + replicate)
  if (config$family == 'supergraph') {
    design <- capture_step(supergraph_design(config))
    # Finish sampling and testing for each graph before moving to the next graph.
    answer <- evaluate_supergraphs(design, config, replicate)
  } else {
    design <- capture_step(if (config$family == 't3')
      t3_design(
        config$n, config$p, config$bandwidth, config$t3_lag_coefficients,
        config$t3_df, config$t3_noise_scale
      ) else
      gaussian_design(config$n, config$p, config$bandwidth, config$precision_offdiag))
    cache <- list(enabled = FALSE)
    if (design$status == 'ok') {
      if (config$response_family == 'gaussian') cache <- prepare_de_lasso(config, replicate, design$value$X)
      sampling <- capture_step(mscrt_sample(design$value$X, design$value$graph, config$target,
        M = config$M, L = config$L, reduced = TRUE
      ))
    } else sampling <- design
    answer <- if (config$response_family == 'binomial') evaluate_binary(design, sampling, config, replicate) else
      evaluate_gaussian(design, sampling, cache, config, replicate)
  }
  rows[[replicate]] <- answer$results
  write_report_csv(do.call(rbind, rows), file.path(out, 'results.csv'))
  message('Figure ', panel, ': replicate ', replicate, ' / ', config$repetitions)
}

writeLines(capture.output(sessionInfo()), file.path(out, 'sessionInfo.txt'))
