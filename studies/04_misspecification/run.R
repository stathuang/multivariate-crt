# Figure S8 and Table S3: Gaussian, t5, and t3 covariate errors.
source('R/setup.R')
study <- '04_misspecification'
source(file.path('studies', study, 'settings.R'))
source(file.path('studies', study, 'functions.R'))
out <- output_directory(study, 'simulation')
save_settings(config, out)
RNGkind(config$rng_kind, config$normal_kind, config$sample_kind)
set.seed(config$root_seed)
seed <- .Random.seed
# The paper's formal stream starts after one check stream and 60 pilot streams.
# Advance those streams without rerunning the pilot experiments.
for (i in seq_len(1L + 12L * config$pilot_repetitions)) seed <- parallel::nextRNGStream(seed)
cells <- cells_for(config)
for (i in seq_len(nrow(cells))) {
  rows <- list()
  for (b in seq_len(config$repetitions)) {
    next_seed <- parallel::nextRNGStream(seed)
    if (!quick || b <= 2L) {
      answer <- one_replication(seed, cells$n[i], cells$distribution[i], config)
      rows[[length(rows) + 1L]] <- data.frame(
        cell_id = i, distribution = cells$distribution[i],
        n = cells$n[i], repetition = b, status = 'ok', T0 = answer$T0, p_value = answer$p_value,
        exceedances = answer$exceedances, ties = answer$ties, reject = answer$reject
      )
    }
    seed <- next_seed
  }
  write_report_csv(do.call(rbind, rows), file.path(out, sprintf('cell-%02d-%s-n%d.csv', i, cells$distribution[i], cells$n[i])))
  message('Completed ', cells$distribution[i], ', n=', cells$n[i])
}

writeLines(capture.output(sessionInfo()), file.path(out, 'sessionInfo.txt'))
