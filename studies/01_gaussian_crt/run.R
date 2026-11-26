# Table 1 (size) and Figure 1 (power). Run from the repository root.
source('R/setup.R')
study <- '01_gaussian_crt'
kind <- argument('--experiment', 'size')
stopifnot(kind %in% c('size', 'power'))
source(file.path('studies', study, paste0(kind, '_settings.R')))
for (m in c('oracle_common', 'oracle_size', 'oracle_power', 'lw'))
  source(file.path('studies', study, 'R', paste0(m, '.R')))
out <- output_directory(study, kind)
save_settings(config, out)
all_rows <- list()
for (i in seq_len(nrow(settings))) {
  cfg <- settings[i, ]
  coefs <- if (kind == 'size') size_coefficients(cfg$s, cfg$t, i) else power_coefficients(cfg$s, cfg$t, i)
  theta_grid <- if (kind == 'size') 0 else config$theta
  for (j in seq_along(theta_grid)) {
    for (b in seq_len(config$repetitions)) {
      if (kind == 'size') {
        value <- capture_step(size_replication_lw(b, cfg, coefs,
          setting_index = i,
          seed_base = config$master_seed, M = config$M, lambda_ridge = config$lambda_ridge
        ))
        rows <- expand.grid(replication = b, statistic = c('raw score', 'poly2 score', 'ridge gain'), method = config$methods)
        rows$p_value <- if (value$status == 'ok') as.vector(value$value$pvalue) else NA_real_
        rows$observed_statistic <- if (value$status == 'ok') rep(value$value$statistic, length(config$methods)) else NA_real_
      } else {
        value <- capture_step(power_replication_lw(b, cfg, coefs,
          theta = theta_grid[j], theta_index = j,
          setting_index = i, seed_base = config$master_seed, M = config$M, lambda_ridge = config$lambda_ridge
        ))
        rows <- data.frame(replication = b, method = config$methods, theta = theta_grid[j])
        rows$p_value <- if (value$status == 'ok') unname(value$value$pvalue[config$methods]) else NA_real_
        rows$observed_statistic <- if (value$status == 'ok') value$value$statistic else NA_real_
      }
      for (name in names(cfg)) rows[[name]] <- cfg[[name]]
      rows$status <- value$status; rows$error <- value$error; rows$warnings <- value$warnings
      rows <- mark_lw_failures(rows, value)
      all_rows[[length(all_rows) + 1L]] <- rows
      if (b %% 100L == 0L) message(kind, ": ", cfg$setting, ", theta=", theta_grid[j], ", replicate ", b, " / ", config$repetitions)
    }
    write_report_csv(do.call(rbind, all_rows), file.path(out, 'results.csv'))
    message(kind, ': ', cfg$setting, ', theta=', theta_grid[j], ' completed')
  }
}

writeLines(capture.output(sessionInfo()), file.path(out, 'sessionInfo.txt'))
