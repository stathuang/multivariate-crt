# Table S4 and Figures S9-S10, including the corresponding rate/paired tables.
source('R/setup.R')
study <- '05_graphical_power'
kind <- argument('--experiment', 'direction')
stopifnot(kind %in% c('direction', 'independent', 'connected'))
cfg <- jsonlite::fromJSON(file.path('studies', study, paste0(kind, '_settings.json')))
RNGkind(cfg$rng[1], cfg$rng[2], cfg$rng[3])
out <- output_directory(study, kind)
save_settings(cfg, out)
M <- as.integer(cfg$M); N <- M + 1L
cases <- if (kind == 'independent') names(cfg$cases) else cfg$cases
methods <- cfg$methods
source(file.path('studies', study, paste0(kind, '.R')))
rows <- list()
if (kind %in% c('direction', 'independent')) {
  reps <- if (quick) 2L else if (kind == 'direction') cfg$direction_repetitions_per_n else cfg$repetitions_per_n
  for (k in seq_along(cfg$sample_sizes)) {
    n <- as.integer(cfg$sample_sizes[k])
    for (b in seq_len(reps)) {
      seed <- as.integer(cfg$seed_base + 10000 * k + b)
      id <- sprintf('n%d-r%04d', n, b)
      r <- if (kind == 'direction') one(n, b, seed, id) else one(n, b, seed)
      if (kind == 'direction') {
        pv <- apply(r$original, 2, rank_p)
        row <- data.frame(
          dataset_id = id, n = n, replicate = b, seed = seed, case = cases[1:3],
          pvalue = pv, reject = pv <= cfg$alpha, T0 = r$original[1, ], Tbar = colMeans(r$original[-1, , drop = FALSE]),
          gap = r$original[1, ] - colMeans(r$original[-1, , drop = FALSE])
        )
      } else row <- do.call(rbind, lapply(methods, function(me) data.frame(
        dataset_id = id, n = n, replicate = b, seed = seed, case = cases, method = me, pvalue = r$pvalue[, me],
        reject = r$reject[, me], T0 = r$scores[1, , me], Tbar = colMeans(r$scores[-1, , me]),
        gap = r$scores[1, , me] - colMeans(r$scores[-1, , me])
      )))
      rows[[length(rows) + 1L]] <- row
    }
    write_report_csv(do.call(rbind, rows), file.path(out, 'tests.csv'))
    message(kind, ': n=', n, ' complete')
  }
} else {
  reps <- if (quick) 2L else cfg$replicates
  for (k in seq_len(reps)) {
    a <- one_replicate(k)
    base <- test_rows(a)
    # Two uncentred controls use the very same copies and quadratic weights.
    controls <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
      C <- a$c[, , i]
      scores <- list(uncentered_pool = colSums(C * solve(a$Dpool, C)), uncentered_identity = colSums(C^2))
      do.call(rbind, lapply(names(scores), function(me) data.frame(
        case = cells$case[i], method = me,
        replicate = k, data_seed = a$data_seed, sampler_seed = a$sampler_seed, b = cells$b[i], r = cells$r[i],
        pvalue = rank_p(scores[[me]]), reject = rank_p(scores[[me]]) <= cfg$alpha
      )))
    }))
    rows[[k]] <- rbind(base, controls)
    if (k %% 25L == 0L || k == reps) {
      write_report_csv(do.call(rbind, rows), file.path(out, 'tests.csv'))
      message('Connected targets: ', k, ' / ', reps)
    }
  }
}

writeLines(capture.output(sessionInfo()), file.path(out, 'sessionInfo.txt'))
