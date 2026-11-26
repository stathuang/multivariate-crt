# Deterministic analyses of saved c/D components. No data or copy generation.
options(warn = 2, digits = 17)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
analysis_dir <- normalizePath(args[[1]], mustWork = TRUE)
manifest <- jsonlite::fromJSON(file.path(analysis_dir, 'INPUTS.json'))
source_dir <- manifest$source_directory
out_dir <- file.path(analysis_dir, 'results')
stopifnot(dir.exists(out_dir), !file.exists(file.path(out_dir, 'CHECK.json')))
cfg <- jsonlite::fromJSON(file.path(source_dir, 'run/CONFIG.json'))
stopifnot(cfg$replicates == 2500L, cfg$t == 8L, cfg$M == 99L,
          cfg$alpha == 0.05)
rng_exists <- exists('.Random.seed', envir = .GlobalEnv, inherits = FALSE)
rng_before <- if (rng_exists) get('.Random.seed', envir = .GlobalEnv) else NULL
methods <- c('original', 'uncentered_pool', 'centered_pool',
             'uncentered_identity', 'centered_identity')
original_methods <- c('original', 'centered_pool', 'centered_identity')
cells <- expand.grid(r = cfg$signal, b = cfg$background_sign)
cells$case <- sprintf('b%s_r%02d', ifelse(cells$b < 0, 'neg', 'pos'),
                      round(100 * cells$r))
rank_p <- function(x) (1 + sum(x[-1L] >= x[1L])) / length(x)
files <- sort(list.files(file.path(source_dir, 'run/records'),
                         pattern = '^batch-[0-9]{2}[.]rds$', full.names = TRUE))
stopifnot(length(files) == 25L)
tests <- vector('list', cfg$replicates)
seen <- integer()
pool_error <- original_score_error <- new_score_error <- 0
record_rank_mismatch <- independent_rank_mismatch <- 0L
minimum_new_rank_spacing <- Inf

for (file in files) {
  records <- readRDS(file)
  for (a in records) {
    k <- a$replicate
    stopifnot(k %in% seq_len(cfg$replicates), !(k %in% seen),
              identical(dim(a$c), c(8L, 100L, 8L)),
              identical(dim(a$D), c(8L, 8L, 100L)),
              identical(dimnames(a$c)[[3L]], cells$case),
              a$data_seed == cfg$data_seed_base + k,
              a$sampler_seed == cfg$sampler_seed_base + k,
              all(is.finite(a$c)), all(is.finite(a$D)))
    seen <- c(seen, k)
    K <- cfg$M + 1L
    pool <- apply(a$D, c(1L, 2L), mean)
    pool_error <- max(pool_error, max(abs(pool - a$Dpool)))
    stopifnot(min(eigen(pool, symmetric = TRUE, only.values = TRUE)$values) > 0)
    C <- matrix(a$c, nrow = cfg$t)
    delta <- C
    for (ca in seq_len(nrow(cells))) {
      idx <- ((ca - 1L) * K + 1L):(ca * K)
      delta[, idx] <- sweep(C[, idx], 1L, rowMeans(C[, idx]), '-')
    }
    score <- array(NA_real_, c(K, nrow(cells), length(methods)),
                   dimnames = list(0:cfg$M, cells$case, methods))
    score[, , 'uncentered_pool'] <- matrix(colSums(C * solve(pool, C)), K)
    score[, , 'centered_pool'] <- matrix(colSums(delta * solve(pool, delta)), K)
    score[, , 'uncentered_identity'] <- matrix(colSums(C^2), K)
    score[, , 'centered_identity'] <- matrix(colSums(delta^2), K)
    for (j in seq_len(K)) {
      z <- a$c[, j, ]
      score[j, , 'original'] <- colSums(z * solve(a$D[, , j], z))
    }
    original_score_error <- max(original_score_error,
      max(abs(score[, , original_methods] - a$scores[, , original_methods])))
    pv <- apply(score, c(2L, 3L), rank_p)
    record_rank_mismatch <- record_rank_mismatch +
      sum(pv[, original_methods] != a$pvalue[, original_methods])
    # Separate numerical forms for both added controls: Cholesky and crossprod.
    alt_pool <- matrix(colSums(forwardsolve(t(chol(pool)), C)^2), K)
    alt_identity <- matrix(vapply(seq_len(ncol(C)), function(j)
      drop(crossprod(C[, j])), 0.0), K)
    new_score_error <- max(new_score_error,
      max(abs(alt_pool - score[, , 'uncentered_pool'])),
      max(abs(alt_identity - score[, , 'uncentered_identity'])))
    independent_rank_mismatch <- independent_rank_mismatch +
      sum(apply(alt_pool, 2L, rank_p) != pv[, 'uncentered_pool']) +
      sum(apply(alt_identity, 2L, rank_p) != pv[, 'uncentered_identity'])
    for (me in c('uncentered_pool', 'uncentered_identity')) {
      for (ca in seq_len(nrow(cells))) {
        minimum_new_rank_spacing <- min(minimum_new_rank_spacing,
          min(abs(score[-1L, ca, me] - score[1L, ca, me])))
      }
    }
    stopifnot(all(is.finite(score)), all(score >= -1e-12))
    z <- expand.grid(case = cells$case, method = methods, stringsAsFactors = FALSE)
    z$replicate <- k
    z$data_seed <- a$data_seed
    z$sampler_seed <- a$sampler_seed
    z$b <- cells$b[match(z$case, cells$case)]
    z$r <- cells$r[match(z$case, cells$case)]
    z$pvalue <- as.vector(pv)
    z$reject <- z$pvalue <= cfg$alpha
    tests[[k]] <- z
  }
  cat(sprintf('Read and analysed saved batch %s; %d replications.\n',
              basename(file), length(seen)))
  flush.console()
}
tests <- do.call(rbind, tests)
row.names(tests) <- NULL
key <- c('replicate', 'case', 'method')
stopifnot(identical(sort(seen), seq_len(cfg$replicates)),
          nrow(tests) == 100000L, !anyDuplicated(tests[key]),
          all(table(tests$case, tests$method) == cfg$replicates))
old_tests <- read.csv(file.path(source_dir, 'run/tests.csv'))
replayed <- merge(tests[tests$method %in% original_methods, ], old_tests,
                  by = key, all = TRUE, suffixes = c('_calc', '_saved'))
stopifnot(nrow(replayed) == 60000L, !anyNA(replayed),
          all(replayed$pvalue_calc == replayed$pvalue_saved),
          all(replayed$reject_calc == replayed$reject_saved),
          all(replayed$data_seed_calc == replayed$data_seed_saved),
          all(replayed$sampler_seed_calc == replayed$sampler_seed_saved))

get_reject <- function(me, b, r) {
  z <- tests[tests$method == me & tests$b == b & tests$r == r, ]
  z <- z[order(z$replicate), ]
  stopifnot(identical(z$replicate, seq_len(cfg$replicates)))
  as.numeric(z$reject)
}
paired_stats <- function(d) {
  B <- length(d)
  value <- mean(d)
  se <- sd(d) / sqrt(B)
  count_se <- sqrt((sum(d^2) - B * value^2) / ((B - 1) * B))
  stopifnot(abs(se - count_se) < 1e-14)
  data.frame(replications = B, difference = value, paired_mcse = se,
             sum_difference = sum(d), sum_squared_difference = sum(d^2))
}
summary <- do.call(rbind, lapply(seq_len(nrow(cells)), function(ca)
  do.call(rbind, lapply(methods, function(me) {
    y <- get_reject(me, cells$b[ca], cells$r[ca])
    rate <- mean(y)
    data.frame(case = cells$case[ca], b = cells$b[ca], r = cells$r[ca],
               method = me, replications = length(y), rejects = sum(y),
               rate = rate, mcse = sqrt(rate * (1 - rate) / length(y)))
  }))))
old_summary <- read.csv(file.path(source_dir, 'run/summary.csv'))
summary_join <- merge(summary[summary$method %in% original_methods, ],
                      old_summary, by = c('case', 'b', 'r', 'method'),
                      suffixes = c('_calc', '_saved'))
stopifnot(nrow(summary_join) == 24L,
          all(summary_join$rejects_calc == summary_join$rejects_saved),
          max(abs(summary_join$rate_calc - summary_join$rate_saved)) < 1e-14)

weights <- list(pooled_Gram = c('centered_pool', 'uncentered_pool'),
                unit_weight = c('centered_identity', 'uncentered_identity'))
fixed <- do.call(rbind, lapply(seq_len(nrow(cells)), function(ca)
  do.call(rbind, lapply(names(weights), function(w) {
    pair <- weights[[w]]
    d <- get_reject(pair[1L], cells$b[ca], cells$r[ca]) -
         get_reject(pair[2L], cells$b[ca], cells$r[ca])
    cbind(data.frame(case = cells$case[ca], b = cells$b[ca], r = cells$r[ca],
                     weight = w, method_a = pair[1L], method_b = pair[2L]),
          paired_stats(d), centered_only = sum(d == 1), uncentered_only = sum(d == -1))
  }))))
signs <- do.call(rbind, lapply(cfg$signal, function(r)
  do.call(rbind, lapply(methods, function(me) {
    d <- get_reject(me, 1, r) - get_reject(me, -1, r)
    cbind(data.frame(r = r, method = me, b_a = 1, b_b = -1), paired_stats(d),
          positive_only = sum(d == 1), negative_only = sum(d == -1))
  }))))
interactions <- do.call(rbind, lapply(cfg$signal, function(r)
  do.call(rbind, lapply(names(weights), function(w) {
    pair <- weights[[w]]
    d <- (get_reject(pair[1L], 1, r) - get_reject(pair[1L], -1, r)) -
         (get_reject(pair[2L], 1, r) - get_reject(pair[2L], -1, r))
    cbind(data.frame(r = r, weight = w, centered = pair[1L], uncentered = pair[2L]),
          paired_stats(d))
  }))))
stopifnot(nrow(summary) == 40L, nrow(fixed) == 16L,
          nrow(signs) == 20L, nrow(interactions) == 8L,
          pool_error < 1e-12, original_score_error < 1e-10,
          new_score_error < 1e-10, record_rank_mismatch == 0L,
          independent_rank_mismatch == 0L)
rng_unchanged <- identical(rng_exists,
  exists('.Random.seed', envir = .GlobalEnv, inherits = FALSE)) &&
  (!rng_exists || identical(rng_before, get('.Random.seed', envir = .GlobalEnv)))
stopifnot(rng_unchanged)
write.csv(tests, file.path(out_dir, 'tests.csv'), row.names = FALSE)
write.csv(summary, file.path(out_dir, 'summary.csv'), row.names = FALSE)
write.csv(fixed, file.path(out_dir, 'fixed_weight_differences.csv'), row.names = FALSE)
write.csv(signs, file.path(out_dir, 'background_sign_differences.csv'), row.names = FALSE)
write.csv(interactions, file.path(out_dir, 'centering_background_interactions.csv'), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out_dir, 'session.txt'))
check <- list(completed_at = as.character(Sys.time()), replications = length(seen),
  methods = methods, settings = nrow(cells), recomputed_pvalues = nrow(tests),
  added_control_pvalues = 40000L, original_pvalues_replayed = nrow(replayed),
  original_pvalue_mismatches = sum(replayed$pvalue_calc != replayed$pvalue_saved),
  original_reject_mismatches = sum(replayed$reject_calc != replayed$reject_saved),
  record_rank_mismatches = record_rank_mismatch,
  independent_control_rank_mismatches = independent_rank_mismatch,
  maximum_pool_error = pool_error, maximum_original_score_error = original_score_error,
  maximum_independent_control_score_error = new_score_error,
  minimum_added_control_observed_copy_spacing = minimum_new_rank_spacing,
  fixed_weight_comparisons = nrow(fixed), background_sign_comparisons = nrow(signs),
  centered_background_interactions = nrow(interactions),
  new_sampling = FALSE, rng_state_unchanged = rng_unchanged,
  source_files_written = FALSE, manuscript_changed = FALSE,
  scope = 'Follow-up analyses of existing saved replications; no new data or copies.')
jsonlite::write_json(check, file.path(out_dir, 'CHECK.json'), pretty = TRUE,
                     auto_unbox = TRUE, digits = 17)
cat(jsonlite::toJSON(check, pretty = TRUE, auto_unbox = TRUE, digits = 17), '\n')
