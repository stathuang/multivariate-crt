bank <- readRDS("seeds.rds")
paths <- list.files("results/formal", pattern = "^cell-", full.names = TRUE)
stopifnot(length(paths) == 12L)
verified <- list()
all_ids <- character()
all_streams <- integer()
summary_saved <- read.csv("results/summary.csv")
for (i in seq_along(paths)) {
  result <- readRDS(file.path(paths[i], "results.rds"))
  rows <- result$replications
  T <- result$scores
  stopifnot(identical(result$config, config), nrow(rows) == 1000L,
    all(rows$repetition == 1:1000), ncol(T) == 200L, nrow(T) == 1000L)
  all_ids <- c(all_ids, paste(rows$cell_id, rows$repetition, sep = ":"))
  all_streams <- c(all_streams, rows$stream_index)
  ok <- rows$status == "ok"
  # Independently reconstruct every p-value from full-precision saved scores.
  exceedances <- rowSums(T[ok, -1, drop = FALSE] >= T[ok, 1])
  p <- (1 + exceedances) / 200
  stopifnot(identical(as.numeric(p), as.numeric(rows$p_value[ok])),
    identical(as.integer(exceedances), as.integer(rows$exceedances[ok])),
    all(rows$reject[ok] == (p <= 0.05)), all(p >= 0.005 & p <= 1),
    max(abs(200 * p - round(200 * p))) < 1e-12,
    all(rows$rank_D[ok] == 13L), all(rows$d[ok] == rows$n[ok] - 13L))
  expected_entries <- bank$mapping[rows$stream_index, ]
  stopifnot(all(expected_entries$phase == "formal"),
            all(expected_entries$seed_key == rows$seed_key),
            all(expected_entries$n == rows$n),
            all(expected_entries$distribution == rows$distribution))
  # Reproduce three predetermined trials per cell; no new formal data added.
  replay_error <- 0
  for (b in c(1L, 500L, 1000L)) {
    if (!ok[b]) next
    j <- rows$stream_index[b]
    replay <- one_replication(bank$streams[[j]], rows$n[b], rows$distribution[b], config)
    replay_error <- max(replay_error, abs(c(replay$T0, replay$Tm) - T[b, ]))
    stopifnot(identical(c(replay$T0, replay$Tm), as.numeric(T[b, ])))
  }
  ss <- summary_saved[summary_saved$cell_id == rows$cell_id[1], ]
  N <- sum(ok)
  k <- sum(p <= 0.05)
  z <- qnorm(0.975)
  estimate <- k / N
  denom <- 1 + z^2 / N
  centre <- (estimate + z^2 / (2 * N)) / denom
  half <- z / denom * sqrt(estimate * (1 - estimate) / N + z^2 / (4 * N^2))
  stopifnot(ss$completed == N, ss$rejections == k, ss$failed == sum(!ok),
    abs(ss$rejection_rate - estimate) < 1e-14,
    abs(ss$wilson_lower - (centre - half)) < 1e-14,
    abs(ss$wilson_upper - (centre + half)) < 1e-14)
  verified[[i]] <- data.frame(cell_id = rows$cell_id[1], attempted = nrow(rows),
     verified_p_values = N, failures = sum(!ok), replayed = sum(ok[c(1, 500, 1000)]),
     maximum_replay_error = replay_error, summary_verified = TRUE)
}
stopifnot(length(all_ids) == 12000L, !anyDuplicated(all_ids), !anyDuplicated(all_streams),
          setequal(all_streams, which(bank$mapping$phase == "formal")))
combined <- read.csv("results/replications.csv")
stopifnot(nrow(combined) == 12000L)
verification <- do.call(rbind, verified)
write.csv(verification, "checks/formal_verification.csv", row.names = FALSE)
print(verification, row.names = FALSE)
cat("All 12,000 attempts accounted for; p-values and summaries reconstructed; 36 exact replays.\n")
