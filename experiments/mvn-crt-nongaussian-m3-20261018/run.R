#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !args[1L] %in% c("prepare", "check", "pilot", "formal", "summarize", "verify")) {
  stop("Usage: Rscript run.R prepare|check|pilot|formal|summarize|verify")
}
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
task_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1L])))
setwd(task_dir)
Sys.setenv(TZ = "Asia/Singapore", OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
options(stringsAsFactors = FALSE, digits = 17)
source("config.R")
source("R/core.R")
source("R/reference.R")
phase <- args[1L]
logfile <- file.path("logs", paste0(phase, "-", format(Sys.time(), "%Y%m%d-%H%M%S"), ".log"))
sink(logfile, split = TRUE)
cat("Phase:", phase, "\nStarted:", singapore_time(), "\nDirectory:", task_dir, "\n")

if (phase == "prepare") {
  if (file.exists("seeds.rds")) stop("Seed bank exists; refusing to replace it")
  stopifnot(config$s == 12L, config$t == 8L, config$M == 199L,
            config$repetitions == 1000L, config$lambda == 1, config$alpha == 0.05)
  RNGkind(config$rng_kind, config$normal_kind, config$sample_kind)
  set.seed(config$root_seed)
  seed <- .Random.seed
  check_seed <- seed
  seed <- parallel::nextRNGStream(seed)
  cells <- cells_for(config)
  mapping <- list()
  streams <- list()
  index <- 1L
  for (stage in c("pilot", "formal")) {
    reps <- if (stage == "pilot") config$pilot_repetitions else config$repetitions
    for (c in seq_len(nrow(cells))) {
      for (b in seq_len(reps)) {
        mapping[[index]] <- data.frame(phase = stage, cell_id = c,
          distribution = cells$distribution[c], n = cells$n[c], repetition = b,
          stream_index = index, seed_key = paste(seed, collapse = ":"))
        streams[[index]] <- seed
        seed <- parallel::nextRNGStream(seed)
        index <- index + 1L
      }
    }
  }
  mapping <- do.call(rbind, mapping)
  stopifnot(!anyDuplicated(mapping$seed_key), nrow(mapping) == 12060L)
  saveRDS(list(check_seed = check_seed, streams = streams, mapping = mapping,
               config = config, prepared_at = singapore_time()), "seeds.rds", version = 3)
  write.csv(mapping, "seeds.csv", row.names = FALSE)
  write.csv(cells, "cells.csv", row.names = FALSE)
  writeLines(capture.output(sessionInfo()), "checks/sessionInfo.txt")
  dput(config, "checks/config_frozen.R")
  cat("Saved 12,060 unique streams: 60 pilot, 12,000 formal; separate check stream.\n")
}

if (phase == "check") source("R/checks.R")

if (phase %in% c("pilot", "formal")) {
  bank <- readRDS("seeds.rds")
  stopifnot(identical(bank$config, config))
  if (!file.exists("checks/PASS.rds")) stop("Run and pass checks first")
  pass <- readRDS("checks/PASS.rds")
  stopifnot(identical(pass$code_md5, tools::md5sum(c("config.R", "R/core.R", "R/reference.R", "R/checks.R", "run.R"))))
  if (phase == "formal" && !file.exists("results/pilot/timing.rds")) stop("Pilot required")
  outdir <- file.path("results", phase)
  if (dir.exists(outdir)) stop("Output exists; refusing to overwrite or silently rerun")
  dir.create(outdir)
  file.copy(c("config.R", "run.R"), outdir)
  dir.create(file.path(outdir, "R"))
  file.copy(c("R/core.R", "R/reference.R", "R/checks.R"), file.path(outdir, "R"))
  writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
  failures_file <- file.path(outdir, "failures.csv")
  warnings_file <- file.path(outdir, "warnings.csv")
  write.csv(data.frame(cell_id = integer(), repetition = integer(), stream_index = integer(),
                       message = character(), call = character()), failures_file, row.names = FALSE)
  write.csv(data.frame(cell_id = integer(), repetition = integer(), stream_index = integer(),
                       message = character()), warnings_file, row.names = FALSE)
  selected <- which(bank$mapping$phase == phase)
  start <- Sys.time()
  start_cpu <- proc.time()
  for (cid in seq_len(12L)) {
    positions <- selected[bank$mapping$cell_id[selected] == cid]
    entry <- bank$mapping[positions, ]
    cdir <- file.path(outdir, sprintf("cell-%02d-%s-n%d", cid, entry$distribution[1], entry$n[1]))
    dir.create(cdir)
    B <- nrow(entry)
    scores <- matrix(NA_real_, B, config$M + 1L,
                     dimnames = list(NULL, c("T0", paste0("T", seq_len(config$M)))))
    rows <- vector("list", B)
    cell_start <- Sys.time()
    cat("Cell", cid, entry$distribution[1], "n =", entry$n[1], "started", singapore_time(), "\n")
    for (b in seq_len(B)) {
      j <- positions[b]
      warnings_seen <- character()
      tick <- proc.time()[["elapsed"]]
      ans <- tryCatch(withCallingHandlers(
        one_replication(bank$streams[[j]], entry$n[b], entry$distribution[b], config),
        warning = function(w) {
          warnings_seen <<- c(warnings_seen, conditionMessage(w))
          invokeRestart("muffleWarning")
        }), error = function(e) e)
      base <- entry[b, c("cell_id", "distribution", "n", "repetition", "stream_index", "seed_key")]
      if (inherits(ans, "error")) {
        rows[[b]] <- cbind(base, status = "failed", T0 = NA_real_, p_value = NA_real_,
           exceedances = NA_integer_, ties = NA_integer_, reject = NA,
           d = NA_integer_, rank_D = NA_integer_, min_eigenvalue = NA_real_,
           max_eigenvalue = NA_real_, response_ss = NA_real_,
           elapsed_seconds = proc.time()[["elapsed"]] - tick,
           error = conditionMessage(ans))
        failure <- data.frame(cell_id = cid, repetition = b, stream_index = j,
          message = conditionMessage(ans), call = paste(deparse(conditionCall(ans)), collapse = " "))
        write.table(failure, failures_file, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
      } else {
        scores[b, ] <- c(ans$T0, ans$Tm)
        rows[[b]] <- cbind(base, status = "ok", T0 = ans$T0, p_value = ans$p_value,
           exceedances = ans$exceedances, ties = ans$ties, reject = ans$reject,
           d = ans$d, rank_D = ans$rank_D, min_eigenvalue = ans$min_eigenvalue,
           max_eigenvalue = ans$max_eigenvalue, response_ss = ans$response_ss,
           elapsed_seconds = proc.time()[["elapsed"]] - tick, error = "")
      }
      if (length(warnings_seen)) for (w in warnings_seen) {
        write.table(data.frame(cell_id = cid, repetition = b, stream_index = j, message = w),
          warnings_file, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
      }
      # Append each attempt immediately; a failed trial is never replaced.
      write.table(rows[[b]], file.path(cdir, "replications.csv"), sep = ",", row.names = FALSE,
                  col.names = b == 1L, append = b > 1L)
      if (b %% 100L == 0L) {
        saveRDS(list(completed = b, scores = scores[seq_len(b), , drop = FALSE]),
                file.path(cdir, "checkpoint.rds"))
        cat("  completed", b, "/", B, "\n")
      }
    }
    table <- do.call(rbind, rows)
    # RDS retains full double precision for rank verification and reproduction.
    saveRDS(list(replications = table, scores = scores, config = config,
                 elapsed_seconds = as.numeric(difftime(Sys.time(), cell_start, units = "secs"))),
            file.path(cdir, "results.rds"))
    con <- gzfile(file.path(cdir, "scores.csv.gz"), "wt")
    write.csv(data.frame(repetition = seq_len(B), scores), con, row.names = FALSE)
    close(con)
    cat("Cell", cid, "finished; success", sum(table$status == "ok"),
        "failures", sum(table$status != "ok"), "seconds",
        as.numeric(difftime(Sys.time(), cell_start, units = "secs")), "\n")
  }
  timing <- list(start = format(start, tz = "Asia/Singapore", usetz = TRUE),
    end = singapore_time(), wall_seconds = as.numeric(difftime(Sys.time(), start, units = "secs")),
    cpu = proc.time() - start_cpu, workers = 1L, threads = 1L)
  saveRDS(timing, file.path(outdir, "timing.rds"))
  dput(timing, file.path(outdir, "timing.txt"))
  if (phase == "pilot") cat("Conservative serial full-run time estimate (pilot x 200):",
    timing$wall_seconds * config$repetitions / config$pilot_repetitions, "seconds\n")
  print(timing)
}

if (phase == "summarize") source("R/summarize.R")
if (phase == "verify") source("R/verify.R")
cat("Finished:", singapore_time(), "\n")
sink()
