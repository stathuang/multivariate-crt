# Optional entry: Rscript experiments/readable-experiments/run.R E04 --trial
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || !args[1L] %in% c(paste0("E", sprintf("%02d", c(3:15, 19:20, 25))), "E23-E24") ||
    length(setdiff(args[-1L], "--trial")))
  stop("Usage: Rscript run.R E03|...|E15|E19|E20|E23-E24|E25 [--trial]")
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
setwd(dirname(normalizePath(script, mustWork = TRUE)))
options(readable.trial = "--trial" %in% args)
source(paste0(args[1L], ".R"), local = TRUE)
