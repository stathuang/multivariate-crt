config <- list(id = "oracle_plugin_power_curve",
  family = "oracle_benchmark", response_family = "gaussian",
  master_seed = 20260429L, repetitions = 400L, M = 199L, lambda_ridge = 1,
  workers = 1L, theta = seq(0, 1.5, by = .25))
if (quick) config$repetitions <- 2L  # Retain every setting, signal level and copy.
alpha <- .05
settings <- data.frame(setting = c("LD-tight", "LD-moderate", "p>n target"),
  n = c(60L, 80L, 80L), s = c(45L, 50L, 60L), t = c(10L, 10L, 40L))
settings$p <- settings$s + settings$t
settings$d <- settings$n - settings$s - 1L
settings$panel <- sprintf("n=%d, s=%d, t=%d, d=%d", settings$n, settings$s, settings$t, settings$d)
config$settings <- settings
config$methods <- c("Oracle CRT", "Plug-in MLE CRT", "Plug-in CRT (d.f. adjustment)", "MVN-CRT", "LW plug-in CRT (D.10)")
# The first actual replicate seed, used by the common configuration guard.
config$seed_base <- config$master_seed + 1010000L
if (.Platform$OS.type == "windows" && config$workers > 1L) {
  warning("Parallel execution is unavailable with mclapply on Windows; using one worker.")
  config$workers <- 1L
}
if (any(!is.finite(c(config$repetitions, config$M, config$workers))) ||
    any(c(config$repetitions, config$M, config$workers) < 1) ||
    any(c(config$repetitions, config$M, config$workers) != floor(c(config$repetitions, config$M, config$workers))))
  stop("Repetitions, randomizations and workers must be positive integers")

