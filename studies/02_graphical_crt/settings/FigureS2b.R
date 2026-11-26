# Figure S2b: full paper settings.
config <- list(
  id = "main_h_l_m",
  family = "main",
  setting = "h_l_m",
  n = 80L,
  p = 120L,
  # Test columns 1:8 conditional on the rest. theta scales the nonlinear block in columns 1:10.
  # Columns 9:10 are conditioning variables; changing target does not move the signal block.
  target = c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L),
  theta = c(0L, 1L, 2L, 3L, 5L, 7L),
  repetitions = 400L,
  M = 100L,
  L = 3L,
  p_rule = "discrete",
  response_family = "gaussian",
  statistics = c("RF", "RF-D", "RF-RR"),
  baselines = c("De-Lasso", "CGM", "dCRT (L1)", "dCRT (RF)"),
  seed_base = 400000L,
  bandwidth = 6L,
  precision_offdiag = 0.2,
  coefficient_policy = "redraw_each_replicate",
  # Comparator execution and diagnostic settings; see docs/CACHE_CANDIDATE.md for De-Lasso.
  de_lasso_policy = "one-nodewise-fit-per-design-v1",
  de_lasso_seed_policy = "separate-nodewise-response-seeds-v1",
  de_lasso_suppress_unused_group_testing = TRUE,
  cgm_diagnostic_policy = "cgm_retained_count_guard_v1"
)
config$methods <- c(config$statistics, config$baselines)
alpha <- 0.05
plot_methods <- c("RF", "RF-D", "RF-RR", "De-Lasso", "dCRT (RF)")

