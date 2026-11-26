# Figure S7b: full paper settings.
config <- list(
  id = "t3_h_l_w",
  family = "t3",
  setting = "h_l_w",
  n = 80L,
  p = 120L,
  target = c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L),
  theta = c(0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5),
  repetitions = 400L,
  M = 100L,
  L = 3L,
  p_rule = "discrete",
  response_family = "gaussian",
  statistics = c("LM-L1-R-SST", "LM-L1-R-SSR"),
  baselines = c("De-Lasso", "CGM", "dCRT (L1)", "dCRT (RF)"),
  seed_base = 1400000L,
  bandwidth = 6L,
  precision_offdiag = 0.2,  # Retained config metadata; X is generated from t innovations.
  t3_lag_coefficients = c(0.2, 0.2, -0.2, -0.2, 0.2, 0.2),
  t3_df = 3,
  t3_noise_scale = 0.7,
  t3_initialization = "manuscript_lag_order",
  coefficient_policy = "redraw_each_replicate",
  # Comparator execution and diagnostic settings; see docs/CACHE_CANDIDATE.md for De-Lasso.
  de_lasso_policy = "one-nodewise-fit-per-design-v1",
  de_lasso_seed_policy = "separate-nodewise-response-seeds-v1",
  de_lasso_suppress_unused_group_testing = TRUE,
  cgm_diagnostic_policy = "cgm_retained_count_guard_v1"
)
config$methods <- c(config$statistics, config$baselines)
alpha <- 0.05
plot_methods <- c("LM-L1-R-SST", "LM-L1-R-SSR", "De-Lasso", "dCRT (L1)")

