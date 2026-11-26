# Figure S4b: full paper settings.
config <- list(
  id = "main_h_gl_m",
  family = "main",
  setting = "h_gl_m",
  n = 80L,
  p = 120L,
  # Test columns 1:8 conditional on the rest; theta scales transformed columns 1:10.
  # Columns 9:10 are conditioning variables; changing target does not move the signal block.
  target = c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L),
  theta = c(0L, 1L, 2L, 3L, 5L, 7L),
  repetitions = 800L,
  M = 100L,
  L = 3L,
  p_rule = "discrete",
  response_family = "binomial",
  statistics = c("RF", "RF-D", "RF-RR"),
  baselines = c("CGM", "dCRT (L1)", "dCRT (RF)"),
  seed_base = 800000L,
  bandwidth = 6L,
  precision_offdiag = 0.2,
  coefficient_policy = "redraw_each_replicate",
  # Preserve the default CGM filter and record its typed failures.
  cgm_diagnostic_policy = "cgm_retained_count_guard_v1"
)
config$methods <- c(config$statistics, config$baselines)
alpha <- 0.05
plot_methods <- c("RF", "RF-D", "RF-RR", "CGM", "dCRT (RF)")

