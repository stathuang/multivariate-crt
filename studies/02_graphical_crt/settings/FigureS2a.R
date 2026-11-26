# Figure S2a: full paper settings.
config <- list(
  id = "main_l_l_m",
  family = "main",
  setting = "l_l_m",
  n = 50L,
  p = 20L,
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
  baselines = c("F-test", "dCRT (L1)", "dCRT (RF)"),
  seed_base = 300000L,
  bandwidth = 6L,
  precision_offdiag = 0.2,
  coefficient_policy = "redraw_each_replicate"
)
config$methods <- c(config$statistics, config$baselines)
alpha <- 0.05
plot_methods <- c("RF", "RF-D", "RF-RR", "F-test", "dCRT (RF)")

