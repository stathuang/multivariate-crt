# Figure S1a: full paper settings.
config <- list(
  id = "main_l_l_w",
  family = "main",
  setting = "l_l_w",
  n = 50L,
  p = 20L,
  target = c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L),
  theta = c(0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5),
  repetitions = 400L,
  M = 100L,
  L = 3L,
  p_rule = "discrete",
  response_family = "gaussian",
  statistics = c("LM-SST", "LM-SSR"),
  baselines = c("F-test", "dCRT (L1)", "dCRT (RF)"),
  seed_base = 100000L,
  bandwidth = 6L,
  precision_offdiag = 0.2,
  coefficient_policy = "redraw_each_replicate"
)
config$methods <- c(config$statistics, config$baselines)
alpha <- 0.05
plot_methods <- c("LM-SST", "LM-SSR", "F-test", "dCRT (L1)")

