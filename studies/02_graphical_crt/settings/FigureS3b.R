# Figure S3b: full paper settings.
config <- list(
  id = "main_h_gl_w",
  family = "main",
  setting = "h_gl_w",
  n = 80L,
  p = 120L,
  target = c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L),
  theta = c(0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5),
  repetitions = 400L,
  M = 100L,
  L = 3L,
  p_rule = "discrete",
  response_family = "binomial",
  statistics = c("GLM-L1-D", "GLM-L1-R-SST", "RF"),
  baselines = c("CGM", "dCRT (L1)", "dCRT (RF)"),
  seed_base = 600000L,
  bandwidth = 6L,
  precision_offdiag = 0.2,
  coefficient_policy = "redraw_each_replicate",
  # Preserve the default CGM filter and record its typed failures.
  cgm_diagnostic_policy = "cgm_retained_count_guard_v1"
)
config$methods <- c(config$statistics, config$baselines)
alpha <- 0.05
plot_methods <- c("GLM-L1-D", "GLM-L1-R-SST", "RF", "CGM", "dCRT (L1)")

