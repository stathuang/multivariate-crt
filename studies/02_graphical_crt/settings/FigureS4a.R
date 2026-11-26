# Figure S4a: full paper settings.
config <- list(
  id = "main_l_gl_m",
  family = "main",
  setting = "l_gl_m",
  n = 50L,
  p = 20L,
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
  baselines = c("Chi-square", "dCRT (L1)", "dCRT (RF)"),
  seed_base = 700000L,
  bandwidth = 6L,
  precision_offdiag = 0.2,
  coefficient_policy = "redraw_each_replicate",
  # Keep the verified likelihood-supremum calculation and its numerical certificates.
  binary_likelihood_policy = "binary_face_exact_dual_v3"
)
config$methods <- c(config$statistics, config$baselines)
alpha <- 0.05
plot_methods <- c("RF", "RF-D", "RF-RR", "Chi-square", "dCRT (RF)")

