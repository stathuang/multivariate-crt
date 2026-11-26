# Figure S3a: full paper settings.
config <- list(
  id = "main_l_gl_w",
  family = "main",
  setting = "l_gl_w",
  n = 50L,
  p = 20L,
  target = c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L),
  theta = c(0L, 1L, 2L, 3L, 4L, 6L, 8L),
  repetitions = 400L,
  M = 100L,
  L = 3L,
  p_rule = "discrete",
  response_family = "binomial",
  statistics = c("GLM-Dev", "RF"),
  baselines = c("Chi-square", "dCRT (L1)", "dCRT (RF)"),
  seed_base = 500000L,
  bandwidth = 6L,
  precision_offdiag = 0.2,
  coefficient_policy = "redraw_each_replicate",
  # Keep the verified likelihood-supremum calculation and its numerical certificates.
  binary_likelihood_policy = "binary_face_exact_dual_v3"
)
config$methods <- c(config$statistics, config$baselines)
alpha <- 0.05
plot_methods <- c("GLM-Dev", "RF", "Chi-square", "dCRT (L1)")

