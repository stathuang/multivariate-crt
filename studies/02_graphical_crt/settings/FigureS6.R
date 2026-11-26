# Figure S6: full paper settings.
config <- list(
  id = "supergraph_h_l_w",
  family = "supergraph",
  setting = "h_l_w",
  n = 80L,
  p = 120L,
  target = 1:8,
  theta = c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5),
  repetitions = 400L,
  M = 100L,
  L = 3L,
  p_rule = "discrete",
  response_family = "gaussian",
  statistics = c("LM-L1-R-SST", "LM-L1-R-SSR"),
  baselines = character(),
  seed_base = 1000000L,
  bandwidth = 6L,                    # True Gaussian design graph.
  precision_offdiag = 0.2,
  graph_bandwidths = c(6L, 12L, 18L), # Sampling graphs, in execution order.
  coefficient_policy = "redraw_each_replicate"
)
alpha <- 0.05

