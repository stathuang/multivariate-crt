# Fixed before checks, pilot, and formal simulation. No external packages.
config <- list(
  experiment = "mvn-crt-nongaussian-m3-20261018",
  distributions = c("Gaussian", "t5", "t3"),
  sample_sizes = c(50L, 100L, 200L, 400L),
  s = 12L, t = 8L, rho = 0.5,
  lambda = 1, M = 199L, alpha = 0.05, repetitions = 1000L,
  pilot_repetitions = 5L,
  rng_kind = "L'Ecuyer-CMRG", normal_kind = "Inversion",
  sample_kind = "Rejection", root_seed = 20261018L,
  workers = 1L,
  method = "Gaussian MVN-CRT; exact conditional distribution of ridge gain",
  response = "rowSums(Z)/sqrt(12) + (Z[,1]^2 - 1)/sqrt(2) + epsilon",
  target_mean = "(Z[,1:8] + Z[,2:9])/sqrt(2)",
  error = "G for Gaussian; sqrt((nu-2)/V_i)*G_i for t_nu; one V_i per row",
  background_fit = "intercept and 12 original Z columns; no quadratic terms"
)
