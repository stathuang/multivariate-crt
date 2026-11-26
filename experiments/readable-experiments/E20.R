# E20: stock graph goodness of fit; retain the complete 1,000-copy calculation.
# 1. Parameters ---------------------------------------------------------------
experiment_id <- "E20"
config <- list(id = "stock_graph_gof", family = "application", response_family = "gaussian",
               repetitions = 1L, seed_base = 0L, observed_seed = 2023L, M = 1000L, L = 3L)
# 2. Load inputs and save sources before computing ----------------------------
experiment_dir <- if (file.exists("Experiments.Rproj")) getwd() else
  file.path(getwd(), "experiments", "readable-experiments")
experiment_dir <- normalizePath(experiment_dir, mustWork = TRUE)
source(file.path(experiment_dir, "R/run_tools.R"), local = TRUE)
loaded_sources <- load_experiment_sources(experiment_dir, environment(), stock = TRUE)
input <- read_prepared_application(file.path(experiment_dir, "data/prepared/stock.rds"))
stock <- input$data
config$data_sha256 <- input$sha256
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
context <- capture_run_context(config, loaded_sources)
output_dir <- new_output_directory(experiment_dir, experiment_id, FALSE)
run <- start_run_record(config, list(context = context), loaded_sources, output_dir, FALSE, .05)
set.seed(config$observed_seed)

# 3. Compute the two analytical diagnostics -----------------------------------
X <- stock$X
graph <- stock$graph
vv <- stock_VV(X, graph)
dp <- stock_DP(X, graph)

# 4. Generate graph-model copies and evaluate the sum of local F statistics ----
copies <- mscrt_sample(X, graph, target = seq_len(ncol(X)), M = config$M,
                       L = config$L, reduced = FALSE)
observed <- stock_F_sum(X, graph)
randomized <- numeric(config$M)
for (m in seq_len(config$M)) {
  randomized[m] <- stock_F_sum(copies$X_copies[[m]], graph)
  if (m %% 100L == 0L) message("Evaluated ", m, " / ", config$M, " copies")
}
f_sum_pvalue <- rank_pvalue(observed, randomized, rule = "discrete")
result <- list(pvalues = c(VV = vv, DP = dp, F_sum = f_sum_pvalue),
               observed = observed, randomized = randomized, M = config$M, L = config$L)
.stock_diagnostic_pvalues(result$pvalues)

# 5. Save the original table style and machine-readable diagnostics -----------
write_stock_gof_table(result, file.path(output_dir, "Stock_GoF.tex"))
saveRDS(result, file.path(output_dir, "result.rds"))
write_report_csv(data.frame(method = names(result$pvalues), pvalue = unname(result$pvalues)),
                 file.path(output_dir, "results.csv"))
write_report_csv(data.frame(copy = 0:config$M, statistic = c(observed, randomized)),
                 file.path(output_dir, "statistics.csv"))
# Retain the first few copies for direct formula inspection without duplicating the whole array.
saveRDS(list(X = X, graph = graph, copies = copies$X_copies[1:3]), file.path(output_dir, "example.rds"))
run$rng_after <- .Random.seed
run$state <- "completed"
saveRDS(run, file.path(output_dir, "run.rds"))
print(result$pvalues)
message("Saved E20 outputs in: ", output_dir)
