# E19: display the fitted stock graph; no simulation is required.
# 1. Parameters ---------------------------------------------------------------
experiment_id <- "E19"
config <- list(id = "stock_graph", family = "application", response_family = "gaussian",
               repetitions = 1L, seed_base = 0L, observed_seed = 1L)
# 2. Load inputs and record this rendering ------------------------------------
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

# 3. Draw the accepted graph and sector legend --------------------------------
plot_stock_graph(stock, file.path(output_dir, "stocks_G.pdf"))
# 4. Save the graph data and rendering identity --------------------------------
saveRDS(list(graph = stock$graph, stocks = stock$stocks, sectors = stock$sectors),
        file.path(output_dir, "graph.rds"))
run$rng_after <- .Random.seed
run$state <- "completed"
saveRDS(run, file.path(output_dir, "run.rds"))
message("Saved E19 graph: ", file.path(output_dir, "stocks_G.pdf"))
