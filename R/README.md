# Shared R functions

- `sampler.R`: Gaussian graphical hub-and-spoke sampling and rank p-values.
- `data.R`, `binary_data.R`, and `extended_designs.R`: simulation data generators.
- `linear_tests.R`, `gaussian_tests.R`, and `binary_tests.R`: the G-CRT statistics in Table S1.
- `dcrt.R`, `comparators.R`, and `hdi_cache.R`: comparison methods and their fitted components.
- `experiment.R`, `binary_experiment.R`, and `supergraph_experiment.R`: one-replication evaluation in the specified method order.
- `binary_likelihood.R`, `multinomial_likelihood.R`, and `bcr_tests.R`: numerical likelihood calculations and diagnostics.
- `stock_data.R` and `bcr_data.R`: application input preparation.
- `application_diagnostics.R`: stock graph and goodness-of-fit calculations.
- `cgm_diagnostics.R` and `result_diagnostics.R`: method diagnostics and failure records.
- `plot_helpers.R`: shared plot styles.
- `setup.R` and `run_helpers.R`: entry-point options, output paths, warning/error handling, and CSV export.

The Gaussian comparison and graphical-power studies keep their specialised mathematical functions beside their settings. The stock and breast-cancer studies retain their application-specific statistic and baseline implementations because their forest sizes, response conventions, and recorded random-number streams differ.
