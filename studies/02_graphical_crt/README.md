# Graphical CRT simulations

These experiments compare G-CRT statistics and the available competing methods under the Gaussian graphical designs in Section 7.2 and Supplement E.

| `--figure` | Design |
|---|---|
| `S1a`, `S1b` | Gaussian linear, low/high dimensions |
| `S2a`, `S2b` | Gaussian nonlinear, low/high dimensions |
| `S3a`, `S3b` | Binary response, first design pair |
| `S4a`, `S4b` | Binary response, second design pair |
| `S5a`, `S5b` | One target variable |
| `S6` | True graph and graphs with additional edges |
| `S7a`, `S7b` | Multivariate t3 covariates |

Each panel has an explicit file in `settings/`. The experiments use 100 copies and three sweeps. S4 uses 800 replications; the other panels use 400. Rejections use `p < 0.05`. Failures remain in the raw CSV. The saved high-dimensional logistic CGM cell at signal 3.5 has 398 valid results from 400 attempts.

Table S1 defines the statistics implemented in `R/linear_tests.R`, `R/gaussian_tests.R`, and `R/binary_tests.R`. The sampler is in `R/sampler.R`. The runner preserves the per-replication method order; for S6 it finishes all tests for one graph before sampling from the next graph.

## Run the full calculations

From the repository root:

```sh
Rscript --vanilla studies/02_graphical_crt/run.R --figure S1a
Rscript --vanilla studies/02_graphical_crt/run.R --figure S1b
Rscript --vanilla studies/02_graphical_crt/run.R --figure S2a
Rscript --vanilla studies/02_graphical_crt/run.R --figure S2b
Rscript --vanilla studies/02_graphical_crt/run.R --figure S3a
Rscript --vanilla studies/02_graphical_crt/run.R --figure S3b
Rscript --vanilla studies/02_graphical_crt/run.R --figure S4a
Rscript --vanilla studies/02_graphical_crt/run.R --figure S4b
Rscript --vanilla studies/02_graphical_crt/run.R --figure S5a
Rscript --vanilla studies/02_graphical_crt/run.R --figure S5b
Rscript --vanilla studies/02_graphical_crt/run.R --figure S6
Rscript --vanilla studies/02_graphical_crt/run.R --figure S7a
Rscript --vanilla studies/02_graphical_crt/run.R --figure S7b
```

Add `--quick` to a calculation command for fewer replications. The output is `results/generated/02_graphical_crt/` for full runs, or `results/quick/02_graphical_crt/` for quick runs. No saved paper result is overwritten.

## Summarize and plot

After all full calculations in this study:

```sh
.venv/bin/python scripts/summarize.py 02_graphical_crt --source generated
Rscript --vanilla scripts/figures.R 02_graphical_crt generated
```

To use the saved paper records instead, omit `--source generated` and the final `generated` argument. The source is always recorded in the output directory. Study 03 has a numerical table and no figure.
