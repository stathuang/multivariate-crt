# Stock returns

Section 8.1 uses two stock groups and a five-trading-day DIA response.

| `--experiment` | Paper results |
|---|---|
| `observed-T1`, `observed-T2` | Table S11 |
| `null-T1`, `null-T2` | Table S10 |
| `noise-05`, `noise-10`, `noise-15`, `noise-20`, `noise-25` | Figure 2 |
| `graph` | Figure S11 |
| `gof` | Table S9 |

Prepare `data/processed/stock.rds` first; see [data instructions](../../data/README.md). The input has 92 windows and 103 stocks. Observed tests use 1,000 copies; null and noise experiments use 400 copies and 400 replications. Random forests use 1,000 trees. The graph goodness-of-fit calculation uses 1,000 copies, including in quick mode.

To fetch and prepare only the stock inputs, run from the repository root:

```sh
python3 scripts/fetch_data.py data/raw --dataset stock
Rscript --vanilla scripts/prepare_data.R raw data/raw stock
```

`settings.json` fixes method order and seeds. The fitted-null experiments condition on the prepared data and fitted null models. `R/stock_transform.R` supports the independent reconstruction of the transformed stock matrix from saved prices.

## Run the full calculations

From the repository root:

```sh
Rscript --vanilla studies/06_stock_returns/run.R --experiment observed-T1
Rscript --vanilla studies/06_stock_returns/run.R --experiment observed-T2
Rscript --vanilla studies/06_stock_returns/run.R --experiment null-T1
Rscript --vanilla studies/06_stock_returns/run.R --experiment null-T2
Rscript --vanilla studies/06_stock_returns/run.R --experiment noise-05
Rscript --vanilla studies/06_stock_returns/run.R --experiment noise-10
Rscript --vanilla studies/06_stock_returns/run.R --experiment noise-15
Rscript --vanilla studies/06_stock_returns/run.R --experiment noise-20
Rscript --vanilla studies/06_stock_returns/run.R --experiment noise-25
Rscript --vanilla studies/06_stock_returns/run.R --experiment graph
Rscript --vanilla studies/06_stock_returns/run.R --experiment gof
```

Add `--quick` to a calculation command for fewer replications. The output is `results/generated/06_stock_returns/` for full runs, or `results/quick/06_stock_returns/` for quick runs. No saved paper result is overwritten.

## Summarize and plot

After all full calculations in this study:

```sh
.venv/bin/python scripts/summarize.py 06_stock_returns --source generated
Rscript --vanilla scripts/figures.R 06_stock_returns generated
```

To use the saved paper records instead, omit `--source generated` and the final `generated` argument. The source is always recorded in the output directory. Study 03 has a numerical table and no figure.
