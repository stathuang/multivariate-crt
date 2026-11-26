# Joint and coordinate tests

Table S2 compares the coordinate test, ridge, MaxCoefT, and their Bonferroni combination under G-CRT and oracle sampling.

`settings.json` defines 14 cells: two designs and seven signal settings, each with 500 replications, 999 copies, and three sweeps. The sample size is 80, with 120 covariates and eight targets. Replication identifiers run from 10000 to 10499.

`run.py` first writes `ridge.jsonl`, then computes MaxCoefT on the same datasets and writes `maxcoeft.jsonl`. The two stages can also be run separately with `--step ridge` and `--step max`. Data hashes ensure that the stages use identical datasets. Each method has its fixed random-number stream; the response Lasso cross-validation has a separate seed.

The seven displayed p-values use `p <= 0.05`. `scripts/population.py` also reconstructs the 500 background datasets used for the numerical population approximations discussed with this comparison.

## Run the full calculations

From the repository root:

```sh
.venv/bin/python studies/03_joint_tests/run.py
```

Add `--quick` to a calculation command for fewer replications. The output is `results/generated/03_joint_tests/` for full runs, or `results/quick/03_joint_tests/` for quick runs. No saved paper result is overwritten.

## Summarize and plot

After all full calculations in this study:

```sh
.venv/bin/python scripts/summarize.py 03_joint_tests --source generated
Rscript --vanilla scripts/figures.R 03_joint_tests generated
```

To use the saved paper records instead, omit `--source generated` and the final `generated` argument. The source is always recorded in the output directory. Study 03 has a numerical table and no figure.
