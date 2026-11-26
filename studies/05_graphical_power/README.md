# Graphical power examples

Supplement F illustrates the direction effect and compares the original and centred statistics.

| `--experiment` | Paper results | Formal calculation |
|---|---|---|
| `direction` | Table S4 | 500 datasets per sample size; 999 copies |
| `independent` | Figure S9; Tables S5–S6 | 100 datasets per sample size; 19 copies |
| `connected` | Figure S10; Tables S7–S8 | 2,500 datasets; 99 copies |

The first two experiments use sample sizes 40, 120, and 360. The connected experiment uses 80 observations, eight targets, and 112 background covariates. All use one sweep per direction and `p <= 0.05`.

The connected calculation evaluates all five statistics on the same copies, including the two uncentred controls. Figure S10 shows the three statistics selected in the paper; the tables retain all five. The direction and independent-target examples keep their distinct copy counts.

## Run the full calculations

From the repository root:

```sh
Rscript --vanilla studies/05_graphical_power/run.R --experiment direction
Rscript --vanilla studies/05_graphical_power/run.R --experiment independent
Rscript --vanilla studies/05_graphical_power/run.R --experiment connected
```

Add `--quick` to a calculation command for fewer replications. The output is `results/generated/05_graphical_power/` for full runs, or `results/quick/05_graphical_power/` for quick runs. No saved paper result is overwritten.

## Summarize and plot

After all full calculations in this study:

```sh
.venv/bin/python scripts/summarize.py 05_graphical_power --source generated
Rscript --vanilla scripts/figures.R 05_graphical_power generated
```

To use the saved paper records instead, omit `--source generated` and the final `generated` argument. The source is always recorded in the output directory. Study 03 has a numerical table and no figure.
