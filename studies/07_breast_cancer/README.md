# Breast cancer

Section 8.2 reports the four observed p-values. Table S12 and Figure S12 report fitted-null rejection rates and the LRT p-value histogram.

Prepare `data/processed/bcr_erplus60.rds` first; see [data instructions](../../data/README.md). The analysis uses 209 ER-positive patients, three background probes, and 18 target probes.

To fetch and prepare only the breast-cancer inputs, run from the repository root:

```sh
python3 scripts/fetch_data.py data/raw --dataset bcr
Rscript --vanilla scripts/prepare_data.R raw data/raw bcr
```

`observed_settings.json` and `null_settings.json` specify the two calculations. The full null experiment has 400 replications. The quick null run uses replications 1 and 74; the second exercises the random-forest empty-class repair.

The runner preserves the original sequence of random draws for the four methods. It then recomputes the multinomial likelihood-ratio test with the certified numerical implementation and repairs an RF failure using its saved RNG state. These steps are needed to reproduce the adopted numerical results. Each replication saves the likelihood diagnostics. The dCRT comparator uses the original numeric response codes and Gaussian nuisance fits.

The target/background selection and transformations precede these tests; the null study holds the chosen sets fixed.

## Run the full calculations

From the repository root:

```sh
Rscript --vanilla studies/07_breast_cancer/run.R --experiment observed
Rscript --vanilla studies/07_breast_cancer/run.R --experiment null
```

Add `--quick` to a calculation command for fewer replications. The output is `results/generated/07_breast_cancer/` for full runs, or `results/quick/07_breast_cancer/` for quick runs. No saved paper result is overwritten.

## Summarize and plot

After all full calculations in this study:

```sh
.venv/bin/python scripts/summarize.py 07_breast_cancer --source generated
Rscript --vanilla scripts/figures.R 07_breast_cancer generated
```

To use the saved paper records instead, omit `--source generated` and the final `generated` argument. The source is always recorded in the output directory. Study 03 has a numerical table and no figure.
