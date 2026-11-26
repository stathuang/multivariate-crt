# Covariate model misspecification

Figure S8 and Table S3 assess the Gaussian MVN-CRT under Gaussian, t5, and t3 covariate errors.

`settings.R` contains all parameters: sample sizes 50, 100, 200, and 400; 12 background covariates; eight targets; 1,000 replications per cell; and 199 randomised copies.

Each t observation has one shared scale across its target coordinates. The ridge penalty is one. P-values are conservative rank p-values. The summary uses Wilson intervals at nominal levels 0.005 through 0.100.

`run.R` directly reproduces the formal experiment. Its random streams advance past the one check stream and 60 pilot streams used before the original formal calculation. The pilots do not have to be rerun. A quick run retains the first two formal replications in each of the 12 cells.

## Run the full calculations

From the repository root:

```sh
Rscript --vanilla studies/04_misspecification/run.R
```

Add `--quick` to a calculation command for fewer replications. The output is `results/generated/04_misspecification/` for full runs, or `results/quick/04_misspecification/` for quick runs. No saved paper result is overwritten.

## Summarize and plot

After all full calculations in this study:

```sh
.venv/bin/python scripts/summarize.py 04_misspecification --source generated
Rscript --vanilla scripts/figures.R 04_misspecification generated
```

To use the saved paper records instead, omit `--source generated` and the final `generated` argument. The source is always recorded in the output directory. Study 03 has a numerical table and no figure.
