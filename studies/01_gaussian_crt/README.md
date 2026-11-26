# Gaussian CRT comparisons

Table 1 and Figure 1 compare the oracle, two regression plug-in CRTs, MVN-CRT, and the Ledoit–Wolf plug-in CRT.

| Calculation | Settings | Output |
|---|---|---|
| `size` | 3 dimension settings; 5,000 replications; 999 copies | Table 1 |
| `power` | Same dimensions; 400 replications per signal; 199 copies | Figure 1 |

The dimension triples `(n, s, t)` are `(60,45,10)`, `(80,50,10)`, and `(80,60,40)`. Settings are in `size_settings.R` and `power_settings.R`. Rejections use `p <= 0.05`; reported Monte Carlo standard errors use the binomial convention. Raw size results also retain the two other score statistics evaluated in the original experiment.

## Run the full calculations

From the repository root:

```sh
Rscript --vanilla studies/01_gaussian_crt/run.R --experiment size
Rscript --vanilla studies/01_gaussian_crt/run.R --experiment power
```

Add `--quick` to a calculation command for fewer replications. The output is `results/generated/01_gaussian_crt/` for full runs, or `results/quick/01_gaussian_crt/` for quick runs. No saved paper result is overwritten.

## Summarize and plot

After all full calculations in this study:

```sh
.venv/bin/python scripts/summarize.py 01_gaussian_crt --source generated
Rscript --vanilla scripts/figures.R 01_gaussian_crt generated
```

To use the saved paper records instead, omit `--source generated` and the final `generated` argument. The source is always recorded in the output directory. Study 03 has a numerical table and no figure.
