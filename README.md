# Multivariate conditional randomisation tests

Replication code for **Multivariate conditional randomisation tests with unknown covariate parameters**, by Xiaotong Lin, Jie Xie, Fangqiao Tian, and Dongming Huang.

Start with the paper result you want to reproduce. Each study contains its settings and a calculation entry point. Shared R functions are in `R/`. Saved paper results are in `results/reference/`; new calculations go to `results/generated/`.

| Study | Paper results | Instructions |
|---|---|---|
| Gaussian CRT comparisons | Table 1; Figure 1; Section 7.1 | [01_gaussian_crt](studies/01_gaussian_crt/README.md) |
| Graphical CRT simulations | Figures S1–S7; Section 7.2 and Supplement E | [02_graphical_crt](studies/02_graphical_crt/README.md) |
| Joint and coordinate tests | Table S2; Supplement E.5 | [03_joint_tests](studies/03_joint_tests/README.md) |
| Covariate model misspecification | Figure S8; Table S3 | [04_misspecification](studies/04_misspecification/README.md) |
| Graphical power examples | Figures S9–S10; Tables S4–S8; Supplement F | [05_graphical_power](studies/05_graphical_power/README.md) |
| Stock returns | Figure 2; Figure S11; Tables S9–S11; Section 8.1 | [06_stock_returns](studies/06_stock_returns/README.md) |
| Breast cancer | Figure S12; Table S12; Section 8.2 p-values | [07_breast_cancer](studies/07_breast_cancer/README.md) |

Table S1 defines test statistics; its implementations are identified in the graphical study README.

## Set up

Run commands from the repository root. The recorded environment uses Python 3.12 and R 4.3.1.

```sh
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
Rscript --vanilla scripts/check_R_environment.R
```

See [environment/README.md](environment/README.md) for R installation and recorded package versions. The scripts do not install packages during an analysis.

## Recompute an experiment

For example, the following commands regenerate Table 1 and Figure 1. The power experiment uses 400 replications per setting and signal level:

```sh
Rscript --vanilla studies/01_gaussian_crt/run.R --experiment size
Rscript --vanilla studies/01_gaussian_crt/run.R --experiment power
.venv/bin/python scripts/summarize.py 01_gaussian_crt --source generated
Rscript --vanilla scripts/figures.R 01_gaussian_crt generated
```

Each study README lists its complete commands. [reproduce.py](reproduce.py) provides an optional entry point for running every calculation in a study:

```sh
.venv/bin/python reproduce.py 01_gaussian_crt
# Run all seven studies, after preparing application inputs:
.venv/bin/python reproduce.py all
```

These are full simulations and can take substantial time. Add `--quick` to check an entry point with fewer datasets while retaining the formal copy count and method settings. Quick outputs go to `results/quick/` and are not used as paper results. A calculation refuses to overwrite an existing nonempty output directory.

## Rebuild figures and tables from saved paper results

This faster route recalculates summaries from saved replication records; it does not rerun the simulations or model fits.

```sh
.venv/bin/python scripts/summarize.py all
Rscript --vanilla scripts/figures.R all
.venv/bin/python scripts/population.py
```

Outputs appear under `results/tables/reference/` and `results/figures/reference/`, with paper labels such as `Table1.tex` and `FigureS1a.pdf`. Each table also has a CSV with labelled columns and row groups; Table S8 has separate files for its two panels. Summary CSV files retain method names, counts, and uncertainty estimates. The TeX files contain numerical rows for the corresponding paper table. The summary command checks saved results against the paper's displayed numerical precision. Plot styling can differ from the paper, and the stock graph layout can vary with the igraph version.

## Application inputs and verification

See [data/README.md](data/README.md) for inputs and preparation. The fixed DIA price snapshot is included in `data/snapshots/`; the other six inputs are downloaded from pinned sources. Downloaded inputs under `data/raw/`, processed application objects, software environments, and newly generated results are excluded from Git. Saved replication results contain p-values and diagnostics, not individual patient or price matrices.

[tests/README.md](tests/README.md) describes the checks. [CREDITS.md](CREDITS.md) identifies external code and data sources.
