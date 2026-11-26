# Replication code for multivariate conditional randomisation tests

This repository contains the replication code for **Multivariate conditional randomisation tests with unknown covariate parameters**. Its scope is the main text and supplementary material actually loaded by the author's v25 `Main.tex`, checked on 20 November 2026. It covers 14 figures (20 PDF assets), 12 numerical tables, one table defining statistics, and the associated numerical statements. The unloaded group-selection extension is excluded.

Figure and table references in this README and RESULT_MAP.md were updated on 26 November 2026. The computational validation and data provenance remain unchanged.

Start with [RESULT_MAP.md](RESULT_MAP.md) for a result-to-command map and [VALIDATION.md](VALIDATION.md) for the checks actually performed. The application instructions below give the exact data sources and preparation commands. [DCRT_ADAPTATIONS.md](DCRT_ADAPTATIONS.md) records the comparator source and executable reconstruction steps. [THIRD_PARTY.md](THIRD_PARTY.md) separates review use from public redistribution. This repository contains no patient-level data, market-price matrices, credentials, or installed software libraries.

## Working with the repository

The source files are tracked individually in Git. The repository is private. Figshare sharing is deferred; the separate review archive contains all exact application inputs.

After cloning, edit the files directly, stage the intended changes, commit, and push to `origin`. New commits provide file-level history and diffs. `SHA256SUMS` identifies the prepared source snapshot; refresh it when preparing a new distribution.

## Quick reconstruction from the adopted computation records

Run from the repository root. Python 3.12 and R 4.3.1 were used for validation. The Python environment was freshly installed; the R library was checked against the recorded versions.

```sh
python3.12 -m venv .venv
.venv/bin/python -m pip install -r environment/requirements.txt
Rscript --vanilla scripts/install-dependencies.R
.venv/bin/python -B scripts/rebuild.py
Rscript --vanilla scripts/render.R
.venv/bin/python -B scripts/population.py
```

`rebuild.py` recomputes counts, denominators, standard errors, exact or Wilson confidence intervals, and paired differences from the saved replication-level records. It checks every displayed numerical table value against the frozen manuscript table at its printed precision. It also verifies the numerical statements specified in RESULT_MAP. `render.R` builds all 20 figures and exports the data used in each plot. `population.py` regenerates the E.5 background samples and computes the population approximations and the E.8 scan norm. It verifies each generated dataset against its recorded SHA-256.

Outputs are in `outputs/rebuilt/` and `outputs/figures/`. Repeating these reconstruction commands replaces those generated files. `reference/` contains adopted calculation records. `expected/` contains the manuscript's table fragments and figure assets solely for comparison. The renderer does not copy these PDFs as its output. Reconstructed plot styling can differ from the manuscript; plotted numbers and conventions are preserved.

## R installation

`environment/R-dependencies.csv` records 149 direct and transitive package versions, dependency edges, and upstream license fields. The direct packages include MASS, Matrix, glmnet, randomForest, ggplot2, scales, digest, hdi, SIHR, lpSolve, gmp, Rmpfr, jsonlite, igraph, huge, nnet, mvtnorm, CVXR, rags2ridges, glasso, bestNormalize, VGAM, and ECOSolveR. SIHR 2.0.1 is required by the checked CGM diagnostic adapter.

For an existing library, the command above checks all pins and fails on a mismatch. For a separate library, obtain the exact source archives, then install them explicitly:

```sh
Rscript --vanilla scripts/install-dependencies.R --download --source-dir environment/sources
Rscript --vanilla scripts/install-dependencies.R --install --source-dir environment/sources --library environment/R-library
R_LIBS_USER="$PWD/environment/R-library" Rscript --vanilla scripts/install-dependencies.R
```

The downloader requests the pinned CRAN sources and Bioconductor 3.17 sources for BiocGenerics, graph, and RBGL. It stops if a source is unavailable. Network download and a fresh build of the complete R dependency tree were **not** performed during this validation. R packages are not bundled. Source installation requires a C/C++/Fortran toolchain and relevant system libraries, including GMP and MPFR. Exact pins alone do not fix the operating system, BLAS, compiler, or numerical solver. Retain the generated session information when comparing results on another machine. Set `R_LIBS_USER` for subsequent analysis commands when using the separate library.

## Run an experiment

The dispatcher prints its command and creates a new directory with a command record and log. No experiment installs packages or downloads data. It does not overwrite the frozen records.

```sh
.venv/bin/python -B scripts/run.py gaussian-size --trial
.venv/bin/python -B scripts/run.py E03 --trial
.venv/bin/python -B scripts/run.py nongaussian --trial
.venv/bin/python -B scripts/run.py matched --trial
.venv/bin/python -B scripts/run.py maxcoeft --trial
.venv/bin/python -B scripts/run.py connected --trial
```

Remove `--trial` to run the full protocol. Trial outputs have reduced repetitions or copy counts where the source specifies them; they are excluded from manuscript evidence. The MaxCoefT pilot retains the full 999-copy protocol for 28 datasets. The connected trial runs an independent implementation check on one full dataset. Direction and independent-target centering use their fixed full protocols and reject `--trial`. Stock goodness of fit likewise has a fixed 1,000-copy protocol.

The full matched comparison has two stages:

```sh
.venv/bin/python -B scripts/run.py matched --workers 4
.venv/bin/python -B scripts/run.py maxcoeft --workers 4 --matched-raw experiments/d5-joint-comparison-20261007/runs/RUN_NAME/raw.jsonl
```

Replace `RUN_NAME` with the directory printed by the first command. Without `--matched-raw`, MaxCoefT reuses the bundled adopted baseline p-values, as in the original analysis. The supplied input must have exactly the formal 14-by-500 grid. Both stages regenerate the same datasets and check their identity. The new parameter changes the input location only.

## Applications

The original sources are [Yahoo Finance](https://finance.yahoo.com/) and [NCBI GEO GSE2034](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE2034). The historical stock and clinical files are available from the [fixed GGMtest revision](https://github.com/xiaotong0201/GGMtest/tree/35a3647910792d39e27986b05e07b04b0a6f6461/GGM-Test-Version0/Application). The following layout is required:

| Files below `data/raw/` | Obtain from |
|---|---|
| `Stock/DataForExp.RData`, `Stock/Pre-X.RData`, `Stock/list_stocks.csv` | `Stock/` in the fixed GGMtest revision above |
| `Breast/GSE2034-PI.csv`, `Breast/ImpVar.csv` | `Breast/` in the same revision |
| `Breast/GSE2034_series_matrix.txt` | [GEO processed matrix](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE2nnn/GSE2034/matrix/GSE2034_series_matrix.txt.gz), decompressed |
| `Stock/DIA.json` | Exact complete Yahoo snapshot downloaded on 23 October 2026; included in the review package. Its public archive URL will be added when the author chooses the release location. |

For this repository, the explicit downloader obtains the six public files and accepts the exact DIA snapshot as a local input. It checks the size and SHA-256 of every file against `data/required-inputs.json`. The historical inputs in GGMtest and the current GEO matrix were verified against those checksums; a new Yahoo query is not substituted for the adopted DIA snapshot.

```sh
python3 scripts/fetch-data.py data/raw --dia-source /path/to/DIA.json
Rscript --vanilla scripts/prepare-data.R raw data/raw
```

When a public URL for the exact DIA file exists, `--dia-url URL` can replace `--dia-source`. To check files already present without any download, run `python3 scripts/fetch-data.py data/raw --check`. The **review package already contains all seven files**: its README gives the direct preparation command and requires no data download. Installation of the R/Python dependencies remains a separate step in either package.

Preparation writes only inside this package and refuses to overwrite existing prepared objects. It reconstructs the 92-block stock analysis and the 209-patient ER-positive breast-cancer analysis. [DATA.md](DATA.md) records the transformation and endpoint details. Then run:

```sh
.venv/bin/python -B scripts/run.py bcr-observed
.venv/bin/python -B scripts/run.py bcr-null --trial
.venv/bin/python -B scripts/run.py stock --scenario stock_three_stat_replay_T1 --trial
.venv/bin/python -B scripts/run.py stock-graph
.venv/bin/python -B scripts/run.py stock-gof
```

The BCR runner reproduces the adopted original random-number sequence, then applies the certified likelihood calculation and the random-forest empty-level repair. It does not rerun all methods under a changed random-number sequence. Its trial includes null replicates 1 and 74, including a repaired RF case. Full null mode uses all 400 replicates. The stock runner uses the complete DIA response and 1,000-tree forests. See RESULT_MAP for all nine stock scenarios. Stock scenarios run serially by default; the underlying engine can be used for parallel campaigns after checking local resources.

## Rebuild from new completed runs

`adopt-run.py` imports a specified completed run into a **separate** rebuild input directory. It initially copies the bundled reference and records every replacement in `INPUTS.json`. Until every relevant result has been replaced, this is explicitly a mixture of new and adopted computations. Trial inputs fail the formal replication-count checks.

```sh
.venv/bin/python scripts/adopt-run.py direction outputs/DIRECTION_RUN/gcrt-direction-centering-m999-20261008
REPLICATION_REFERENCE=outputs/fresh-reference .venv/bin/python -B scripts/rebuild.py
REPLICATION_REFERENCE=outputs/fresh-reference Rscript --vanilla scripts/render.R
```

Use the result-containing directory, not merely the dispatcher log directory. Other import kinds are listed by `--help`. For connected tests, use the `controls/results` directory containing all five methods; for non-Gaussian MVN tests, use `results/formal`. Import both `matched` and `maxcoeft` to replace both stages of Table S2. Rebuild checks retain the frozen manuscript comparison, so a changed computation is reported instead of silently relabelled as the submitted result.

## Cost and provenance

The full graphical and application suites are expensive. The historical complete stock campaign required about 49.5 process-hours (about 4.3 wall-hours with 12 workers). The current validation did not rerun that campaign. A single observed stock replay retains all 1,000 copies and 1,000 trees and can take several minutes. The direction and centering examples are much smaller; their complete calculations were rerun. Use the trial records in VALIDATION before scheduling full calculations.

`provenance/original-sources.json` identifies the local source of each copied file and its original hash. `provenance/portability.patch` records adaptations relative to those sources. `provenance/manuscript-inventory.json` identifies the actual loaded manuscript files, assets, labels, and hashes. Dates and old directory names identify provenance; they are not additional manuscript versions to execute. `SHA256SUMS` identifies the prepared source snapshot. Subsequent changes are tracked by Git commits. The public DIA archive location and third-party redistribution status remain identified in THIRD_PARTY.md; neither changes the numerical implementation supplied for review.
