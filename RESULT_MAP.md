# Results and executable commands

This map covers the main text and supplements loaded by v25 `Main.tex`. Figure and table numbers were aligned with the manuscript compiled on 26 November 2026; supplementary figures and tables use the S prefix. The adopted calculations and the frozen provenance inventory remain those verified on 20 November 2026. Source labels identify the same results across these versions. Commands below follow `python -B scripts/run.py` unless a different script is stated. `python` means the configured Python 3.12 environment. All plot and table summaries can be rebuilt without application data by `scripts/rebuild.py` and `scripts/render.R`.

## Every submitted figure and table

| Result | Source label / location | Full computation | Adopted result input / reconstructed output |
|---|---|---|---|
| Table 1 | `tab:joint-crt-benchmark-size`, SimulationStudies | `gaussian-size` | `reference/gaussian/E01/results.csv`; `outputs/rebuilt/gaussian-size.csv` |
| Figure 1 | `fig:oracle-plugin-power-curve`, SimulationStudies | `gaussian-power` | `reference/gaussian/E02/results.csv`; `outputs/rebuilt/gaussian-power.csv` |
| Figure 2 | `fig: stock noise`, RealDataIllustrations | `stock --scenario stock_noise_05`, then `_10`, `_15`, `_20`, `_25` | `reference/stock/stock_noise_*.csv`; `outputs/rebuilt/stock-noise.csv` |
| Table S1 | `tab:mscrt-statistic-catalogue`, AppendixSimulationDetails | Definitions, not a Monte Carlo result; code map below | Linear, binary, and RF statistic modules |
| Figures S1–S2 | `fig: app linear regression`, `fig: app nonlinear regression` | `E03`, `E04`; `E05`, `E06` | `reference/graphical/main_{l,h}_l_{w,m}.csv`; `outputs/rebuilt/graphical.csv` |
| Figures S3–S4 | `fig: app logistic regression`, `fig: app nonlinear binary` | `E07`, `E08`; `E09`, `E10` | `reference/graphical/main_{l,h}_gl_{w,m}.csv`; `outputs/rebuilt/graphical.csv` |
| Table S2 | `tab:matched-direct-comparison`, JointAndCoordinatewiseComparison | `matched`, then `maxcoeft --matched-raw PATH` | `experiments/d5-joint-comparison-20261007/results/complete/raw.jsonl` plus `reference/maxcoeft/raw.jsonl`; `outputs/rebuilt/matched.csv` |
| Figure S5 | `fig: linear regression-univariate`, AppendixAdditionalSimulations | `E12`, `E13` | `reference/graphical/univariate_{l,h}_l_w.csv` |
| Figure S6 | `fig: varying graphsize`, AppendixAdditionalSimulations | `E11` | `reference/graphical/supergraph_h_l_w.csv` |
| Figure S7 | `fig: linear regression-t3`, AppendixAdditionalSimulations | `E14`, `E15` | `reference/graphical/t3_{l,h}_l_w.csv` |
| Figure S8, Table S3 | `fig:mvn-nongaussian-rejection-rates`, `tab:mvn-nongaussian-size`, AppendixAdditionalSimulations | `nongaussian` | `reference/nongaussian/cell-*.csv`; `outputs/rebuilt/nongaussian.csv`, `nongaussian-size.tex` |
| Table S4 | `tab:direction-mechanism`, AppendixGCRTComparisons | `direction` | `reference/direction/direction-tests.csv`; `outputs/rebuilt/direction.csv` |
| Figure S9, Tables S5–S6 | `fig:centering-comparison`, `tab:centering-rates`, `tab:centering-paired`, AppendixGCRTComparisons | `centering` | `reference/centering/tests.csv`; `outputs/rebuilt/centering*.csv` |
| Figure S10, Tables S7–S8 | `fig:connected-centering`, `tab:connected-centering-rates`, `tab:connected-centering-paired`, AppendixConnectedCentering | `connected` | `reference/connected/tests.csv`; `outputs/rebuilt/connected*.csv` |
| Figure S11 | `fig: sp100_G`, AppendixRealDataIllustrations | `stock-graph` | Public graph summary `reference/stock/graph.rds`; `outputs/figures/Application/Stock/stocks_G.pdf` |
| Table S9 | `tab: GoF stock`, AppendixRealDataIllustrations | `stock-gof` | `reference/stock/Stock_GoF.csv`, `gof.rds`; `outputs/rebuilt/Stock_GoF.tex` |
| Table S10 | `tab: fund error check`, AppendixRealDataIllustrations | `stock --scenario stock_typeI_T1`, then `stock_typeI_T2` | `reference/stock/stock_typeI_*.csv`; `outputs/rebuilt/stock-null.csv` |
| Table S11 | `tab: fund p-value`, AppendixRealDataIllustrations | `stock --scenario stock_three_stat_replay_T1`, then `stock_three_stat_replay_T2` | `reference/stock/stock_three_stat_replay_*.csv`; `outputs/rebuilt/stock-observed.csv` |
| Table S12, Figure S12 | `tab:appendix BCR type 1`, `fig: BCR pvalue lrt`, AppendixRealDataIllustrations | `bcr-null` | `reference/bcr/UNRESTRICTED_RESULTS.csv`; `outputs/rebuilt/bcr-null.csv` |
| BCR observed p-values in Section 8.2 | RealDataIllustrations | `bcr-observed` | `reference/bcr/observed.csv`; `outputs/rebuilt/bcr-observed.csv` |

The graph, histogram, and multipanel displays explain why 14 figure environments require 20 PDF files. `outputs/figures/MANIFEST.csv` lists every regenerated asset and the exact plotted CSV beside it. Numerical `.tex` outputs are compact rows for inspection; `expected/` preserves the manuscript formatting separately. The computed table columns follow the corresponding manuscript table order.

## Protocols and comparison methods

- **Gaussian joint CRT benchmark:** E01 uses three `(n,p,t)` settings `(60,45,10)`, `(80,50,10)`, and `(80,60,40)`, 5,000 repetitions, and 999 randomised copies. E02 uses 400 repetitions per signal setting and 199 copies. Methods are the oracle, maximum-likelihood plug-in, degrees-of-freedom plug-in, MVN-CRT, and Ledoit–Wolf plug-in CRT. Both use `p <= 0.05`; their SE is the binomial population formula. The later Ledoit–Wolf extension is included here, rather than the earlier four-method export.
- **Main graphical simulations:** E03/E04 are low/high dimensional linear responses; E05/E06 nonlinear continuous responses; E07/E08 logistic responses; E09/E10 nonlinear binary responses. E03–E08, E11–E15 use 400 repetitions; E09/E10 use 800. Copy count is 100 and sweep count is 3. Comparators include the F-test, numerically checked chi-squared likelihood ratio test, de-sparsified Lasso (`hdi`), CGM (`SIHR`), and Lasso/RF dCRT. All original random-number ordering and per-method warning/error records are retained. These summaries use `p < 0.05` and sample-indicator SE. Two CGM outputs in the high-dimensional logistic simulation at theta 3.5 failed; that method's denominator there is 398, with the two failed attempts preserved. No replacement p-value is supplied.
- **Matched comparison:** two connected precision designs, 14 cells, 500 independent datasets per cell, replication IDs 10000–10499, `n=80`, `p=120`, eight targets, `M=999`, and `L=3`. Background ridge penalty is `0.2n`; target penalty is one. Both designs use the same standardized group signal. The seven displayed tests are coordinatewise G-CRT with Bonferroni correction, joint G-CRT with ridge, MaxCoefT, and their two-test Bonferroni combination, and the three corresponding oracle joint tests. Lasso background cross-validation for MaxCoefT has its own recorded R seed. Formal baseline p-values and new MaxCoefT p-values share verified data hashes. Historical unreported exploratory statistics in the base runner are retained as computation provenance, but do not create additional submitted results.
- **Non-Gaussian MVN-CRT:** Gaussian, standardized t5, and standardized t3 errors; sample sizes 50, 100, 200, and 400; 12 background and eight target variables; 1,000 repetitions per cell; `M=199`; nominal levels 0.005 through 0.100. A common scale mixture applies within each multivariate t row. The explicit seed stream begins after the separately recorded 60 pilot datasets. Confidence intervals are Wilson intervals. The pilot and formal outputs remain distinct.
- **Direction:** 500 datasets at each of `n=40,120,360`, three responses on each dataset, `M=999`, `L=1`; 4,500 test rows. Frozen simulated prefixes verify the original random streams.
- **Independent-target centering:** 100 datasets at each of three sample sizes, four responses and four statistics, `M=19`, `L=1`; 4,800 test rows. The manuscript uses this 19-copy study. The separate 999-copy centering experiment is not substituted. Pointwise intervals are exact binomial intervals; differences and their Monte Carlo SEs use paired rejection indicators.
- **Connected-target centering:** `n=80`, eight dependent targets, 112 background variables, `M=99`, `L=1`; 2,500 shared datasets and eight response settings. The calculation stores sufficient-statistic checks, score vectors, and Gram matrices. An additional deterministic pass evaluates the two uncentred fixed-weight controls on those same stored quantities. Tables use all five methods (100,000 p-values); the figure uses the original three methods. This distinction is intentional.
- **Stocks:** 92 complete five-trading-day windows for 103 stocks, with two fixed target groups. The accepted X matrix and graph are unchanged; the DIA response uses the complete 461-date snapshot. Observed tests use 1,000 copies; fitted-null and noise calculations use 400 copies and 400 repetitions. All RF branches use 1,000 trees. The nine scenario names are the two `stock_three_stat_replay_*`, two `stock_typeI_*`, and five `stock_noise_*` entries in the table. Saved fitted-null models, seeds, folds, and tree counts are carried through the frozen engine configurations.
- **Breast cancer:** 209 ER-positive patients, the first 60 probes in the original 76-probe list in their exact list order, three selected background and 18 selected target probes. Transformation is fitted on all 286 patients before ER-positive restriction. The endpoint and selected hypothesis remain exploratory; see DATA. Four hundred unrestricted fitted-null datasets are retained. The RF repair changes only the failed RF results for replicates 74, 306, 375, 383, and 390; it does not condition on every class being observed. Certified likelihood values replace the earlier numerical LRT calculation. Existing dCRT p-values retain the original RNG path.

## Numerical statements beyond tables

`rebuild.py` and `population.py` verify or expose the following calculations. Percentages are converted from the unrounded records before comparison.

| Manuscript statement | Calculation / output |
|---|---|
| Section 7.1 oracle/MVN maximum power difference 0.023 | Maximum absolute cell difference 0.0225 in `gaussian-power.csv` |
| E.5 maximum G-CRT/oracle difference 4.4 percentage points | Seven-method matched grid; maximum paired-method cell difference 0.044 |
| E.5 dense theta 0.75 power comparison, 39.6/82.8% in A and 75.8/10.0% in B | Table S2 values and `matched.csv` |
| E.5 conditional standardized effects 3.514–3.722 (A), 0.934–1.198 (B), group variance 0.5625 | Regenerated precision matrices and response coefficients; `population.json` |
| E.5 leading gain differences 0.0408 (A), 0.0945 (B) | Exact original 500 background samples per design; `population.json`; unrounded values 0.0407600449497 and 0.0944721688355 |
| E.5 combined-test changes -5.8/+5.8 points, paired SE 1.1/1.4 points | Paired rejection indicators in the saved matched records; `matched-paired.csv` |
| Non-Gaussian sizes 4.8–6.0%, Gaussian 4.4–5.5% | 12,000 p-values, 240 level-specific Wilson intervals; `nongaussian.csv` |
| Direction, independent-centering, and connected-centering claims | Full cell means/rates, paired differences, and exact intervals in `direction.csv`, `centering*.csv`, and `connected*.csv` |
| E.8 scan norm about 0.0598 | Precision-induced operator norm 0.0598150221680; `population.json` |
| Stock graph maximum/median degree 30/19 and edge density about 18.4% | `scripts/render.R`, saved graph; `stock-graph-summary.csv` |
| Stock transformation exponents 1.4 (MMM), 1.5 (MO) | `Rscript scripts/check-stock-transform.R INPUT/Stock`; complete transformed X matrix compared with the preserved snapshot |
| Stock noise sigma 0.1: 94.3% vs 72.8% | 377/400 vs 291/400; `stock-noise.csv` |
| Stock observed and null results, including 0.75–7.25% CRT rates and 26/15% de-Lasso rates | Tables S10–S11, `stock-null.csv`, and `stock-observed.csv` |
| BCR observed 9.58e-6, 0.018, 0.697, 0.407 | `bcr-observed.csv` in LRT, RF, dCRT(L1), dCRT(RF) order |
| BCR log-likelihood bound gaps at most 1e-5 | `reference/bcr/likelihood-certificates.csv`: 401 saved bound summaries and source hashes; `bcr-adopted.R` recomputes the likelihood fits. Full patient-derived certificate arrays remain outside the public bundle. |
| BCR null 40%, 7.25%, 3.25%, 3.25% | 160, 29, 13, 13 rejections among 400; `bcr-null.csv` |

Analytic constants in theorem statements are mathematical definitions or bounds, not simulation output. The original question of proof correctness is outside this packaging task. Numerical examples used to explain the theory are included above.

## Statistic and method definitions in code

Paths below are relative to `experiments/readable-experiments/R/`, except where specified.

| Table S1 name / method | Source |
|---|---|
| LM-SST, LM-SSR, and their Lasso/distillation/residual variants | `linear_tests.R`, `gaussian_tests.R` |
| RF, RF-D, RF-RR | `gaussian_tests.R`, `application_tests.R` |
| GLM-Dev, GLM-L1-D, GLM-L1-R-SST | `binary_tests.R`, `binary_likelihood.R` |
| G-CRT hub-and-spoke sampler | `sampler.R`; stock engine `R/core.R` |
| Multinomial LRT with numerical verification | `multinomial_likelihood.R`, `bcr_tests.R` |
| dCRT (L1), dCRT (RF) | `dcrt.R`; frozen original/application baseline modules |
| De-Lasso and CGM diagnostics | `comparators.R`, `hdi_cache.R`, `cgm_diagnostics.R` |
| Matched ridge / coordinatewise test | `experiments/d5-joint-comparison-20261007/src/mscrt/core.py` |
| MaxCoefT and fixed Lasso background | `experiments/e5-max-coefficient-t-20261103/src/statistic.py`, `background.R` |
| MVN-CRT and Ledoit–Wolf plug-in | `experiments/appendix-b-lw-20261001/R/` |

Historical paths in manifests identify their source runs. Executable adapters use this package's own files. Never execute a command embedded in a historical manifest as a setup instruction.
