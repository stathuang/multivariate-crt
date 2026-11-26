# dCRT source, historical changes, and reconstruction

The distilled conditional randomization test (dCRT) comparator follows the code attributed to Molei Liu in [Distillation-CRT](https://github.com/moleibobliu/Distillation-CRT), associated with Liu, Katsevich, Janson, and Ramdas, *Fast and Powerful Conditional Randomization Testing via Distillation*. This document distinguishes upstream updates, the historical GGMtest adaptation, and the restricted implementations used by this manuscript.

## What changed upstream, and when?

The public commit history checked on 20 November 2026 contains 23 upstream commits. Its most recent commit is dated **2 January 2024**. The two latest changes are:

- [24 June 2023](https://github.com/moleibobliu/Distillation-CRT/commit/7383560eab7d0024794dc01220fd0f5eddb01742): an empty-selection return was added, with an unrelated comment-encoding change.
- [2 January 2024](https://github.com/moleibobliu/Distillation-CRT/commit/dddc8cc4027c236233cb9f0afb65290fb622be15): `doParallel` and `core_num` options, a parallel implementation, and an example were added; the serial implementation was reorganised.

The adopted historical file was already present in [GGMtest on 21 December 2023](https://github.com/xiaotong0201/GGMtest/blob/2b9495b77b3ff4d879e86d8d83d99595ca97b44c/GGM-Test-Version0/CRT/dCRT_function.R). It matches the later GGMtest file and our preserved historical source, apart from the final newline. Thus, the January 2024 upstream update followed this documented use in GGMtest. There is no evidence of an upstream update during the 2026 revision work. This establishes a documented-use date, not the date the research project began.

The random-forest difference described below is already present when comparing the adopted file with upstream's **June 2020** source. It was not introduced by the January 2024 parallel update. We do not know who originally made each local edit or its exact date; the file comparison establishes the edits themselves.

## Historical adaptation relative to the 2021 upstream functions

| Change | Adopted behaviour and scope |
|---|---|
| Rename `dCRT` to `d_CRT` | Separates the upstream-style function from the project's group-test wrapper. |
| Define `p = ncol(X)` inside `Creat_condition_gaussian` | Makes the conditional-fit dimension local instead of relying on an outer variable. |
| For the no-interaction RF branch, replace `X[, -indx] - mean_X` by `X[, -indx]` | Fits the response forest to the original nuisance predictors. The target residual remains `X[, indx] - mean_X`. This is the adopted historical behaviour, not a newly introduced change or a claim that every upstream branch is incorrect. |
| Add the empty-selection return and remove progress prints | Affects the general legacy interface. The manuscript uses no screening and nonempty targets. |
| Add `drop = FALSE` in the known-covariance conditional mean | Preserves matrix dimensions; the adopted estimated-covariance branch does not use this path. |

An executable patch is supplied as `provenance/dcrt-2021-to-adopted.patch`. It operates on the two extracted function definitions, avoiding unrelated generators and estimators in the original script. To reconstruct the historical functions from the pinned 2021 source:

```sh
mkdir -p outputs/dcrt-source
curl --fail --location 'https://raw.githubusercontent.com/moleibobliu/Distillation-CRT/ebab7dcffda9d4af4355120952bed70de484a8b2/functions.R' -o outputs/dcrt-source/upstream.R
Rscript --vanilla scripts/extract-dcrt-functions.R outputs/dcrt-source/upstream.R outputs/dcrt-source/upstream-functions.R
patch -o outputs/dcrt-source/adopted-functions.R outputs/dcrt-source/upstream-functions.R provenance/dcrt-2021-to-adopted.patch
Rscript --vanilla scripts/check-dcrt-adaptation.R outputs/dcrt-source/adopted-functions.R
```

These are optional source-reconstruction commands. Ordinary simulations and application analyses use the supplied comparator implementation and do not download or patch dCRT. Use a new output directory when repeating extraction, because the extractor refuses to overwrite files. The extraction and patch were checked with R 4.3.1; another R version can format `deparse()` output differently.

## Restricted implementations supplied for this paper

The package retains only `CDR = 'No'`, `k = 0`, `MC_free = TRUE`, `Sigma_X = NULL`, and `Gen_X = NULL`. Here `CDR` controls screening, `k` controls response interactions, and `MC_free` selects the analytic tail calculation. The last two arguments select the estimated Gaussian conditional distribution.

For each target coordinate, the implementation preserves the historical order: cross-validate and refit the conditional Lasso; fit the response model; compute the residual-product statistic and its Gaussian tail probability. It then combines coordinate p-values by Bonferroni, `min(1, length(target) * min(pvalues))`. No function sets a random seed internally.

| Implementation | Response fits and additions |
|---|---|
| `experiments/readable-experiments/R/dcrt.R`, `dcrt_test()` | Gaussian Lasso, binomial Lasso, or a 100-tree regression forest; optional residual output for the existing oracle comparison. |
| `experiments/bcr-original/R/baselines.R`, `.baseline_dcrt()` | Same restricted calculation; a 100-tree regression forest with the numeric response used by the historical comparator. |
| `experiments/stock-rf1000/engine/R/baselines.R`, `.baseline_dcrt()` | Same calculation with 1,000 trees, as specified for the adopted stock campaign. The forest wrapper records diagnostics without consuming random numbers. |

The comparator's regression forest must not be confused with the separate classification forest used as a breast-cancer G-CRT statistic. Input checks, package-qualified calls, and numerical tail evaluation make the restricted code easier to run without changing the selected protocol.

`scripts/check-dcrt-adaptation.R` compares all three implementations against the historical functions using Gaussian, binomial, and forest response fits, two targets in a specified order, and the same seed. It checks both p-values and the final random-number state. Its nine synthetic cases supplement the application and manuscript-result validations in `VALIDATION.md`; they are not an exhaustive equivalence proof for all unused upstream branches.

## Attribution and distribution

The submission package supplies the adopted implementation and a historical source copy for inspection, so reviewers need not obtain dCRT separately. The GitHub candidate supplies the same implemented comparisons and the adaptation instructions. It has not been uploaded, and the source patch does not itself grant third-party permission. No explicit upstream software licence was found in the checked sources. Preserve attribution and resolve public redistribution separately from providing a runnable review copy.
