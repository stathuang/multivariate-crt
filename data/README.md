# Application data

Simulation studies generate their own data. The stock and breast-cancer studies use seven input files, listed with exact checksums in `required-inputs.json`.

| File under `data/raw/` | Source |
|---|---|
| `Stock/DataForExp.RData`, `Stock/Pre-X.RData`, `Stock/list_stocks.csv` | Stock files in the [fixed GGMtest revision](https://github.com/xiaotong0201/GGMtest/tree/35a3647910792d39e27986b05e07b04b0a6f6461/GGM-Test-Version0/Application) |
| `Breast/GSE2034-PI.csv`, `Breast/ImpVar.csv` | Breast files in the same revision |
| `Breast/GSE2034_series_matrix.txt` | [NCBI GEO processed expression matrix](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE2nnn/GSE2034/matrix/GSE2034_series_matrix.txt.gz), decompressed |
| `Stock/DIA.json` | Complete Yahoo Finance snapshot acquired on 23 October 2026, included at [snapshots/DIA.json](snapshots/DIA.json) |

Download and verify the six public files, copy the included DIA snapshot, and prepare both applications:

```sh
python3 scripts/fetch_data.py data/raw
Rscript --vanilla scripts/prepare_data.R raw data/raw
Rscript --vanilla scripts/check_stock_transform.R data/raw/Stock
```

Each application can also be prepared independently. For breast cancer:

```sh
python3 scripts/fetch_data.py data/raw --dataset bcr
Rscript --vanilla scripts/prepare_data.R raw data/raw bcr
```

For stocks, use `--dataset stock` and the final argument `stock` instead. Omitting the dataset selection retains the default of preparing both applications. The breast-cancer route does not require any stock files, and the stock route does not require any breast-cancer files.

The downloader checks sizes and SHA-256 hashes. Add `--check` to verify existing files without copying or downloading anything. Optional `--dia-source /path/to/DIA.json` or `--dia-url HTTPS_URL` overrides the bundled location; the supplied file must still match the recorded hash. The downloader does not substitute a current Yahoo query for the saved snapshot. Only the fixed DIA snapshot is included in Git; `data/raw/` is a local input directory.

Preparation writes the selected objects, `data/processed/stock.rds` and/or `data/processed/bcr_erplus60.rds`, and records their hashes in `data/processed/IMPORT.csv`. Separate runs retain both manifest entries. It refuses to overwrite existing processed objects. The stock object contains 92 complete five-day windows, 103 stocks, the fitted graph, two target sets, and fitted null models. The transform check independently reconstructs the stock matrix from saved prices.

The breast-cancer preparation uses the processed GEO expression matrix. It selects the first 60 probes in the original ordered marker list, applies order normalisation across all 286 patients, and then retains the 209 ER-positive patients. The fixed analysis has three background probes, 18 targets, and response counts 152, 52, and 5. The marker list comes from Wang et al. (2005), Lancet 365:671–679, Table 3, [doi:10.1016/S0140-6736(05)17947-1](https://doi.org/10.1016/S0140-6736(05)17947-1).

The clinical input contains relapse/follow-up time and a brain-relapse flag, without a separate brain-relapse date. Category 2 is relapse within 48 months with that flag. The earliest conversion of clinical records to this CSV has not been reconstructed. Target selection depends on the response; the null simulations keep the selected target and background sets fixed.
