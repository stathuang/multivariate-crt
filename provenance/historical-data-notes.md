> Historical source note. Its commands refer to the earlier study tree. Use [DATA.md](../DATA.md) for this distribution and the corrected DIA response.

# Application input data

`public/` contains a factual probe-identifier list and its citation, suitable for inclusion with the code. `external-cache/` contains locally supplied stock and clinical/expression inputs and every derived application matrix/model. **Exclude external-cache from any public source bundle.** Public availability of a source is not a verified redistribution licence for these exact cached files. No patient-level data or market-price derivatives are declared redistributable by this package.

`INPUT_MANIFEST.csv` records SHA-256 hashes and sizes of six imported inputs and three prepared RDS datasets. It contains portable paths. The import source is the author's preserved historical application collection, captured on 2026-09-27; the exact source-directory location is supplied explicitly by the caller and is never embedded in code. Original source files are read-only. Preparation writes only to this package's `data/` directory. Analysis uses only the cached inputs, never a sibling reference directory.

## Explicit local preparation

From the package directory, in R:

```r
source("R/applications.R")
prepare_application_data(source_dir = "/path/to/supplied/Application", data_dir = "data")
```

The supplied directory must contain these exact files:

- `Stock/DataForExp.RData`: transformed 92-by-103 covariate matrix, saved DIA response, graph, two target groups.
- `Stock/Pre-X.RData`: saved transformed matrix, stock sectors, dates and ancillary source snapshot.
- `Stock/list_stocks.csv`: 103 ordered stock symbols and sectors.
- `Breast/GSE2034_series_matrix.txt`: GEO processed expression matrix; no CEL-level reconstruction is claimed.
- `Breast/GSE2034-PI.csv`: supplied clinical table, including brain-relapse flag.
- `Breast/ImpVar.csv`: 76 probe identifiers in Wang Table 3 order.

Existing caches are not overwritten unless the caller explicitly sets `overwrite=TRUE`. Preparation records hashes, preserves warnings, and produces `stock.rds`, `bcr_historical.rds`, and `bcr_erplus60.rds` in `external-cache/prepared/`. RDS files use format version 3 and are portable across supported R platforms. No function downloads data or installs dependencies. A public code bundle without external-cache fails explicitly at application analysis until these inputs are provided; there is no automatic reference-path fallback.

## Source and availability limitations

Stock legacy source code names Yahoo Finance through tidyquant/quantmod, with the covariate interval 2020-09-01 to 2022-06-30 and a DIA response. Its saved transformed input is the accepted reproducibility starting point. The fresh-download path, corporate-action adjustments, historical symbol universe, and exact response/covariate date alignment have not been rebuilt. The stored `Pre-X` covariates are checked exactly against `DataForExp`, and sector-weighted glasso is recomputed. Market-data redistribution terms have not been established; obtain a lawful exact input snapshot separately.

Breast expression source: [NCBI GEO GSE2034](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE2034). The supplied matrix identifies accession GSE2034, PubMed 15721472, public date 2005-02-23, and last series update 2018-08-10. The CSV clinical-table conversion history is not fully documented; its exact bytes are hashed, and fields are checked. The old script's endpoint uses follow-up/relapse month and the brain-relapse flag. A separate brain-event date is not available in this CSV. The exact downloaded source of the CSV and its conversion to the local format remain provenance limitations.

Original marker source: Wang et al. (2005), *Gene-expression profiles to predict distant metastasis of lymph-node-negative primary breast cancer*, Lancet 365:671–679, doi:10.1016/S0140-6736(05)17947-1, Table 3, p.677. [Primary article PDF, institutional mirror](https://delta-apache-vm.cs.tau.ac.il/~nin/Courses/NC06/Wang2005.pdf). Checked 2026-09-27. All 60 ER-positive IDs match the first 60 local `ImpVar.csv` rows in exact order. `public/wang2005_erplus60.csv` records those identifiers only; it contains no expression or patient observations. Original article text/PDF is not repackaged.
