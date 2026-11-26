# Verification

The expected table fragments are comparison fixtures from the manuscript. The analysis does not use their values to compute results.

```sh
.venv/bin/python scripts/summarize.py all
Rscript --vanilla scripts/figures.R all
.venv/bin/python tests/check_replays.py
```

The first command recomputes all 12 numerical tables from saved replication records and checks every displayed value at the paper's printed precision. It also checks counts, paired comparisons, and the saved likelihood certificates. The second regenerates the 20 figure assets and exports their plotted data.

`check_replays.py` compares available quick-run outputs against the corresponding saved replications. It reports exactly which calculations were present and checked. It rejects empty result files and removes the previous verification verdict before each check. It tests p-values, status, dataset identity where available, and selected statistic values. A quick run checks reproducibility of those calculations; it is not a fresh full Monte Carlo run.

The CSV regression check recreates the connected-target results in the 17-digit format used by new calculations, then verifies Tables S4–S8 against the paper. It uses saved records and does not rerun simulations:

```sh
.venv/bin/python tests/check_summary_roundtrip.py
```

The empty-result regression check confirms that empty CSV and JSONL files fail, valid records still pass, and a failed check cannot leave a stale passing verdict:

```sh
.venv/bin/python tests/check_replays_empty.py
```
