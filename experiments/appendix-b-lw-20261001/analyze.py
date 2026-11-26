"""Check completed exports and summarize paired outcomes; standard-library only.

Usage: python3 analyze.py E01_OUTPUT_DIRECTORY E02_OUTPUT_DIRECTORY
"""
import csv
import json
import math
from pathlib import Path
import statistics
import sys


def read_csv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle))


def check_run(directory, power):
    cells = {}
    total = 0
    repeats, copies = (400, 199) if power else (5000, 999)
    with (directory / "results.csv").open() as handle:
        for row in csv.DictReader(handle):
            key = (row["setting"], row["theta"] if power else row["statistic"], row["method"])
            rep = int(row["replication"])
            cell = cells.setdefault(key, {})
            assert rep not in cell, ("Duplicate", key, rep)
            assert row["status"] == "ok" and not row["error"] and not row["warnings"], row
            p = float(row["p_value"])
            assert 1 / (copies + 1) <= p <= 1, row
            assert abs(p * (copies + 1) - round(p * (copies + 1))) < 1e-9, row
            cell[rep] = int(p <= 0.05)
            total += 1
    assert len(cells) == (105 if power else 45)
    for cell in cells.values():
        assert set(cell) == set(range(1, repeats + 1))
    for row in read_csv(directory / "counts.csv"):
        key = (row["setting"], row["theta"] if power else row["statistic"], row["method"])
        values = list(cells[key].values())
        rate = statistics.mean(values)
        reported = float(row["rejection_rate"] if power else row["empirical_size"])
        assert abs(reported - rate) < 1e-15
        assert abs(float(row["se"]) - math.sqrt(rate * (1 - rate) / repeats)) < 1e-15
        assert row["complete"] == "TRUE" and int(row["failed"]) == int(row["missing"]) == 0
    return cells, total


def wilson(values):
    n, estimate, z = len(values), statistics.mean(values), 1.959963984540054
    divisor = 1 + z * z / n
    center = (estimate + z * z / (2 * n)) / divisor
    half = z * math.sqrt(estimate * (1 - estimate) / n + z * z / (4 * n * n)) / divisor
    return [max(0.0, center - half), min(1.0, center + half)]


def main():
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    size_dir, power_dir = map(Path, sys.argv[1:])
    size, size_count = check_run(size_dir, False)
    power, power_count = check_run(power_dir, True)
    root = Path(__file__).resolve().parent
    comparisons = []
    for (setting, theta, method), reps in power.items():
        if method != "LW plug-in CRT (D.10)":
            continue
        for baseline in ["Oracle CRT", "MVN-CRT"]:
            other = power[setting, theta, baseline]
            differences = [reps[i] - other[i] for i in sorted(reps)]
            comparisons.append(dict(setting=setting, theta=float(theta), comparison="LW minus " + baseline,
                                    n=len(differences), difference=statistics.mean(differences),
                                    paired_mc_se=statistics.stdev(differences) / math.sqrt(len(differences))))
    with (root / "paired-power-comparisons.csv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(comparisons[0]))
        writer.writeheader()
        writer.writerows(comparisons)
    null = []
    for (setting, statistic, method), reps in size.items():
        if statistic == "ridge gain":
            values = list(reps.values())
            null.append(dict(setting=setting, method=method, n=len(values), rejections=sum(values),
                             rate=statistics.mean(values), wilson_95=wilson(values)))
    result = dict(size_rows=size_count, power_rows=power_count,
                  checks="Unique complete replicates, p-value grids, errors/warnings, counts and MC SEs passed.",
                  null_ridge=null)
    (root / "analysis.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
