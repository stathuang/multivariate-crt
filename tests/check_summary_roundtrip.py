#!/usr/bin/env python3
"""Check that 17-digit simulation CSVs still produce the paper's power tables."""
import csv
import importlib.util
from pathlib import Path
import shutil
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def main():
    spec = importlib.util.spec_from_file_location("summarize", ROOT / "scripts/summarize.py")
    summarize = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(summarize)
    reference = ROOT / "results/reference/05_graphical_power"
    with tempfile.TemporaryDirectory(prefix="crt-csv-regression-") as directory:
        temporary = Path(directory)
        for part in ["direction", "independent", "connected"]:
            (temporary / part).mkdir()
            shutil.copy2(reference / part / "tests.csv", temporary / part / "tests.csv")
        # Recreate the r column emitted by R's sprintf("%.17g") writer.
        with (reference / "connected/tests.csv").open(newline="") as source:
            reader = csv.DictReader(source)
            rows = list(reader)
            fieldnames = reader.fieldnames
        for row in rows:
            row["r"] = format(float(row["r"]), ".17g")
        with (temporary / "connected/tests.csv").open("w", newline="") as output:
            writer = csv.DictWriter(output, fieldnames=fieldnames, quoting=csv.QUOTE_ALL)
            writer.writeheader()
            writer.writerows(rows)
        summarize.INPUT = temporary
        summarize.OUT = temporary / "tables"
        summarize.OUT.mkdir()
        summarize.VERIFY = True
        summarize.summarize_graphical_power()
        assert len(summarize.TABLES) == 5
        assert sum(check["items"] for check in summarize.CHECKS) == 348
        assert all(check["status"] == "passed" for check in summarize.CHECKS)
        print("PASS: Tables S4–S8, 348 displayed values, with 17-digit CSV parameters.")


if __name__ == "__main__":
    main()
