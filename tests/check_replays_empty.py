#!/usr/bin/env python3
"""Exercise the replay checker on empty and valid result files in temporary copies."""
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
CSV = Path("02_graphical_crt/FigureS1b/results.csv")
JSONL = Path("03_joint_tests/ridge.jsonl")


def main():
    with (ROOT / "results/reference" / CSV).open() as source:
        header, row = next(source), next(source)
    with (ROOT / "results/reference" / JSONL).open() as source:
        record = next(source)
    cases = [
        ("header-only CSV", {CSV: header}, False),
        ("zero-byte CSV", {CSV: ""}, False),
        ("empty JSONL", {JSONL: ""}, False),
        ("empty CSV alongside valid JSONL", {CSV: header, JSONL: record}, False),
        ("empty JSONL alongside valid CSV", {CSV: header + row, JSONL: ""}, False),
        ("valid CSV", {CSV: header + row}, True),
        ("valid JSONL", {JSONL: record}, True),
        ("valid CSV and JSONL", {CSV: header + row, JSONL: record}, True),
    ]
    with tempfile.TemporaryDirectory(prefix="crt-empty-replay-") as directory:
        for index, (name, outputs, passes) in enumerate(cases):
            root = Path(directory) / str(index)
            (root / "tests").mkdir(parents=True)
            shutil.copy2(ROOT / "tests/check_replays.py", root / "tests/check_replays.py")
            for path, text in {CSV: header + row, JSONL: record}.items():
                target = root / "results/reference" / path
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_text(text)
            for path, text in outputs.items():
                target = root / "results/quick" / path
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_text(text)
            verdict = root / "results/quick/verification.json"
            verdict.write_text(json.dumps({"status": "passed", "checks": ["stale"]}))
            result = subprocess.run(
                [sys.executable, str(root / "tests/check_replays.py")],
                cwd=root, capture_output=True, text=True, timeout=30,
            )
            assert (result.returncode == 0) == passes, (name, result.stdout, result.stderr)
            if passes:
                checked = json.loads(verdict.read_text())
                assert checked["status"] == "passed"
                assert len(checked["checks"]) == len(outputs)
                assert all(c.get("rows", c.get("datasets", 0)) > 0 for c in checked["checks"])
            else:
                assert "PASS:" not in result.stdout, (name, result.stdout)
                assert not verdict.exists(), (name, "stale passing verdict survived")
                assert "Empty quick results:" in result.stderr or "EmptyDataError" in result.stderr, result.stderr
            print(f"PASS: {name}")


if __name__ == "__main__":
    main()
