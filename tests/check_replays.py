"""Compare available quick calculations with the saved paper replications."""

from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
QUICK = ROOT / "results/quick"
REFERENCE = ROOT / "results/reference"
checks = []

# A failed check must not leave a previous passing verdict.
(QUICK / "verification.json").unlink(missing_ok=True)


def compare(path, keys, values, status=False):
    fresh = QUICK / path
    if not fresh.exists():
        return
    new = pd.read_csv(fresh).fillna({"case": "null"})
    assert not new.empty, f"Empty quick results: {path}"
    old = pd.read_csv(REFERENCE / path).fillna({"case": "null"})
    assert not new.duplicated(keys).any(), str(path)
    z = new.merge(old, on=keys, validate="one_to_one", suffixes=("_new", "_paper"))
    assert len(z) == len(new), (str(path), len(z), len(new))
    errors = {}
    for name in values:
        a, b = z[name + "_new"], z[name + "_paper"]
        assert np.array_equal(a.isna(), b.isna()), (str(path), name, "missingness")
        ok = a.notna()
        error = float(np.max(np.abs(a[ok] - b[ok]))) if ok.any() else 0.0
        assert np.allclose(a[ok], b[ok], atol=1e-10, rtol=1e-10), (
            str(path),
            name,
            error,
        )
        errors[name] = error
    if status:
        # The numerical LRT certificate is a successful fit; older observed rows used 'ok'.
        normalize = lambda s: s.replace({"numerically_certified": "ok"})
        assert (normalize(z.status_new) == normalize(z.status_paper)).all(), (
            str(path),
            "status",
        )
    checks.append(dict(file=str(path), rows=len(new), maximum_differences=errors))


for kind in ["size", "power"]:
    keys = ["setting", "replication", "method"] + (
        ["statistic"] if kind == "size" else ["theta"]
    )
    compare(Path("01_gaussian_crt") / kind / "results.csv", keys, ["p_value"], True)
for path in sorted((QUICK / "02_graphical_crt").glob("*/results.csv")):
    compare(
        path.relative_to(QUICK),
        ["scenario", "replicate", "theta", "graph", "method"],
        ["pvalue"],
        True,
    )
for path in sorted((QUICK / "04_misspecification/simulation").glob("cell-*.csv")):
    compare(
        path.relative_to(QUICK),
        ["distribution", "n", "repetition"],
        ["p_value", "T0"],
        True,
    )
for part in ["direction", "independent", "connected"]:
    keys = (
        ["replicate", "case"]
        + ([] if part == "direction" else ["method"])
        + ([] if part == "connected" else ["n"])
    )
    compare(
        Path("05_graphical_power") / part / "tests.csv",
        keys,
        ["pvalue"] + ([] if part == "connected" else ["T0", "Tbar", "gap"]),
    )
for path in sorted((QUICK / "06_stock_returns").glob("*/results.csv")):
    compare(
        path.relative_to(QUICK),
        (
            ["method"]
            if path.parent.name == "gof"
            else ["scenario", "replicate", "method"]
        ),
        ["pvalue"],
        path.parent.name != "gof",
    )
for part in ["observed", "null"]:
    compare(
        Path("07_breast_cancer") / part / "results.csv",
        ["replicate", "method"],
        ["pvalue"],
        True,
    )
identity = lambda r: tuple(
    r[k] for k in ["design", "arrangement", "pattern", "theta", "rep"]
)
for name in ["ridge.jsonl", "maxcoeft.jsonl"]:
    path = Path("03_joint_tests") / name
    if not (QUICK / path).exists():
        continue
    new = [json.loads(s) for s in (QUICK / path).read_text().splitlines()]
    assert new, f"Empty quick results: {path}"
    old = {
        identity(r): r
        for r in map(json.loads, (REFERENCE / path).read_text().splitlines())
    }
    assert len({identity(r) for r in new}) == len(new)
    for r in new:
        saved = old[identity(r)]
        assert r["data_sha256"] == saved["data_sha256"], (name, identity(r), "data")
        for key, value in r["pvalues"].items():
            assert abs(value - saved["pvalues"][key]) < 1e-12, (name, identity(r), key)
    checks.append(
        dict(
            file=str(path),
            datasets=len(new),
            data_hashes="identical",
            pvalues="matched",
        )
    )
assert checks, "No quick calculations found; see study README commands."
(QUICK / "verification.json").write_text(
    json.dumps(dict(status="passed", checks=checks), indent=2) + "\n"
)
print(
    "PASS:", len(checks), "available calculation files compared with the paper records."
)
for result in checks:
    print(result)
