#!/usr/bin/env python3
"""Table S2: coordinate, ridge, MaxCoefT, and Bonferroni comparisons."""
import os

for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[key] = "1"
import argparse, json, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "python"))
import numpy as np
from inputs import CONFIG, DESIGNS, generator, graph_design, array_sha
from mscrt.core import (
    RidgeStatistic,
    conditional_parameters,
    coordinate_pvalues,
    g_copies,
    oracle_copies,
    rank_pvalue,
)
import max_coefficient


def ridge_record(design, cell, rep):
    cfg = CONFIG
    spec = DESIGNS[design]
    n, p, t = cfg["n"], cfg["p"], cfg["t"]
    target = np.arange(t)
    pattern, theta = cell["pattern"], cell["theta"]
    omega, sigma, graph = graph_design(
        p, t, "connected", spec["bandwidth"], spec["offdiag"]
    )
    s, b, v = conditional_parameters(omega, target)
    direction = np.ones(t) / np.sqrt(t) if pattern != "sparse" else np.eye(t)[0]
    direction = direction / np.sqrt(direction @ v @ direction)
    gamma = np.zeros(len(s))
    gamma[:4] = [1.0, -0.7, 0.5, 0.3]
    key = ("connected", pattern, float(theta), int(rep))
    rng = generator(spec["seed"], "data", *key)
    x = rng.normal(size=(n, p)) @ np.linalg.cholesky(sigma).T
    y = (
        x[:, s] @ gamma
        + theta * ((x[:, target] - x[:, s] @ b.T) @ direction)
        + rng.normal(size=n)
    )
    stat = RidgeStatistic.fit(
        y, x[:, s], n * cfg["nuisance_penalty"], cfg["target_penalty"]
    )
    coord, _ = coordinate_pvalues(x, target, graph, stat, omega)
    pv = {"CoordG": coord}
    observed = stat(x[:, target])
    for sampler in ("Oracle", "G"):
        if sampler == "Oracle":
            copies = oracle_copies(
                x, target, omega, cfg["M"], generator(spec["seed"], "oracle", *key)
            )
        else:
            copies = g_copies(
                x,
                target,
                graph,
                cfg["M"],
                cfg["L"],
                generator(spec["seed"], "g", cfg["L"], *key),
            )[:, :, target]
        pv[sampler + "Ridge"] = rank_pvalue(observed, stat(copies))
    return dict(
        design=design,
        arrangement="connected",
        pattern=pattern,
        theta=theta,
        rep=rep,
        design_seed=spec["seed"],
        data_stream=["data", *key],
        data_sha256=array_sha(x, y),
        M=cfg["M"],
        L=cfg["L"],
        pvalues=pv,
        reject={k: int(v <= cfg["alpha"]) for k, v in pv.items()},
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--step", choices=["all", "ridge", "max"], default="all")
    a = parser.parse_args()
    out = ROOT / "results" / ("quick" if a.quick else "generated") / "03_joint_tests"
    out.mkdir(parents=True, exist_ok=True)
    base = out / "ridge.jsonl"
    target = out / "maxcoeft.jsonl"
    (out / "settings.json").write_text(json.dumps(CONFIG, indent=2) + "\n")
    if a.step in ("all", "ridge"):
        stop = CONFIG["rep_start"] + 2 if a.quick else CONFIG["rep_stop"]
        with base.open("x") as f:
            for d in CONFIG["designs"]:
                for cell in CONFIG["cells"]:
                    for rep in range(CONFIG["rep_start"], stop):
                        f.write(
                            json.dumps(ridge_record(d, cell, rep), sort_keys=True)
                            + "\n"
                        )
                        f.flush()
                    print("Ridge:", d, cell, flush=True)
    if a.step in ("all", "max"):
        rows = [json.loads(s) for s in base.read_text().splitlines()]
        expected = 28 if a.quick else 7000
        if len(rows) != expected:
            raise ValueError(f"Expected {expected} complete ridge records")
        for folder in ("scores", "records", "failures", "work"):
            (out / folder).mkdir(exist_ok=True)
        max_coefficient.init_worker(str(out))
        try:
            with target.open("x") as f:
                for i, row in enumerate(rows, 1):
                    result = max_coefficient.one_replication(row)
                    if not result["ok"]:
                        raise RuntimeError(result)
                    f.write(
                        (out / "records" / (result["name"] + ".json"))
                        .read_text()
                        .replace("\n", "")
                        + "\n"
                    )
                    f.flush()
                    if i % 25 == 0:
                        print("MaxCoefT:", i, "/", len(rows), flush=True)
        finally:
            max_coefficient._R.terminate()


if __name__ == "__main__":
    main()
