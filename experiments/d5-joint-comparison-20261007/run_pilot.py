#!/usr/bin/env python3
"""Targeted statistic-geometry experiment for Appendix D.5.

This script reuses the retained v17 data/copy generators without modifying them.
It writes no output files: all summaries are emitted as JSON on stdout.

Goal
----
Compare the existing ridge-quadratic joint statistic with two max-type joint
statistics while keeping the data-generating model and randomization mechanism
fixed.

1. MaxCorr: maximum absolute correlation between the common background-ridge
   residual of Y and each background-ridge residualized target column.
2. MaxPartial: for each target j, use the same partial-residual construction as
   the graph-based coordinatewise baseline, convert the squared partial
   correlation to its local Beta-tail p-value, and take -log(min_j p_j).
   The rank CRT calibrates this maximum jointly; there is no Bonferroni factor
   inside the joint statistic.

The observed target matrix and every randomized target matrix receive exactly
the same scoring rule.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import beta

HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE
V17_CODE = HERE / "src"
sys.path.insert(0, str(V17_CODE))

from mscrt.core import (  # noqa: E402
    RidgeStatistic,
    coordinate_pvalues,
    g_copies,
    oracle_copies,
    rank_pvalue,
)
from mscrt.experiments import generator, graph_design  # noqa: E402


DESIGNS = {
    "A": {"seed": 27092611, "bandwidth": 6, "offdiag": 0.20},
    "B": {"seed": 27092612, "bandwidth": 1, "offdiag": -0.45},
}


def _batch_residual(values: np.ndarray, design: np.ndarray) -> np.ndarray:
    """Residualize one vector per batch against a full-rank batch design."""
    q, r = np.linalg.qr(design, mode="reduced")
    diag = np.abs(np.diagonal(r, axis1=-2, axis2=-1))
    if np.any(diag < 1e-10):
        raise RuntimeError("Unexpected rank deficiency in batch projection.")
    coef = np.einsum("bnk,bn->bk", q, values)
    return values - np.einsum("bnk,bk->bn", q, coef)


def _broadcast_fixed_columns(x: np.ndarray, cols: np.ndarray, batch: int) -> np.ndarray:
    if len(cols) == 0:
        return np.empty((batch, x.shape[0], 0), dtype=float)
    return np.broadcast_to(x[:, cols][None, :, :], (batch, x.shape[0], len(cols)))


def maxcorr_scores(xt: np.ndarray, statistic: RidgeStatistic) -> np.ndarray:
    """Max absolute background-ridge residual correlation."""
    xt = np.asarray(xt, float)
    if xt.ndim == 2:
        xt = xt[None, :, :]
    z = np.einsum("ij,bjk->bik", statistic.h, xt)
    y = statistic.yres
    num = np.abs(np.einsum("bnt,n->bt", z, y))
    den = np.sqrt(np.sum(z * z, axis=1)) * max(np.linalg.norm(y), 1e-30)
    corr = np.divide(num, den, out=np.zeros_like(num), where=den > 1e-30)
    return np.max(corr, axis=1)


def maxpartial_scores(
    xt: np.ndarray,
    x_observed: np.ndarray,
    target: np.ndarray,
    graph: np.ndarray,
    statistic: RidgeStatistic,
) -> np.ndarray:
    """Max-type score matching the coordinatewise partial-residual geometry.

    For each candidate and target j:
      * residualize the common response residual on the other candidate target
        residuals;
      * partial both X_j and that response vector on j's graph neighbors;
      * convert the squared partial correlation to the same local Beta-tail
        quantity used by the coordinatewise baseline.

    The joint score is -log(min local p). The rank CRT, rather than Bonferroni,
    calibrates the maximum across coordinates.
    """
    xt = np.asarray(xt, float)
    if xt.ndim == 2:
        xt = xt[None, :, :]
    batch, n, t = xt.shape
    if not np.array_equal(target, np.arange(t)):
        raise ValueError("This experiment assumes target columns are reordered first.")

    z = np.einsum("ij,bjk->bik", statistic.h, xt)
    ones = np.ones((batch, n, 1), dtype=float)
    yres = np.broadcast_to(statistic.yres[None, :], (batch, n))
    local_p = np.ones((batch, t), dtype=float)

    for j in range(t):
        others = np.arange(t) != j
        response_design = np.concatenate((ones, z[:, :, others]), axis=2)
        r_y = _batch_residual(yres, response_design)

        nb = np.flatnonzero(graph[j])
        nb_t = nb[nb < t]
        nb_s = nb[nb >= t]
        parts = [ones]
        if len(nb_t):
            parts.append(xt[:, :, nb_t])
        if len(nb_s):
            parts.append(_broadcast_fixed_columns(x_observed, nb_s, batch))
        d = np.concatenate(parts, axis=2)

        r_j = _batch_residual(xt[:, :, j], d)
        u_j = _batch_residual(r_y, d)
        dot = np.einsum("bn,bn->b", r_j, u_j)
        den = np.einsum("bn,bn->b", r_j, r_j) * np.einsum("bn,bn->b", u_j, u_j)
        w = np.divide(dot * dot, den, out=np.zeros(batch), where=den > 1e-26)
        w = np.clip(w, 0.0, 1.0)
        df = n - d.shape[2]
        if df <= 1:
            local_p[:, j] = 1.0
        else:
            local_p[:, j] = beta.sf(w, 0.5, (df - 1) / 2)

    min_p = np.min(local_p, axis=1)
    return -np.log(np.maximum(min_p, 1e-300))


def summarize_binary(values: np.ndarray) -> dict:
    values = np.asarray(values, float)
    n = len(values)
    rate = float(values.mean())
    mcse = float(values.std(ddof=1) / np.sqrt(n)) if n > 1 else 0.0
    return {"reps": n, "rate": rate, "mcse": mcse, "rejections": int(values.sum())}


def paired(a: np.ndarray, b: np.ndarray) -> dict:
    d = np.asarray(a, float) - np.asarray(b, float)
    n = len(d)
    return {
        "difference": float(d.mean()),
        "paired_mcse": float(d.std(ddof=1) / np.sqrt(n)) if n > 1 else 0.0,
    }


def run_design(
    design_name: str,
    reps: int,
    m: int,
    sweeps: int,
    patterns: list[str],
    thetas: list[float],
) -> dict:
    spec = DESIGNS[design_name]
    n, p, t = 80, 120, 8
    target = np.arange(t)
    omega, sigma, graph = graph_design(p, t, "connected", spec["bandwidth"], spec["offdiag"])
    chol = np.linalg.cholesky(sigma)

    # conditional mean parameters, rederived directly to keep the new script local
    ott = omega[np.ix_(target, target)]
    s = np.setdiff1d(np.arange(p), target)
    v = np.linalg.solve(ott, np.eye(t))
    b = -np.linalg.solve(ott, omega[np.ix_(target, s)])
    gamma = np.zeros(len(s))
    gamma[:4] = np.array([1.0, -0.7, 0.5, 0.3])

    cells: list[tuple[str, float]] = [("null", 0.0)]
    for pattern in patterns:
        for theta in thetas:
            if theta > 0:
                cells.append((pattern, float(theta)))

    result = {
        "design": design_name,
        "seed": spec["seed"],
        "bandwidth": spec["bandwidth"],
        "offdiag": spec["offdiag"],
        "M": m,
        "L": sweeps,
        "reps": reps,
        "target_arrangement": "connected",
        "cells": {},
    }

    for pattern, theta in cells:
        direction = np.ones(t) / np.sqrt(t) if pattern != "sparse" else np.eye(t)[0]
        direction = direction / np.sqrt(direction @ v @ direction)

        records: dict[str, list[bool]] = {
            "CoordG-Bonf": [],
            "OracleRidge": [],
            "GRidge-L": [],
            "OracleMaxCorr": [],
            "GMaxCorr-L": [],
            "OracleMaxPartial": [],
            "GMaxPartial-L": [],
        }

        for rep in range(reps):
            key = ("connected", pattern, theta, rep)
            rng = generator(spec["seed"], "data", *key)
            x = rng.normal(size=(n, p)) @ chol.T
            eps = rng.normal(size=n)
            mu = x[:, s] @ b.T
            y = x[:, s] @ gamma + theta * ((x[:, target] - mu) @ direction) + eps

            ridge = RidgeStatistic.fit(y, x[:, s], n * 0.2, 1.0)
            ridge0 = ridge(x[:, target])
            maxcorr0 = float(maxcorr_scores(x[:, target], ridge)[0])
            maxpartial0 = float(maxpartial_scores(x[:, target], x, target, graph, ridge)[0])

            coord_g, _ = coordinate_pvalues(x, target, graph, ridge, omega)
            records["CoordG-Bonf"].append(coord_g <= 0.05)

            oracle = oracle_copies(x, target, omega, m, generator(spec["seed"], "oracle", *key))
            records["OracleRidge"].append(rank_pvalue(ridge0, ridge(oracle)) <= 0.05)
            records["OracleMaxCorr"].append(
                rank_pvalue(maxcorr0, maxcorr_scores(oracle, ridge)) <= 0.05
            )
            records["OracleMaxPartial"].append(
                rank_pvalue(maxpartial0, maxpartial_scores(oracle, x, target, graph, ridge)) <= 0.05
            )

            gc = g_copies(
                x, target, graph, m, sweeps,
                generator(spec["seed"], "g", sweeps, *key),
            )
            gxt = gc[:, :, target]
            records["GRidge-L"].append(rank_pvalue(ridge0, ridge(gxt)) <= 0.05)
            records["GMaxCorr-L"].append(
                rank_pvalue(maxcorr0, maxcorr_scores(gxt, ridge)) <= 0.05
            )
            records["GMaxPartial-L"].append(
                rank_pvalue(maxpartial0, maxpartial_scores(gxt, x, target, graph, ridge)) <= 0.05
            )

        arr = {k: np.asarray(vv, bool) for k, vv in records.items()}
        cell = {
            "methods": {k: summarize_binary(vv) for k, vv in arr.items()},
            "paired_vs_coord": {
                k: paired(vv, arr["CoordG-Bonf"])
                for k, vv in arr.items()
                if k != "CoordG-Bonf"
            },
            "paired_maxpartial_vs_ridge": {
                "Oracle": paired(arr["OracleMaxPartial"], arr["OracleRidge"]),
                "G": paired(arr["GMaxPartial-L"], arr["GRidge-L"]),
            },
        }
        result["cells"][f"{pattern}:{theta:g}"] = cell

    return result


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", choices=["A", "B", "both"], default="both")
    ap.add_argument("--reps", type=int, default=100)
    ap.add_argument("--M", type=int, default=199)
    ap.add_argument("--L", type=int, default=3)
    ap.add_argument("--patterns", default="dense,sparse")
    ap.add_argument("--thetas", default="0.75")
    args = ap.parse_args()

    patterns = [x.strip() for x in args.patterns.split(",") if x.strip()]
    thetas = [float(x) for x in args.thetas.split(",") if x.strip()]
    designs = ["A", "B"] if args.design == "both" else [args.design]

    out = {
        "experiment": "joint-max-statistic-geometry-pilot",
        "note": "Exploratory; no manuscript claim until reviewed.",
        "designs": [
            run_design(d, args.reps, args.M, args.L, patterns, thetas)
            for d in designs
        ],
    }
    print(json.dumps(out, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
