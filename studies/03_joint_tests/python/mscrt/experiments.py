"""Pre-specified paired experiments; raw replicate records are always retained."""

from __future__ import annotations
import argparse, hashlib, json, platform, sys, time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import beta, f, norm


def generator(seed: int, *keys: object) -> np.random.Generator:
    """Stable, method-specific streams, invariant to method order and job sharding."""
    digest = hashlib.sha256(json.dumps([seed, *keys], sort_keys=True).encode()).digest()
    return np.random.default_rng(np.frombuffer(digest[:16], dtype="<u4"))


def graph_design(p: int, t: int, kind: str, bandwidth: int, offdiag: float):
    omega = np.eye(p)
    for j in range(1, bandwidth + 1):
        omega += offdiag * (np.eye(p, k=j) + np.eye(p, k=-j))
    positions = (
        np.arange(t) if kind == "connected" else np.linspace(0, p - 1, t, dtype=int)
    )
    if kind == "dispersed" and np.any(np.diff(positions) <= bandwidth):
        raise ValueError("Dispersed targets must have no within-target edges.")
    order = np.r_[positions, np.setdiff1d(np.arange(p), positions)]
    omega = omega[np.ix_(order, order)]
    if np.linalg.eigvalsh(omega).min() <= 0:
        raise ValueError("Precision is not positive definite.")
    sigma = np.linalg.inv(omega)
    scales = np.sqrt(np.diag(sigma))
    # Standardize each marginal to variance one.
    omega = scales[:, None] * omega * scales[None, :]
    sigma = np.linalg.inv(omega)
    graph = abs(omega) > 1e-12
    np.fill_diagonal(graph, False)
    return omega, sigma, graph
