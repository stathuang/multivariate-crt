"""Candidate-symmetric MaxPartial with observable local diagnostics.

Full-rank arithmetic follows run_pilot.py. Rank-deficient designs use the retained
SVD rank convention, applied separately and identically to every candidate.
"""
from pathlib import Path
import sys
import numpy as np
from scipy.stats import beta

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import run_pilot as reference
from mscrt.core import basis, rank_pvalue


def residual(values, design):
    q, r = np.linalg.qr(design, mode='reduced')
    diag = np.abs(np.diagonal(r, axis1=-2, axis2=-1))
    ranks = np.full(len(values), design.shape[-1], dtype=int)
    good = np.all(diag >= 1e-10, axis=1) & (design.shape[-1] <= design.shape[-2])
    result = values - np.einsum('bnk,bk->bn', q, np.einsum('bnk,bn->bk', q, values))
    for i in np.flatnonzero(~good):
        qi = basis(design[i])
        ranks[i] = qi.shape[1]
        result[i] = values[i] - qi @ (qi.T @ values[i])
    return result, ranks


def maxpartial(xt, x_fixed, target, graph, statistic, details=False):
    xt = np.asarray(xt, float)
    if xt.ndim == 2:
        xt = xt[None]
    batch, n, t = xt.shape
    if not np.array_equal(target, np.arange(t)):
        raise ValueError('Targets must be reordered first.')
    if not np.all(np.isfinite(xt)):
        raise ValueError('Candidate matrices must be finite.')
    z = np.einsum('ij,bjk->bik', statistic.h, xt)
    ones = np.ones((batch, n, 1))
    yres = np.broadcast_to(statistic.yres[None], (batch, n))
    local_p = np.ones((batch, t))
    w_all = np.zeros((batch, t))
    df_all = np.zeros((batch, t), dtype=int)
    for j in range(t):
        a = np.concatenate((ones, z[:, :, np.arange(t) != j]), axis=2)
        ry, _ = residual(yres, a)
        nb = np.flatnonzero(graph[j])
        nt, ns = nb[nb < t], nb[nb >= t]
        parts = [ones]
        if len(nt):
            parts.append(xt[:, :, nt])
        if len(ns):
            parts.append(reference._broadcast_fixed_columns(x_fixed, ns, batch))
        d = np.concatenate(parts, axis=2)
        r, ranks = residual(xt[:, :, j], d)
        u, _ = residual(ry, d)
        dot = np.einsum('bn,bn->b', r, u)
        den = np.einsum('bn,bn->b', r, r)*np.einsum('bn,bn->b', u, u)
        w = np.clip(np.divide(dot*dot, den, out=np.zeros(batch), where=den > 1e-26), 0, 1)
        df = n-ranks
        ok = (df > 1) & (den > 1e-26)
        local_p[ok, j] = beta.sf(w[ok], .5, (df[ok]-1)/2)
        w_all[:, j], df_all[:, j] = w, df
    score = -np.log(np.maximum(local_p.min(axis=1), 1e-300))
    if details:
        return score, dict(local_p=local_p, squared_partial_correlation=w_all, df=df_all)
    return score


def bonferroni_two(p_ridge, p_partial):
    return min(1., 2*min(p_ridge, p_partial))
