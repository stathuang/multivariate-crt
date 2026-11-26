"""Gaussian graphical CRT sampling and the statistics used in Table S2.

No source code from the authors' repository is included. All statistic functions
must apply the same fixed rule to the observation and every randomized copy.
"""

from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from numpy.typing import NDArray
from scipy.special import betainc
from scipy.stats import beta, norm

Array = NDArray[np.float64]


def validate_graph(graph: Array, p: int) -> NDArray[np.bool_]:
    a = np.asarray(graph, dtype=bool)
    if a.shape != (p, p) or not np.array_equal(a, a.T) or np.any(np.diag(a)):
        raise ValueError(
            "Graph must be a symmetric p by p adjacency matrix with zero diagonal."
        )
    return a


def basis(a: Array) -> Array:
    """Rank-revealing orthonormal basis; used outside the fast Gaussian batch path."""
    a = np.asarray(a, float)
    if a.shape[1] == 0:
        return np.zeros((a.shape[0], 0))
    u, s, _ = np.linalg.svd(a, full_matrices=False)
    tol = np.finfo(float).eps * max(a.shape) * (s[0] if len(s) else 0.0)
    return u[:, s > tol]


def _batch_rotate(
    x: Array,
    j: int,
    nb: NDArray,
    rng: np.random.Generator,
    cached_q: Array | None = None,
) -> None:
    """In-place column update for independent branches, with fresh normals."""
    batch, n, _ = x.shape
    # This is the conservative identity convention, not a pseudoinverse
    # extension to rank-deficient n <= degree+1 designs.
    if n <= len(nb) + 1:
        return
    if cached_q is None:
        d = np.concatenate((np.ones((batch, n, 1)), x[:, :, nb]), axis=2)
        q, rr = np.linalg.qr(d, mode="reduced")
        diag = np.abs(np.diagonal(rr, axis1=-2, axis2=-1))
        # Gaussian designs have full rank a.s.; handle unexpected degeneracy safely.
        if np.any(diag < 1e-12):
            for k in range(batch):
                qk = basis(d[k])
                col = x[k, :, j]
                f = qk @ (qk.T @ col)
                r = col - f
                w = rng.normal(size=n)
                w -= qk @ (qk.T @ w)
                a, bnorm = np.linalg.norm(r), np.linalg.norm(w)
                if a > 1e-14 and bnorm > 1e-14:
                    x[k, :, j] = f + a * w / bnorm
            return
        qt = q.swapaxes(-1, -2)
        col = x[:, :, j, None]
        f = (q @ (qt @ col))[..., 0]
        w = rng.normal(size=(batch, n, 1))
        w = (w - q @ (qt @ w))[..., 0]
    else:
        q = cached_q
        col = x[:, :, j]
        f = (col @ q) @ q.T
        w = rng.normal(size=(batch, n))
        w -= (w @ q) @ q.T
    r = x[:, :, j] - f
    a = np.linalg.norm(r, axis=1)
    bnorm = np.linalg.norm(w, axis=1)
    ok = (a > 1e-14) & (bnorm > 1e-14)
    x[ok, :, j] = f[ok] + w[ok] * (a[ok] / bnorm[ok])[:, None]


def g_copies(
    x: Array,
    target: NDArray,
    graph: Array,
    m: int,
    sweeps: int,
    rng: np.random.Generator,
    order: list[int] | None = None,
) -> Array:
    """Hub-and-spoke construction: forward hub and independent branches in full reverse order.

    No convergence assumption or empirical stopping rule is used. Graph, target,
    order and sweeps must be fixed before viewing the inference data.
    """
    x = np.asarray(x, float)
    if x.ndim != 2 or not np.all(np.isfinite(x)) or m < 1 or sweeps < 0:
        raise ValueError("Invalid data, copy count or sweep count.")
    t = np.asarray(target, dtype=int)
    if len(set(t)) != len(t) or len(t) == 0 or np.any((t < 0) | (t >= x.shape[1])):
        raise ValueError("Invalid target indices.")
    a = validate_graph(graph, x.shape[1])
    order = list(map(int, t)) if order is None else list(order)
    if any(j not in t for j in order):
        raise ValueError("Updates may only change target columns.")
    sequence = order * sweeps
    neighbors = {j: np.flatnonzero(a[j]) for j in set(order)}
    cache = {}
    for j, nb in neighbors.items():
        if not set(nb).intersection(t) and x.shape[0] > len(nb) + 1:
            cache[j] = basis(np.column_stack((np.ones(x.shape[0]), x[:, nb])))
    hub = x[None, :, :].copy()
    for j in sequence:
        _batch_rotate(hub, j, neighbors[j], rng, cache.get(j))
    branches = np.repeat(hub, m, axis=0)
    for j in reversed(sequence):
        _batch_rotate(branches, j, neighbors[j], rng, cache.get(j))
    return branches


def conditional_parameters(
    omega: Array, target: NDArray
) -> tuple[NDArray, Array, Array]:
    p = omega.shape[0]
    t = np.asarray(target, int)
    s = np.setdiff1d(np.arange(p), t)
    ott = omega[np.ix_(t, t)]
    v = np.linalg.solve(ott, np.eye(len(t)))
    b = -np.linalg.solve(ott, omega[np.ix_(t, s)])
    return s, b, v


def oracle_copies(
    x: Array, target: NDArray, omega: Array, m: int, rng: np.random.Generator
) -> Array:
    s, b, v = conditional_parameters(omega, target)
    mu = x[:, s] @ b.T
    return (
        mu[None]
        + rng.normal(size=(m, x.shape[0], len(target))) @ np.linalg.cholesky(v).T
    )


def rank_pvalue(observed: float, randomized: Array, decimals: int = 12) -> float:
    """Conservative rank rule, with the SAME rounding applied to all statistics.

    Rounding is part of the statistic and prevents spurious rejection on an
    exactly degenerate fiber due to floating-point QR error.
    """
    z = np.asarray(randomized, float).ravel()
    if z.size < 1 or not np.isfinite(observed) or not np.all(np.isfinite(z)):
        raise ValueError("Statistics must be finite and at least one copy is needed.")
    return float(
        (1 + np.count_nonzero(np.round(z, decimals) >= np.round(observed, decimals)))
        / (len(z) + 1)
    )


@dataclass
class RidgeStatistic:
    """Fixed linear distillation based only on (Y, XS), plus a ridge group gain.

    Target residualization is recomputed by applying the SAME fixed smoother to
    every copy; no observed-target-fitted predictor is silently held fixed.
    """

    h: Array
    yres: Array
    penalty: float = 1.0

    @classmethod
    def fit(
        cls, y: Array, xs: Array, nuisance_penalty: float, penalty: float = 1.0
    ) -> "RidgeStatistic":
        if nuisance_penalty <= 0 or penalty <= 0:
            raise ValueError("Penalties must be positive.")
        n = len(y)
        c = np.eye(n) - np.ones((n, n)) / n
        z = c @ xs
        # Dual ridge inverse remains well-defined when s >= n.
        h = nuisance_penalty * np.linalg.solve(
            z @ z.T + nuisance_penalty * np.eye(n), c
        )
        h = (h + h.T) / 2
        return cls(h, h @ y, penalty)

    def __call__(self, xt: Array) -> Array | float:
        xt = np.asarray(xt, float)
        single = xt.ndim == 2
        if single:
            xt = xt[None]
        z = self.h @ xt
        gram = z.swapaxes(-1, -2) @ z
        q = z.swapaxes(-1, -2) @ self.yres
        coef = np.linalg.solve(gram + self.penalty * np.eye(z.shape[-1]), q[..., None])[
            ..., 0
        ]
        out = np.sum(q * coef, axis=-1) / max(float(self.yres @ self.yres), 1e-30)
        return float(out[0]) if single else out


def spherical_pvalue(xj: Array, d: Array, response: Array) -> float:
    """Exact two-sided local residual-score p-value; no Monte Carlo resolution loss."""
    q = basis(d)
    df = len(xj) - q.shape[1]
    r = xj - q @ (q.T @ xj)
    u = response - q @ (q.T @ response)
    den = float((r @ r) * (u @ u))
    if df <= 1 or den <= 1e-26:
        return 1.0
    w = np.clip(float(r @ u) ** 2 / den, 0, 1)
    return float(beta.sf(w, 0.5, (df - 1) / 2))


def coordinate_pvalues(
    x: Array,
    target: NDArray,
    graph: Array,
    statistic: RidgeStatistic,
    omega: Array | None = None,
) -> tuple[float, float | None]:
    """Exact graph-based coordinates and an oracle with the SAME residual score.

    The coordinate response vector projects the common distilled response off the
    OTHER target columns. It never uses the column being randomized.
    """
    n, p = x.shape
    t = np.asarray(target, int)
    pg, po = [], []
    for j in t:
        others = t[t != j]
        zothers = statistic.h @ x[:, others]
        u0 = basis(np.column_stack((np.ones(n), zothers)))
        r_y = statistic.yres - u0 @ (u0.T @ statistic.yres)
        nb = np.flatnonzero(graph[j])
        d = np.column_stack((np.ones(n), x[:, nb]))
        pg.append(spherical_pvalue(x[:, j], d, r_y))
        if omega is not None:
            q = basis(d)
            u = r_y - q @ (q.T @ r_y)
            # Correct graph implies u.T E[Xj | X_-j] = 0.
            denom = np.linalg.norm(u) / np.sqrt(omega[j, j])
            po.append(
                float(2 * norm.sf(abs(u @ x[:, j]) / denom)) if denom > 1e-13 else 1.0
            )
    return min(1.0, len(t) * min(pg)), (min(1.0, len(t) * min(po)) if po else None)
