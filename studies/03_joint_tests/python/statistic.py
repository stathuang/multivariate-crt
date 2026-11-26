"""Maximum absolute target coefficient t, from one joint OLS fit per candidate."""

import numpy as np


def coefficient_t(candidates, response):
    x = np.asarray(candidates, dtype=float)
    y = np.asarray(response, dtype=float)
    single = x.ndim == 2
    if single:
        x = x[None]
    if x.ndim != 3 or y.shape != (x.shape[1],) or x.shape[2] < 1:
        raise ValueError("Expected candidate matrices and a common response vector")
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        raise ValueError("Nonfinite input")
    count, n, t = x.shape
    df = n - t - 1
    if df <= 0:
        raise ValueError("Positive residual degrees of freedom required")
    # An intercept and all target columns enter the SAME regression.
    design = np.concatenate((np.ones((count, n, 1)), x), axis=2)
    q, r = np.linalg.qr(design, mode="reduced")
    diagonal = np.abs(np.diagonal(r, axis1=-2, axis2=-1))
    # R's glm uses a relative QR tolerance of 1e-11 with its default controls.
    if np.any(diagonal <= 1e-11 * np.linalg.norm(design, axis=1)):
        raise ValueError("Target regression is rank deficient")
    qty = (q.swapaxes(-1, -2) @ y)[..., None]
    coefficients = np.linalg.solve(r, qty)[..., 0]
    residuals = y - (design @ coefficients[..., None])[..., 0]
    variance = np.sum(residuals**2, axis=1) / df
    if np.any(variance <= 0):
        raise ValueError("Zero residual variance")
    inverse_r = np.linalg.solve(r, np.broadcast_to(np.eye(t + 1), r.shape))
    standard_errors = np.sqrt(variance[:, None] * np.sum(inverse_r**2, axis=2))
    values = coefficients[:, 1:] / standard_errors[:, 1:]
    if not np.all(np.isfinite(values)):
        raise ValueError("Nonfinite target coefficient t")
    return values[0] if single else values


def max_coefficient_t(candidates, response):
    values = np.max(np.abs(coefficient_t(candidates, response)), axis=-1)
    return float(values) if np.ndim(values) == 0 else values
