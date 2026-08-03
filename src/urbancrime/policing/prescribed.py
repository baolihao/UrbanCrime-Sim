"""Prescribed police fields, including the budget-normalized snapshot rule."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def budget_normalized_snapshot(
    attractiveness: ArrayLike,
    *,
    weights: ArrayLike,
    budget: float,
    kappa: float,
    cutoff: float,
    positivity_floor: float,
) -> NDArray[np.float64]:
    """Construct the frozen police density used by the snapshot strategy.

    The unscaled profile is ``-log(0.5 * (1 - tanh(kappa*(A-cutoff))))``.
    It is then scaled so that the mass-lumped integral of ``pi`` equals the
    requested total budget.
    """

    A = np.asarray(attractiveness, dtype=float)
    w = np.asarray(weights, dtype=float)
    if A.shape != w.shape:
        raise ValueError("attractiveness and weights must have the same shape")
    if budget < 0 or kappa <= 0 or positivity_floor <= 0:
        raise ValueError("budget must be nonnegative; kappa and floor must be positive")
    if np.any(w < 0) or not np.any(w > 0):
        raise ValueError("integration weights must be nonnegative with positive total mass")
    if budget == 0:
        return np.zeros_like(A)

    attenuation = 0.5 * (1.0 - np.tanh(kappa * (A - cutoff)))
    raw = -np.log(np.clip(attenuation, positivity_floor, 1.0))
    raw_budget = float(np.dot(w, raw))
    if raw_budget <= positivity_floor:
        return np.full_like(A, budget / float(w.sum()))
    return raw * (budget / raw_budget)
