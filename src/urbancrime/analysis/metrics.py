"""Mesh-independent scalar diagnostics."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike


def crime_intensity(A: ArrayLike, rho: ArrayLike, pi: ArrayLike) -> np.ndarray:
    return np.asarray(A, dtype=float) * np.asarray(rho, dtype=float) * np.exp(
        -np.asarray(pi, dtype=float)
    )


def weighted_average(values: ArrayLike, weights: ArrayLike) -> float:
    field = np.asarray(values, dtype=float)
    measure = np.asarray(weights, dtype=float)
    if field.shape != measure.shape:
        raise ValueError("values and weights must have the same shape")
    total = float(measure.sum())
    if total <= 0:
        raise ValueError("weights must have positive sum")
    return float(np.dot(field, measure) / total)
