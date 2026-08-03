"""Mass-constrained instantaneous optimal police allocation."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def optimal_police_density(
    crime_intensity: ArrayLike,
    *,
    weights: ArrayLike,
    budget: float,
    positivity_floor: float,
    root_tolerance: float,
    max_iterations: int,
) -> NDArray[np.float64]:
    """Solve the nodal, mass-lumped Zipkin allocation problem.

    The optimizer is ``pi_i = max(log(f_i/lambda), 0)``.  The multiplier
    ``lambda`` is found from ``sum_i weights_i*pi_i = budget``.
    """

    intensity = np.maximum(np.asarray(crime_intensity, dtype=float), positivity_floor)
    w = np.asarray(weights, dtype=float)
    if intensity.shape != w.shape:
        raise ValueError("crime_intensity and weights must have the same shape")
    if budget < 0 or positivity_floor <= 0 or root_tolerance <= 0:
        raise ValueError("invalid budget, positivity floor, or root tolerance")
    if max_iterations < 1 or np.any(w < 0) or not np.any(w > 0):
        raise ValueError("invalid iteration count or integration weights")
    if budget == 0:
        return np.zeros_like(intensity)

    def residual(multiplier: float) -> float:
        pi = np.maximum(np.log(intensity / multiplier), 0.0)
        return float(np.dot(w, pi) - budget)

    upper = float(np.max(intensity))
    lower = upper
    for _ in range(max_iterations):
        lower *= 0.5
        if residual(lower) >= 0:
            break
    else:
        raise RuntimeError("failed to bracket the optimal-police multiplier")

    for _ in range(max_iterations):
        midpoint = 0.5 * (lower + upper)
        value = residual(midpoint)
        if abs(value) <= root_tolerance * max(1.0, budget):
            lower = upper = midpoint
            break
        if value > 0:
            lower = midpoint
        else:
            upper = midpoint
    multiplier = 0.5 * (lower + upper)
    result = np.maximum(np.log(intensity / multiplier), 0.0)
    error = abs(float(np.dot(w, result)) - budget)
    if error > 10.0 * root_tolerance * max(1.0, budget):
        raise RuntimeError(f"optimal police budget residual is too large: {error}")
    return result
