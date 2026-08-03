"""Delayed adaptive police PDE-ODE closure."""

from __future__ import annotations

from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray


def delayed_police_residual(
    *,
    pi: Any,
    pi_previous: Any,
    H: Any,
    H_previous: Any,
    test_pi: Any,
    test_H: Any,
    A: Any,
    rho: Any,
    tau: Any,
    dt: Any,
    measure: Any = None,
) -> Any:
    """Return the consistent backward-Euler residual for ``pi`` and ``H``."""

    import ufl

    dx = measure if measure is not None else ufl.dx
    police = (
        ufl.inner(pi - pi_previous, test_pi) * dx
        + dt
        * ufl.inner(
            ufl.grad(pi) - 2.0 * (pi / H) * ufl.grad(H),
            ufl.grad(test_pi),
        )
        * dx
    )
    information = (
        ufl.inner(H - H_previous, test_H) * dx
        + dt * ufl.inner((H - rho * A * ufl.exp(-pi)) / tau, test_H) * dx
    )
    return police + information


def mass_lumped_h_update(
    H_previous: ArrayLike,
    A: ArrayLike,
    rho: ArrayLike,
    pi: ArrayLike,
    tau: ArrayLike,
    dt: float,
) -> NDArray[np.float64]:
    """Pointwise backward-Euler update used by the partitioned method."""

    H0, a, r, p, relaxation = np.broadcast_arrays(
        np.asarray(H_previous, dtype=float),
        np.asarray(A, dtype=float),
        np.asarray(rho, dtype=float),
        np.asarray(pi, dtype=float),
        np.asarray(tau, dtype=float),
    )
    if dt <= 0 or np.any(relaxation <= 0):
        raise ValueError("dt and tau must be positive")
    source = r * a * np.exp(-p)
    return (relaxation * H0 + dt * source) / (relaxation + dt)
