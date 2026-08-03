"""Canonical heterogeneous-flux burglary model.

There is intentionally no formulation switch.  Every simulation uses

    A_t = div(eta grad(A - Ast)) - A + rho A exp(-pi) + Ast,
    rho_t = div(grad(rho) - 2 rho/A grad(A))
            - (rho A - q) exp(-pi).

The zero-flux condition associated with the attractiveness equation is
``eta * grad(A - Ast) . n = 0``.
"""

from __future__ import annotations

from typing import Any

from urbancrime.coefficients import BurglaryCoefficients


def attractiveness_flux(A: Any, ast: Any, eta: Any) -> Any:
    """Return the canonical heterogeneous attractiveness flux."""

    import ufl

    return eta * ufl.grad(A - ast)


def burglary_residual(
    *,
    A: Any,
    A_previous: Any,
    rho: Any,
    rho_previous: Any,
    test_A: Any,
    test_rho: Any,
    dt: Any,
    coefficients: BurglaryCoefficients,
    pi: Any,
    measure: Any = None,
) -> Any:
    """Return the backward-Euler weak residual for the two persistent fields."""

    import ufl

    dx = measure if measure is not None else ufl.dx
    attenuation = ufl.exp(-pi)
    eta, ast, q = coefficients.eta, coefficients.ast, coefficients.q

    attractiveness = (
        ufl.inner(A - A_previous, test_A) * dx
        + dt * ufl.inner(attractiveness_flux(A, ast, eta), ufl.grad(test_A)) * dx
        + dt * ufl.inner(A - rho * A * attenuation - ast, test_A) * dx
    )
    criminals = (
        ufl.inner(rho - rho_previous, test_rho) * dx
        + dt
        * ufl.inner(
            ufl.grad(rho) - 2.0 * (rho / A) * ufl.grad(A),
            ufl.grad(test_rho),
        )
        * dx
        + dt * ufl.inner((rho * A - q) * attenuation, test_rho) * dx
    )
    return attractiveness + criminals
