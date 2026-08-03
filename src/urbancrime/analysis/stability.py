"""Homogeneous-state stability and Hopf bifurcation calculations.

This analysis is deliberately restricted to constant coefficients on a square
domain.  It is not a stability solver for spatially heterogeneous fields.  The
zero Neumann mode is excluded because perturbations must preserve total police
mass.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Iterable

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import brentq


@dataclass(frozen=True)
class StabilityParameters:
    eta: float
    ast_bar: float
    q_bar: float
    pi_bar: float
    domain_length: float

    def validate(self) -> None:
        if self.eta <= 0 or self.ast_bar <= 0 or self.q_bar <= 0:
            raise ValueError("eta, ast_bar, and q_bar must be positive")
        if self.pi_bar < 0 or self.domain_length <= 0:
            raise ValueError("pi_bar must be nonnegative and domain_length must be positive")


@dataclass(frozen=True)
class Equilibrium:
    A: float
    rho: float
    pi: float
    H: float


@dataclass(frozen=True)
class HopfPoint:
    mu: float
    tau: float
    omega: float

    @property
    def period(self) -> float:
        return 2.0 * np.pi / self.omega


def equilibrium(parameters: StabilityParameters) -> Equilibrium:
    parameters.validate()
    attenuation = np.exp(-parameters.pi_bar)
    A_bar = parameters.ast_bar + parameters.q_bar * attenuation
    rho_bar = parameters.q_bar / A_bar
    H_bar = rho_bar * A_bar * attenuation
    return Equilibrium(A=A_bar, rho=rho_bar, pi=parameters.pi_bar, H=H_bar)


def equilibrium_residual(parameters: StabilityParameters) -> NDArray[np.float64]:
    steady = equilibrium(parameters)
    attenuation = np.exp(-steady.pi)
    return np.array(
        [
            -steady.A + steady.rho * steady.A * attenuation + parameters.ast_bar,
            -(steady.rho * steady.A - parameters.q_bar) * attenuation,
            steady.rho * steady.A * attenuation - steady.H,
        ],
        dtype=float,
    )


def neumann_eigenvalues_square(domain_length: float, max_mode: int) -> NDArray[np.float64]:
    if domain_length <= 0 or max_mode < 1:
        raise ValueError("domain_length and max_mode must be positive")
    values = {
        round(((m * np.pi / domain_length) ** 2 + (n * np.pi / domain_length) ** 2), 14)
        for m in range(max_mode + 1)
        for n in range(max_mode + 1)
        if m or n
    }
    return np.asarray(sorted(values), dtype=float)


def jacobian(mu: float, tau: float, parameters: StabilityParameters) -> NDArray[np.float64]:
    if mu <= 0 or tau <= 0:
        raise ValueError("mu and tau must be positive")
    steady = equilibrium(parameters)
    attenuation = np.exp(-steady.pi)
    rho_e = steady.rho * attenuation
    A_e = steady.A * attenuation
    rho_A_e = steady.rho * steady.A * attenuation
    return np.array(
        [
            [
                -parameters.eta * mu - 1.0 + rho_e,
                A_e,
                -rho_A_e,
                0.0,
            ],
            [
                2.0 * mu * steady.rho / steady.A - rho_e,
                -mu - A_e,
                0.0,
                0.0,
            ],
            [0.0, 0.0, -mu, 2.0 * mu * steady.pi / steady.H],
            [rho_e / tau, A_e / tau, -rho_A_e / tau, -1.0 / tau],
        ],
        dtype=float,
    )


def no_police_subsystem_jacobian(
    mu: float, parameters: StabilityParameters
) -> NDArray[np.float64]:
    """Return the classical two-field block at ``pi_bar = 0``."""

    if parameters.pi_bar != 0:
        raise ValueError("no-police subsystem requires pi_bar = 0")
    if mu <= 0:
        raise ValueError("mu must be positive")
    steady = equilibrium(parameters)
    return np.array(
        [
            [
                -parameters.eta * mu - 1.0 + steady.rho,
                steady.A,
            ],
            [
                2.0 * mu * steady.rho / steady.A - steady.rho,
                -mu - steady.A,
            ],
        ],
        dtype=float,
    )


def characteristic_coefficients(
    mu: float, tau: float, parameters: StabilityParameters
) -> tuple[float, float, float, float]:
    coefficients = np.real_if_close(np.poly(jacobian(mu, tau, parameters)))
    if np.iscomplexobj(coefficients):
        raise FloatingPointError("characteristic coefficients unexpectedly became complex")
    return tuple(float(value) for value in coefficients[1:])  # type: ignore[return-value]


def hurwitz_determinant(mu: float, tau: float, parameters: StabilityParameters) -> float:
    a3, a2, a1, a0 = characteristic_coefficients(mu, tau, parameters)
    return a1 * a2 * a3 - a1**2 - a0 * a3**2


def _characteristic_affine_in_inverse_tau(
    mu: float, parameters: StabilityParameters
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Return ``alpha, beta`` such that ``a(tau)=alpha+beta/tau``."""

    at_one = np.asarray(characteristic_coefficients(mu, 1.0, parameters))
    at_two = np.asarray(characteristic_coefficients(mu, 2.0, parameters))
    beta = 2.0 * (at_one - at_two)
    alpha = 2.0 * at_two - at_one
    return alpha, beta


def hurwitz_determinant_grid(
    mu: float, tau: NDArray[np.float64], parameters: StabilityParameters
) -> NDArray[np.float64]:
    """Evaluate the determinant on an unchanged tau grid without repeated eigensolves.

    Only one Jacobian row depends on the delay and it is proportional to
    ``1/tau``. Every characteristic coefficient is therefore affine in inverse
    delay; two direct evaluations determine it exactly up to roundoff.
    """

    delays = np.asarray(tau, dtype=float)
    if np.any(delays <= 0):
        raise ValueError("tau values must be positive")
    alpha, beta = _characteristic_affine_in_inverse_tau(mu, parameters)
    coefficients = alpha[:, None] + beta[:, None] / delays[None, :]
    a3, a2, a1, a0 = coefficients
    return a1 * a2 * a3 - a1**2 - a0 * a3**2


def is_hopf_candidate(mu: float, tau: float, parameters: StabilityParameters) -> bool:
    a3, a2, a1, a0 = characteristic_coefficients(mu, tau, parameters)
    return min(a0, a1, a2, a3) > 0.0 and a1 / a3 > 0.0


def critical_delay_for_mode(
    mu: float,
    parameters: StabilityParameters,
    *,
    tau_min: float,
    tau_max: float,
    bracket_samples: int,
    root_tolerance: float,
) -> HopfPoint | None:
    if tau_min <= 0 or tau_max <= tau_min or bracket_samples < 2 or root_tolerance <= 0:
        raise ValueError("invalid critical-delay search parameters")
    grid = np.geomspace(tau_min, tau_max, bracket_samples)
    values = hurwitz_determinant_grid(mu, grid, parameters)
    changes = np.flatnonzero(np.signbit(values[:-1]) != np.signbit(values[1:]))
    candidates: list[HopfPoint] = []
    for index in changes:
        tau = brentq(
            lambda value: hurwitz_determinant(mu, value, parameters),
            float(grid[index]),
            float(grid[index + 1]),
            xtol=root_tolerance,
        )
        if is_hopf_candidate(mu, tau, parameters):
            a3, _, a1, _ = characteristic_coefficients(mu, tau, parameters)
            candidates.append(HopfPoint(mu=mu, tau=tau, omega=float(np.sqrt(a1 / a3))))
    return min(candidates, key=lambda point: point.tau) if candidates else None


def critical_delay(
    parameters: StabilityParameters,
    *,
    max_mode: int,
    tau_min: float,
    tau_max: float,
    bracket_samples: int,
    root_tolerance: float,
) -> HopfPoint | None:
    modes = neumann_eigenvalues_square(parameters.domain_length, max_mode)
    points = [
        point
        for mu in modes
        if (
            point := critical_delay_for_mode(
                float(mu),
                parameters,
                tau_min=tau_min,
                tau_max=tau_max,
                bracket_samples=bracket_samples,
                root_tolerance=root_tolerance,
            )
        )
        is not None
    ]
    return min(points, key=lambda point: point.tau) if points else None


def dominant_growth_rate(
    tau: float, parameters: StabilityParameters, *, max_mode: int
) -> tuple[float, float, float]:
    best = (-np.inf, 0.0, np.nan)
    for mu in neumann_eigenvalues_square(parameters.domain_length, max_mode):
        eigenvalues = np.linalg.eigvals(jacobian(float(mu), tau, parameters))
        value = eigenvalues[int(np.argmax(eigenvalues.real))]
        if value.real > best[0]:
            best = (float(value.real), float(abs(value.imag)), float(mu))
    return best


def transversality(
    point: HopfPoint,
    parameters: StabilityParameters,
    *,
    relative_step: float,
) -> float:
    if relative_step <= 0 or relative_step >= 1:
        raise ValueError("relative_step must lie in (0, 1)")
    delta = point.tau * relative_step
    target = 1j * point.omega
    plus = min(
        np.linalg.eigvals(jacobian(point.mu, point.tau + delta, parameters)),
        key=lambda value: abs(value - target),
    )
    minus = min(
        np.linalg.eigvals(jacobian(point.mu, point.tau - delta, parameters)),
        key=lambda value: abs(value - target),
    )
    return float((plus.real - minus.real) / (2.0 * delta))


def parameter_sweep(
    base: StabilityParameters,
    parameter_name: str,
    values: Iterable[float],
    **critical_delay_options: float | int,
) -> NDArray[np.float64]:
    if parameter_name not in {"eta", "ast_bar", "q_bar", "pi_bar", "domain_length"}:
        raise ValueError(f"unsupported stability parameter: {parameter_name}")
    output: list[float] = []
    for value in values:
        parameters = replace(base, **{parameter_name: float(value)})
        point = critical_delay(parameters, **critical_delay_options)  # type: ignore[arg-type]
        output.append(np.inf if point is None else point.tau)
    return np.asarray(output, dtype=float)
