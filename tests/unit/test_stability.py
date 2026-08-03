import numpy as np

from urbancrime.analysis.stability import (
    StabilityParameters,
    critical_delay,
    dominant_growth_rate,
    equilibrium,
    equilibrium_residual,
    hurwitz_determinant,
    hurwitz_determinant_grid,
    jacobian,
    no_police_subsystem_jacobian,
    transversality,
)


PARAMETERS = StabilityParameters(
    eta=0.3,
    ast_bar=1.0 / 50.0,
    q_bar=1.5,
    pi_bar=0.5,
    domain_length=10.0,
)


def test_equilibrium_uses_effective_source_without_extra_attenuation() -> None:
    steady = equilibrium(PARAMETERS)
    np.testing.assert_allclose(steady.rho, PARAMETERS.q_bar / steady.A)
    np.testing.assert_allclose(equilibrium_residual(PARAMETERS), 0.0, atol=1.0e-13)


def test_hopf_point_crosses_direct_eigenvalue_boundary() -> None:
    point = critical_delay(
        PARAMETERS,
        max_mode=8,
        tau_min=0.1,
        tau_max=100.0,
        bracket_samples=500,
        root_tolerance=1.0e-10,
    )
    assert point is not None
    below = dominant_growth_rate(0.99 * point.tau, PARAMETERS, max_mode=8)[0]
    above = dominant_growth_rate(1.01 * point.tau, PARAMETERS, max_mode=8)[0]
    assert below < 0 < above
    assert abs(transversality(point, PARAMETERS, relative_step=1.0e-5)) > 1.0e-8


def test_vectorized_hurwitz_matches_direct_characteristic_polynomial() -> None:
    delays = np.geomspace(0.03, 300.0, 25)
    vectorized = hurwitz_determinant_grid(0.7, delays, PARAMETERS)
    direct = np.asarray([hurwitz_determinant(0.7, tau, PARAMETERS) for tau in delays])
    np.testing.assert_allclose(vectorized, direct, rtol=1.0e-10, atol=1.0e-10)


def test_zero_police_limit_factorizes_into_classical_subsystem() -> None:
    parameters = StabilityParameters(
        eta=0.3,
        ast_bar=1.0 / 50.0,
        q_bar=1.5,
        pi_bar=0.0,
        domain_length=10.0,
    )
    mu, tau = 0.7, 5.0
    full = np.sort_complex(np.linalg.eigvals(jacobian(mu, tau, parameters)))
    expected = np.sort_complex(
        np.concatenate(
            [
                np.linalg.eigvals(no_police_subsystem_jacobian(mu, parameters)),
                np.asarray([-mu, -1.0 / tau]),
            ]
        )
    )
    np.testing.assert_allclose(full, expected, rtol=1.0e-12, atol=1.0e-12)
    point = critical_delay(
        parameters,
        max_mode=8,
        tau_min=0.1,
        tau_max=100.0,
        bracket_samples=300,
        root_tolerance=1.0e-10,
    )
    assert point is None
