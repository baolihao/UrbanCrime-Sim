import numpy as np

from urbancrime.abm import burglary_probability_from_fields


def test_abm_event_rate_converges_to_pde_reaction_rate() -> None:
    A = 1.7
    rho = 2.5
    pi = 0.4
    dt = 1.0e-4
    agents = 10_000_000
    probability = float(burglary_probability_from_fields(A, pi, dt))
    rng = np.random.default_rng(2026)
    observed_rate = rng.binomial(agents, probability) / (agents * dt) * rho
    pde_rate = rho * A * np.exp(-pi)
    np.testing.assert_allclose(observed_rate, pde_rate, rtol=0.06)
