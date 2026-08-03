import numpy as np

from urbancrime.policing.delayed import mass_lumped_h_update
from urbancrime.policing.optimal import optimal_police_density
from urbancrime.policing.prescribed import budget_normalized_snapshot


def test_snapshot_budget_uses_weights() -> None:
    A = np.array([0.5, 1.0, 2.0])
    weights = np.array([0.1, 0.2, 0.7])
    pi = budget_normalized_snapshot(
        A,
        weights=weights,
        budget=0.6,
        kappa=5.0,
        cutoff=1.5,
        positivity_floor=1.0e-12,
    )
    np.testing.assert_allclose(np.dot(weights, pi), 0.6, rtol=1.0e-12, atol=1.0e-12)


def test_optimal_budget_and_complementarity() -> None:
    intensity = np.array([0.5, 3.0, 1.0, 8.0])
    weights = np.array([0.1, 0.2, 0.3, 0.4])
    pi = optimal_police_density(
        intensity,
        weights=weights,
        budget=0.75,
        positivity_floor=1.0e-15,
        root_tolerance=1.0e-13,
        max_iterations=200,
    )
    np.testing.assert_allclose(np.dot(weights, pi), 0.75, rtol=1.0e-11)
    assert np.all(pi >= 0)
    assert pi[np.argmax(intensity)] == np.max(pi)


def test_mass_lumped_h_update_preserves_positive_state() -> None:
    result = mass_lumped_h_update(
        H_previous=np.array([1.0, 2.0]),
        A=np.array([2.0, 3.0]),
        rho=np.array([0.5, 0.25]),
        pi=np.array([0.1, 0.2]),
        tau=np.array([1.0, 2.0]),
        dt=0.1,
    )
    assert np.all(result > 0)
