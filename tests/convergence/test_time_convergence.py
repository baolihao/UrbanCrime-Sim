from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from urbancrime.config import TimeConfig, load_run_config
from urbancrime.runtime import run_pde


ROOT = Path(__file__).resolve().parents[2]


@pytest.mark.fenics
@pytest.mark.slow
def test_backward_euler_has_first_order_time_convergence() -> None:
    base = load_run_config(ROOT / "configs/smoke/no_police.yaml")
    initial = np.array([base.initial_conditions.A["value"], base.initial_conditions.rho["value"]])
    ast = float(base.model.ast.value)
    q = float(base.model.q.value)

    def right_hand_side(_time: float, state: np.ndarray) -> np.ndarray:
        A, rho = state
        return np.array([-A + rho * A + ast, -(rho * A - q)])

    reference = solve_ivp(
        right_hand_side, (0.0, 0.1), initial, rtol=1.0e-12, atol=1.0e-14
    ).y[:, -1]
    errors = []
    for dt in (0.05, 0.025):
        config = replace(
            base,
            time=TimeConfig(start=0.0, final=0.1, step=dt, save_every=1),
        )
        result = run_pde(config, write_output=False)
        final = np.array(
            [result.final_state.fields["A"][0], result.final_state.fields["rho"][0]]
        )
        errors.append(float(np.linalg.norm(final - reference)))
    assert errors[0] / errors[1] > 1.7
