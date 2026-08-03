import json
from pathlib import Path

import numpy as np
import pytest

from urbancrime.config import load_run_config
from urbancrime.runtime import run_pde


ROOT = Path(__file__).resolve().parents[2]


@pytest.mark.fenics
def test_constant_coefficient_refactor_matches_tagged_legacy_solver() -> None:
    config = load_run_config(ROOT / "configs/regression/legacy_no_police_constant.yaml")
    reference = json.loads(
        (ROOT / "tests/regression/data/legacy_no_police_constant.json").read_text()
    )["final"]
    result = run_pde(config, write_output=False)
    A = result.final_state.fields["A"]
    rho = result.final_state.fields["rho"]
    actual = {
        "A_min": float(A.min()),
        "A_max": float(A.max()),
        "A_mean_dof": float(A.mean()),
        "A_l2_dof": float(np.linalg.norm(A)),
        "rho_min": float(rho.min()),
        "rho_max": float(rho.max()),
        "rho_mean_dof": float(rho.mean()),
        "rho_l2_dof": float(np.linalg.norm(rho)),
    }
    for name, expected in reference.items():
        assert actual[name] == pytest.approx(expected, rel=2.0e-9, abs=2.0e-11)
