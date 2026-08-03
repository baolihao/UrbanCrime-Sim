from pathlib import Path

import numpy as np
import pytest

from urbancrime.config import load_run_config
from urbancrime.runtime import run_pde
from urbancrime.runtime import PDERunner


ROOT = Path(__file__).resolve().parents[2]


@pytest.mark.fenics
def test_implicit_two_field_smoke() -> None:
    config = load_run_config(ROOT / "configs/smoke/no_police.yaml")
    result = run_pde(config, write_output=False)
    assert result.times[-1] == pytest.approx(0.1)
    assert np.all(np.isfinite(result.metrics["crime_average"]))


@pytest.mark.fenics
def test_partitioned_delayed_smoke_conserves_police_mass() -> None:
    config = load_run_config(ROOT / "configs/smoke/delayed.yaml")
    result = run_pde(config, write_output=False)
    np.testing.assert_allclose(result.metrics["police_mass"], 2.0, rtol=1.0e-10, atol=1.0e-10)


@pytest.mark.fenics
def test_state_validation_rejects_negative_rho_and_pi() -> None:
    no_police = PDERunner(
        load_run_config(ROOT / "configs/smoke/no_police.yaml"), write_output=False
    )
    no_police.fields["rho"].x.array[0] = -1.0e-3
    with pytest.raises(FloatingPointError, match="rho became negative"):
        no_police._validate_state()

    delayed = PDERunner(
        load_run_config(ROOT / "configs/smoke/delayed.yaml"), write_output=False
    )
    delayed.fields["pi"].x.array[0] = -1.0e-3
    with pytest.raises(FloatingPointError, match="pi became negative"):
        delayed._validate_state()
