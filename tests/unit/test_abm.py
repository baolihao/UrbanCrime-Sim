from dataclasses import replace

import numpy as np

from urbancrime.abm import BurglaryABM, load_abm_config


def test_delayed_abm_police_transport_conserves_budget() -> None:
    config = load_abm_config("configs/abm/delayed_square.yaml")
    config = replace(
        config,
        grid=replace(config.grid, nx=8, ny=7),
        time=replace(config.time, final=0.01),
    )
    model = BurglaryABM(config)
    model.H *= 1.0 + 0.2 * np.sin(model.x)
    initial = float(model.pi.sum())
    model.step()
    np.testing.assert_allclose(model.pi.sum(), initial, rtol=1.0e-11, atol=1.0e-11)
    assert np.min(model.pi) >= 0.0
