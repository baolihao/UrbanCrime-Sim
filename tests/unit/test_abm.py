from dataclasses import replace
import json

import numpy as np

from urbancrime.abm import BurglaryABM, load_abm_config, run_abm


def _small_delayed_config():
    config = load_abm_config("configs/abm/delayed_square.yaml")
    return replace(
        config,
        grid=replace(config.grid, nx=8, ny=7),
        time=replace(config.time, final=0.2, save_every=2),
        output=replace(config.output, snapshot_times=(0.0, 0.1, 0.2)),
    )


def test_initial_integer_populations_preserve_configured_totals() -> None:
    model = BurglaryABM(_small_delayed_config())
    expected_criminals = int(
        np.ceil(8 * 7 * 0.6 / model.criminal_density_per_agent)
    )
    expected_police = int(np.ceil(8 * 7 * 0.5 / model.police_density_per_agent))
    assert model.criminals.dtype == np.int64
    assert model.police_agents.dtype == np.int64
    assert int(model.criminals.sum()) == expected_criminals
    assert int(model.police_agents.sum()) == expected_police


def test_delayed_abm_police_agents_conserve_integer_budget() -> None:
    model = BurglaryABM(_small_delayed_config())
    model.H *= 1.0 + 0.2 * np.sin(model.x)
    initial = int(model.police_agents.sum())
    for _ in range(5):
        model.step()
        assert int(model.police_agents.sum()) == initial
        assert np.min(model.police_agents) >= 0
        np.testing.assert_array_equal(model.police_agents, model.police_agents.astype(int))


def test_information_update_uses_realized_events() -> None:
    config = _small_delayed_config()
    model = BurglaryABM(config)
    old_information = model.H.copy()
    result = model.step()
    assert config.policing.tau is not None
    expected = old_information + config.time.step / config.policing.tau * (
        result.realized_crime_rate - old_information
    )
    np.testing.assert_allclose(model.H, expected)


def test_total_events_do_not_depend_on_output_frequency(tmp_path) -> None:
    base = _small_delayed_config()
    first = replace(base, time=replace(base.time, save_every=1))
    second = replace(base, time=replace(base.time, save_every=5))
    run_abm(first, output_directory=tmp_path / "every-step")
    run_abm(second, output_directory=tmp_path / "every-five")
    summary_first = json.loads((tmp_path / "every-step" / "summary.json").read_text())
    summary_second = json.loads((tmp_path / "every-five" / "summary.json").read_text())
    assert summary_first["total_events"] == summary_second["total_events"]
    assert summary_first["total_generated"] == summary_second["total_generated"]
    assert summary_first["final_criminals"] == summary_second["final_criminals"]


def test_no_police_case_keeps_pi_zero() -> None:
    config = load_abm_config("configs/abm/no_police_case1.yaml")
    config = replace(
        config,
        grid=replace(config.grid, nx=6, ny=5),
        time=replace(config.time, final=0.04),
    )
    model = BurglaryABM(config)
    for _ in range(2):
        model.step()
    assert not np.any(model.police_agents)
    assert not np.any(model.pi)
    assert not np.any(model.H)
