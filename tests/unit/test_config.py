from pathlib import Path

from urbancrime.abm import load_abm_config
from urbancrime.config import NoPoliceConfig, load_run_config


ROOT = Path(__file__).resolve().parents[2]


def test_no_police_config_is_fully_typed() -> None:
    config = load_run_config(ROOT / "configs/simulations/no_police_square.yaml")
    assert isinstance(config.policing, NoPoliceConfig)
    assert config.model.q.value == 1.0
    assert config.time.steps == 10_000
    assert config.resolved_mapping()["derived"]["time_steps"] == 10_000


def test_smoke_config_has_no_formulation_switch() -> None:
    config = load_run_config(ROOT / "configs/smoke/no_police.yaml")
    assert "formulation" not in config.resolved_mapping()


def test_all_continuum_configs_expose_state_tolerance_and_no_historical_q_names() -> None:
    for path in sorted((ROOT / "configs").rglob("*.yaml")):
        if path.parts[-2] in {"abm", "studies"} or path.name in {
            "policy_comparison.yaml",
            "stability.yaml",
        }:
            continue
        config = load_run_config(path)
        assert config.solver.state_tolerance > 0
        resolved = config.resolved_mapping()
        text = str(resolved)
        assert "Bbar" not in text
        assert "GTfactor" not in text


def test_all_abm_configs_are_typed_and_have_no_hidden_legacy_names() -> None:
    for path in sorted((ROOT / "configs/abm").glob("*.yaml")):
        config = load_abm_config(path)
        assert config.schema_version == 2
        assert config.time.steps > 0
        text = path.read_text(encoding="utf-8")
        assert "Bbar" not in text
        assert "n_total" not in text


def test_m3as_no_noise_figure_configs_have_no_stochastic_initial_data() -> None:
    for case in (2, 3):
        config = load_run_config(
            ROOT / f"configs/simulations/no_police_case{case}_no_noise_square.yaml"
        )
        assert config.initial_conditions.noise is None
        assert config.initial_conditions.A["kind"] == "constant"
        assert config.initial_conditions.rho["kind"] == "constant"


def test_extension_abm_configs_save_the_published_field_panel_times() -> None:
    expected = {
        8: (765.0, 768.0, 771.0, 774.0, 777.0, 780.0),
        9: (744.0, 747.0, 750.0, 753.0, 756.0, 759.0),
    }
    for case, times in expected.items():
        config = load_abm_config(ROOT / f"configs/abm/delayed_case{case}.yaml")
        assert config.output.snapshot_times == times
