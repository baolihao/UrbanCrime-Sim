from pathlib import Path

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
