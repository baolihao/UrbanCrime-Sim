"""Typed configuration for the discrete burglary-and-police ABM.

Configuration values use the nondimensional variables printed in the papers.
The dimensional constants in :class:`ABMScalingConfig` are only used to map
integer agents to the continuum fields ``rho`` and ``pi``.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal, Mapping

import yaml

from urbancrime.config import CoefficientSpec


@dataclass(frozen=True)
class ABMGridConfig:
    nx: int
    ny: int
    spacing: float
    boundary: Literal["no_flow", "periodic"]


@dataclass(frozen=True)
class ABMTimeConfig:
    start: float
    final: float
    step: float
    save_every: int

    @property
    def steps(self) -> int:
        raw = (self.final - self.start) / self.step
        steps = round(raw)
        if abs(raw - steps) > 1.0e-10:
            raise ValueError("ABM interval must be divisible by its time step")
        return steps


@dataclass(frozen=True)
class ABMScalingConfig:
    """Dimensional constants used in the paper's nondimensionalization."""

    diffusivity: float
    omega: float
    theta: float
    beta: float


@dataclass(frozen=True)
class ABMModelConfig:
    eta: float
    ast: CoefficientSpec
    generation_ratio: CoefficientSpec
    arrest_probability: float
    generation_distribution: Literal["bernoulli", "poisson"]


@dataclass(frozen=True)
class ABMInitialConfig:
    dynamic_attractiveness: CoefficientSpec
    criminal_density: CoefficientSpec
    police_density: CoefficientSpec
    information: CoefficientSpec | None
    seed: int


@dataclass(frozen=True)
class ABMPoliceConfig:
    strategy: Literal["none", "delayed"]
    tau: float | None


@dataclass(frozen=True)
class ABMOutputConfig:
    root: Path
    snapshot_times: tuple[float, ...]


@dataclass(frozen=True)
class ABMConfig:
    schema_version: int
    name: str
    grid: ABMGridConfig
    time: ABMTimeConfig
    scaling: ABMScalingConfig
    model: ABMModelConfig
    initial: ABMInitialConfig
    policing: ABMPoliceConfig
    output: ABMOutputConfig


def _map(value: Any, name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise TypeError(f"{name} must be a mapping")
    return value


def _coefficient(raw: Any, name: str, base: Path) -> CoefficientSpec:
    spec = CoefficientSpec.from_mapping(_map(raw, name))
    if spec.path is not None and not spec.path.is_absolute():
        spec = CoefficientSpec(
            kind=spec.kind,
            value=spec.value,
            profile=spec.profile,
            parameters=spec.parameters,
            path=(base / spec.path).resolve(),
            field_name=spec.field_name,
        )
    return spec


def load_abm_config(path: str | Path) -> ABMConfig:
    config_path = Path(path)
    with config_path.open(encoding="utf-8") as stream:
        raw = _map(yaml.safe_load(stream), str(path))
    grid_raw = _map(raw["grid"], "grid")
    time_raw = _map(raw["time"], "time")
    scaling_raw = _map(raw["scaling"], "scaling")
    model_raw = _map(raw["model"], "model")
    initial_raw = _map(raw["initial"], "initial")
    police_raw = _map(raw["policing"], "policing")
    output_raw = _map(raw["output"], "output")
    information_raw = initial_raw.get("information")
    config = ABMConfig(
        schema_version=int(raw["schema_version"]),
        name=str(raw["name"]),
        grid=ABMGridConfig(
            nx=int(grid_raw["nx"]),
            ny=int(grid_raw["ny"]),
            spacing=float(grid_raw["spacing"]),
            boundary=str(grid_raw["boundary"]),  # type: ignore[arg-type]
        ),
        time=ABMTimeConfig(
            start=float(time_raw["start"]),
            final=float(time_raw["final"]),
            step=float(time_raw["step"]),
            save_every=int(time_raw["save_every"]),
        ),
        scaling=ABMScalingConfig(
            diffusivity=float(scaling_raw["diffusivity"]),
            omega=float(scaling_raw["omega"]),
            theta=float(scaling_raw["theta"]),
            beta=float(scaling_raw["beta"]),
        ),
        model=ABMModelConfig(
            eta=float(model_raw["eta"]),
            ast=_coefficient(model_raw["ast"], "model.ast", config_path.parent),
            generation_ratio=_coefficient(
                model_raw["generation_ratio"], "model.generation_ratio", config_path.parent
            ),
            arrest_probability=float(model_raw["arrest_probability"]),
            generation_distribution=str(  # type: ignore[arg-type]
                model_raw["generation_distribution"]
            ),
        ),
        initial=ABMInitialConfig(
            dynamic_attractiveness=_coefficient(
                initial_raw["dynamic_attractiveness"],
                "initial.dynamic_attractiveness",
                config_path.parent,
            ),
            criminal_density=_coefficient(
                initial_raw["criminal_density"], "initial.criminal_density", config_path.parent
            ),
            police_density=_coefficient(
                initial_raw["police_density"], "initial.police_density", config_path.parent
            ),
            information=(
                _coefficient(information_raw, "initial.information", config_path.parent)
                if information_raw is not None
                else None
            ),
            seed=int(initial_raw["seed"]),
        ),
        policing=ABMPoliceConfig(
            strategy=str(police_raw["strategy"]),  # type: ignore[arg-type]
            tau=float(police_raw["tau"]) if police_raw.get("tau") is not None else None,
        ),
        output=ABMOutputConfig(
            root=Path(output_raw["root"]),
            snapshot_times=tuple(float(value) for value in output_raw.get("snapshot_times", ())),
        ),
    )
    _validate(config)
    return config


def _validate(config: ABMConfig) -> None:
    if config.schema_version != 2:
        raise ValueError("Python ABM configurations require schema_version: 2")
    if config.grid.nx < 2 or config.grid.ny < 2 or config.grid.spacing <= 0:
        raise ValueError("ABM grid requires nx, ny >= 2 and positive spacing")
    if config.grid.boundary not in {"no_flow", "periodic"}:
        raise ValueError("ABM boundary must be no_flow or periodic")
    if (
        config.time.step <= 0
        or config.time.final <= config.time.start
        or config.time.save_every < 1
    ):
        raise ValueError("invalid ABM time controls")
    _ = config.time.steps
    if config.time.step > 1:
        raise ValueError("the attractiveness decay update requires nondimensional time.step <= 1")
    scaling = config.scaling
    if min(scaling.diffusivity, scaling.omega, scaling.theta, scaling.beta) <= 0:
        raise ValueError("ABM scaling constants must be positive")
    if not 0 <= config.model.eta <= 1:
        raise ValueError("model.eta must lie in [0, 1]")
    if not 0 <= config.model.arrest_probability <= 1:
        raise ValueError("model.arrest_probability must lie in [0, 1]")
    if config.model.generation_distribution not in {"bernoulli", "poisson"}:
        raise ValueError("invalid criminal generation distribution")
    if config.policing.strategy not in {"none", "delayed"}:
        raise ValueError("Python ABM supports none and delayed policing")
    if config.policing.strategy == "none" and config.model.arrest_probability != 0:
        raise ValueError("arrest_probability must be zero without police")
    if config.policing.strategy == "delayed":
        if config.policing.tau is None or config.policing.tau <= 0:
            raise ValueError("delayed ABM requires positive tau")
        if config.time.step > config.policing.tau:
            raise ValueError("explicit information update requires time.step <= policing.tau")
    for value in config.output.snapshot_times:
        raw_step = (value - config.time.start) / config.time.step
        if (
            value < config.time.start
            or value > config.time.final
            or abs(raw_step - round(raw_step)) > 1e-8
        ):
            raise ValueError("snapshot times must lie on ABM time steps inside the run interval")
