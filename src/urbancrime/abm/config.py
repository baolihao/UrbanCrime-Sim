"""Typed configuration for the Python agent-based model."""

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
class ABMModelConfig:
    eta: float
    omega: float
    theta: float
    ast: CoefficientSpec
    gamma: CoefficientSpec
    generation_distribution: Literal["bernoulli", "poisson"]


@dataclass(frozen=True)
class ABMInitialConfig:
    dynamic_attractiveness: float
    mean_criminals_per_site: float
    seed: int


@dataclass(frozen=True)
class ABMPoliceConfig:
    strategy: Literal["none", "delayed"]
    tau: float | None
    initial_mean_density: float
    diffusivity: float
    diffusion_substeps: int
    information_floor: float


@dataclass(frozen=True)
class ABMOutputConfig:
    root: Path


@dataclass(frozen=True)
class ABMConfig:
    schema_version: int
    name: str
    grid: ABMGridConfig
    time: ABMTimeConfig
    model: ABMModelConfig
    initial: ABMInitialConfig
    policing: ABMPoliceConfig
    output: ABMOutputConfig


def _map(value: Any, name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise TypeError(f"{name} must be a mapping")
    return value


def load_abm_config(path: str | Path) -> ABMConfig:
    with Path(path).open(encoding="utf-8") as stream:
        raw = _map(yaml.safe_load(stream), str(path))
    grid_raw = _map(raw["grid"], "grid")
    time_raw = _map(raw["time"], "time")
    model_raw = _map(raw["model"], "model")
    initial_raw = _map(raw["initial"], "initial")
    police_raw = _map(raw["policing"], "policing")
    output_raw = _map(raw["output"], "output")
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
        model=ABMModelConfig(
            eta=float(model_raw["eta"]),
            omega=float(model_raw["omega"]),
            theta=float(model_raw["theta"]),
            ast=CoefficientSpec.from_mapping(_map(model_raw["ast"], "model.ast")),
            gamma=CoefficientSpec.from_mapping(_map(model_raw["gamma"], "model.gamma")),
            generation_distribution=str(model_raw["generation_distribution"]),  # type: ignore[arg-type]
        ),
        initial=ABMInitialConfig(
            dynamic_attractiveness=float(initial_raw["dynamic_attractiveness"]),
            mean_criminals_per_site=float(initial_raw["mean_criminals_per_site"]),
            seed=int(initial_raw["seed"]),
        ),
        policing=ABMPoliceConfig(
            strategy=str(police_raw["strategy"]),  # type: ignore[arg-type]
            tau=float(police_raw["tau"]) if police_raw.get("tau") is not None else None,
            initial_mean_density=float(police_raw["initial_mean_density"]),
            diffusivity=float(police_raw["diffusivity"]),
            diffusion_substeps=int(police_raw["diffusion_substeps"]),
            information_floor=float(police_raw["information_floor"]),
        ),
        output=ABMOutputConfig(root=Path(output_raw["root"])),
    )
    _validate(config)
    return config


def _validate(config: ABMConfig) -> None:
    if config.grid.nx < 2 or config.grid.ny < 2 or config.grid.spacing <= 0:
        raise ValueError("ABM grid requires nx, ny >= 2 and positive spacing")
    if config.grid.boundary not in {"no_flow", "periodic"}:
        raise ValueError("ABM boundary must be no_flow or periodic")
    if config.time.step <= 0 or config.time.final <= config.time.start or config.time.save_every < 1:
        raise ValueError("invalid ABM time controls")
    _ = config.time.steps
    if not 0 <= config.model.eta <= 1 or config.model.omega <= 0 or config.model.theta < 0:
        raise ValueError("invalid ABM eta, omega, or theta")
    if config.model.generation_distribution not in {"bernoulli", "poisson"}:
        raise ValueError("invalid criminal generation distribution")
    if config.initial.dynamic_attractiveness < 0 or config.initial.mean_criminals_per_site < 0:
        raise ValueError("ABM initial means must be nonnegative")
    if config.policing.strategy not in {"none", "delayed"}:
        raise ValueError("Python ABM supports none and delayed policing")
    if config.policing.strategy == "delayed" and (config.policing.tau is None or config.policing.tau <= 0):
        raise ValueError("delayed ABM requires positive tau")
    if (
        config.policing.initial_mean_density < 0
        or config.policing.diffusivity < 0
        or config.policing.diffusion_substeps < 1
        or config.policing.information_floor <= 0
    ):
        raise ValueError("invalid ABM police controls")
