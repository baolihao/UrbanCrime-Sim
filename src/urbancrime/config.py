"""Typed, serializable configuration for simulations and parameter studies.

Every value that can change a numerical result belongs in a configuration file.
Derived quantities (for example the number of time steps or a total police budget
computed from a mean density) are calculated once and written to the resolved run
metadata.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Literal, Mapping, TypeAlias

import yaml


def _mapping(value: Any, name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise TypeError(f"{name} must be a mapping")
    return value


@dataclass(frozen=True)
class CoefficientSpec:
    """Description of a scalar FEM coefficient.

    YAML supports constants, registered analytic profiles, and fields stored on
    disk. Python callers may additionally pass callables directly to
    :func:`urbancrime.coefficients.as_coefficient`.
    """

    kind: Literal["constant", "analytic", "file"]
    value: float | None = None
    profile: str | None = None
    parameters: Mapping[str, Any] = field(default_factory=dict)
    path: Path | None = None
    field_name: str | None = None

    @classmethod
    def from_mapping(cls, raw: Mapping[str, Any]) -> "CoefficientSpec":
        kind = str(raw.get("kind", ""))
        if kind not in {"constant", "analytic", "file"}:
            raise ValueError(f"unknown coefficient kind: {kind!r}")
        spec = cls(
            kind=kind,  # type: ignore[arg-type]
            value=float(raw["value"]) if raw.get("value") is not None else None,
            profile=raw.get("profile"),
            parameters=dict(raw.get("parameters", {})),
            path=Path(raw["path"]) if raw.get("path") else None,
            field_name=raw.get("field"),
        )
        if spec.kind == "constant" and spec.value is None:
            raise ValueError("constant coefficient requires 'value'")
        if spec.kind == "analytic" and not spec.profile:
            raise ValueError("analytic coefficient requires 'profile'")
        if spec.kind == "file" and spec.path is None:
            raise ValueError("file coefficient requires 'path'")
        return spec


@dataclass(frozen=True)
class ModelConfig:
    eta: CoefficientSpec
    ast: CoefficientSpec
    q: CoefficientSpec

    @classmethod
    def from_mapping(cls, raw: Mapping[str, Any]) -> "ModelConfig":
        return cls(
            eta=CoefficientSpec.from_mapping(_mapping(raw.get("eta"), "model.eta")),
            ast=CoefficientSpec.from_mapping(_mapping(raw.get("ast"), "model.ast")),
            q=CoefficientSpec.from_mapping(_mapping(raw.get("q"), "model.q")),
        )


@dataclass(frozen=True)
class BudgetConfig:
    total: float | None = None
    mean_density: float | None = None

    @classmethod
    def from_mapping(cls, raw: Mapping[str, Any]) -> "BudgetConfig":
        total = float(raw["total"]) if raw.get("total") is not None else None
        mean = float(raw["mean_density"]) if raw.get("mean_density") is not None else None
        if (total is None) == (mean is None):
            raise ValueError("police budget requires exactly one of total or mean_density")
        if total is not None and total < 0:
            raise ValueError("total police budget must be nonnegative")
        if mean is not None and mean < 0:
            raise ValueError("mean police density must be nonnegative")
        return cls(total=total, mean_density=mean)

    def resolve_total(self, domain_area: float) -> float:
        if domain_area <= 0:
            raise ValueError("domain_area must be positive")
        return self.total if self.total is not None else float(self.mean_density) * domain_area


@dataclass(frozen=True)
class NoPoliceConfig:
    strategy: Literal["none"] = "none"


@dataclass(frozen=True)
class PrescribedPoliceConfig:
    strategy: Literal["prescribed"]
    profile: str
    activation_time: float
    budget: BudgetConfig
    positivity_floor: float
    parameters: Mapping[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class OptimalPoliceConfig:
    strategy: Literal["optimal"]
    activation_time: float
    budget: BudgetConfig
    update_each_iteration: bool
    positivity_floor: float
    root_tolerance: float
    root_max_iterations: int


@dataclass(frozen=True)
class DelayedPoliceConfig:
    strategy: Literal["delayed"]
    tau: CoefficientSpec
    initial_mean_density: float
    h_discretization: Literal["consistent", "mass_lumped"]


PoliceConfig: TypeAlias = (
    NoPoliceConfig | PrescribedPoliceConfig | OptimalPoliceConfig | DelayedPoliceConfig
)


def police_config_from_mapping(raw: Mapping[str, Any]) -> PoliceConfig:
    strategy = str(raw.get("strategy", ""))
    if strategy == "none":
        return NoPoliceConfig()
    if strategy == "prescribed":
        result = PrescribedPoliceConfig(
            strategy="prescribed",
            profile=str(raw["profile"]),
            activation_time=float(raw["activation_time"]),
            budget=BudgetConfig.from_mapping(_mapping(raw.get("budget"), "policing.budget")),
            parameters=dict(raw.get("parameters", {})),
            positivity_floor=float(raw["positivity_floor"]),
        )
        if result.activation_time < 0 or result.positivity_floor <= 0:
            raise ValueError("prescribed activation time and positivity floor are invalid")
        return result
    if strategy == "optimal":
        result = OptimalPoliceConfig(
            strategy="optimal",
            activation_time=float(raw["activation_time"]),
            budget=BudgetConfig.from_mapping(_mapping(raw.get("budget"), "policing.budget")),
            update_each_iteration=bool(raw["update_each_iteration"]),
            positivity_floor=float(raw["positivity_floor"]),
            root_tolerance=float(raw["root_tolerance"]),
            root_max_iterations=int(raw["root_max_iterations"]),
        )
        if result.activation_time < 0 or result.positivity_floor <= 0:
            raise ValueError("optimal activation time and positivity floor are invalid")
        if result.root_tolerance <= 0 or result.root_max_iterations < 1:
            raise ValueError("optimal root-solver controls are invalid")
        return result
    if strategy == "delayed":
        result = DelayedPoliceConfig(
            strategy="delayed",
            tau=CoefficientSpec.from_mapping(_mapping(raw.get("tau"), "policing.tau")),
            initial_mean_density=float(raw["initial_mean_density"]),
            h_discretization=str(raw["h_discretization"]),  # type: ignore[arg-type]
        )
        if result.initial_mean_density < 0:
            raise ValueError("initial police mean density must be nonnegative")
        if result.h_discretization not in {"consistent", "mass_lumped"}:
            raise ValueError("h_discretization must be consistent or mass_lumped")
        return result
    raise ValueError(f"unknown policing strategy: {strategy!r}")


@dataclass(frozen=True)
class DomainConfig:
    kind: Literal["rectangle", "polygon"]
    parameters: Mapping[str, Any]

    @classmethod
    def from_mapping(cls, raw: Mapping[str, Any]) -> "DomainConfig":
        kind = str(raw.get("kind", ""))
        if kind not in {"rectangle", "polygon"}:
            raise ValueError(f"unknown domain kind: {kind!r}")
        return cls(kind=kind, parameters={k: v for k, v in raw.items() if k != "kind"})  # type: ignore[arg-type]


@dataclass(frozen=True)
class MeshConfig:
    cell_type: Literal["triangle", "quadrilateral"]
    polynomial_degree: int
    nx: int | None = None
    ny: int | None = None
    characteristic_length: float | None = None
    boundary_length: float | None = None
    interior_length: float | None = None
    transition_width: float | None = None

    @classmethod
    def from_mapping(cls, raw: Mapping[str, Any]) -> "MeshConfig":
        config = cls(
            cell_type=str(raw["cell_type"]),  # type: ignore[arg-type]
            polynomial_degree=int(raw["polynomial_degree"]),
            nx=int(raw["nx"]) if raw.get("nx") is not None else None,
            ny=int(raw["ny"]) if raw.get("ny") is not None else None,
            characteristic_length=(
                float(raw["characteristic_length"])
                if raw.get("characteristic_length") is not None
                else None
            ),
            boundary_length=(
                float(raw["boundary_length"])
                if raw.get("boundary_length") is not None
                else None
            ),
            interior_length=(
                float(raw["interior_length"])
                if raw.get("interior_length") is not None
                else None
            ),
            transition_width=(
                float(raw["transition_width"])
                if raw.get("transition_width") is not None
                else None
            ),
        )
        if config.cell_type not in {"triangle", "quadrilateral"}:
            raise ValueError("mesh.cell_type must be triangle or quadrilateral")
        if config.polynomial_degree < 1:
            raise ValueError("mesh.polynomial_degree must be positive")
        if (config.nx is None) != (config.ny is None):
            raise ValueError("mesh.nx and mesh.ny must be provided together")
        if config.nx is not None and (config.nx < 1 or config.ny is None or config.ny < 1):
            raise ValueError("mesh.nx and mesh.ny must be positive")
        return config


@dataclass(frozen=True)
class TimeConfig:
    start: float
    final: float
    step: float
    save_every: int

    @classmethod
    def from_mapping(cls, raw: Mapping[str, Any]) -> "TimeConfig":
        config = cls(
            start=float(raw["start"]),
            final=float(raw["final"]),
            step=float(raw["step"]),
            save_every=int(raw["save_every"]),
        )
        if config.final <= config.start or config.step <= 0 or config.save_every < 1:
            raise ValueError("invalid time interval, time step, or save cadence")
        return config

    @property
    def steps(self) -> int:
        raw = (self.final - self.start) / self.step
        result = round(raw)
        if abs(raw - result) > 1.0e-10:
            raise ValueError("time interval must be an integer multiple of time.step")
        return int(result)


@dataclass(frozen=True)
class SolverConfig:
    method: Literal["implicit", "partitioned"]
    relative_tolerance: float
    absolute_tolerance: float
    max_iterations: int
    state_tolerance: float
    convergence_criterion: Literal["incremental", "residual"]
    block_order: tuple[str, ...]
    petsc_options: Mapping[str, str]

    @classmethod
    def from_mapping(cls, raw: Mapping[str, Any]) -> "SolverConfig":
        config = cls(
            method=str(raw["method"]),  # type: ignore[arg-type]
            relative_tolerance=float(raw["relative_tolerance"]),
            absolute_tolerance=float(raw["absolute_tolerance"]),
            max_iterations=int(raw["max_iterations"]),
            state_tolerance=float(raw["state_tolerance"]),
            convergence_criterion=str(raw["convergence_criterion"]),  # type: ignore[arg-type]
            block_order=tuple(str(value) for value in raw["block_order"]),
            petsc_options={str(k): str(v) for k, v in raw["petsc_options"].items()},
        )
        if config.method not in {"implicit", "partitioned"}:
            raise ValueError("solver.method must be implicit or partitioned")
        if config.relative_tolerance <= 0 or config.absolute_tolerance <= 0:
            raise ValueError("solver tolerances must be positive")
        if config.max_iterations < 1:
            raise ValueError("solver.max_iterations must be positive")
        if config.state_tolerance <= 0:
            raise ValueError("solver.state_tolerance must be positive")
        if config.convergence_criterion not in {"incremental", "residual"}:
            raise ValueError("invalid Newton convergence criterion")
        if config.method == "partitioned" and not config.block_order:
            raise ValueError("partitioned solver requires an explicit block_order")
        if len(config.block_order) != len(set(config.block_order)):
            raise ValueError("solver.block_order entries must be unique")
        return config


@dataclass(frozen=True)
class InitialConditionConfig:
    A: Mapping[str, Any]
    rho: Mapping[str, Any]
    pi: Mapping[str, Any] | None = None
    H: Mapping[str, Any] | None = None
    noise: Mapping[str, Any] | None = None

    @classmethod
    def from_mapping(cls, raw: Mapping[str, Any]) -> "InitialConditionConfig":
        return cls(
            A=dict(_mapping(raw.get("A"), "initial_conditions.A")),
            rho=dict(_mapping(raw.get("rho"), "initial_conditions.rho")),
            pi=dict(raw["pi"]) if raw.get("pi") else None,
            H=dict(raw["H"]) if raw.get("H") else None,
            noise=dict(raw["noise"]) if raw.get("noise") else None,
        )


@dataclass(frozen=True)
class OutputConfig:
    root: Path
    fields: tuple[str, ...]
    format: str

    @classmethod
    def from_mapping(cls, raw: Mapping[str, Any]) -> "OutputConfig":
        return cls(
            root=Path(raw["root"]),
            fields=tuple(str(item) for item in raw["fields"]),
            format=str(raw["format"]),
        )


@dataclass(frozen=True)
class RunConfig:
    schema_version: int
    name: str
    model: ModelConfig
    policing: PoliceConfig
    domain: DomainConfig
    mesh: MeshConfig
    time: TimeConfig
    initial_conditions: InitialConditionConfig
    solver: SolverConfig
    output: OutputConfig

    @classmethod
    def from_mapping(cls, raw: Mapping[str, Any]) -> "RunConfig":
        return cls(
            schema_version=int(raw["schema_version"]),
            name=str(raw["name"]),
            model=ModelConfig.from_mapping(_mapping(raw.get("model"), "model")),
            policing=police_config_from_mapping(_mapping(raw.get("policing"), "policing")),
            domain=DomainConfig.from_mapping(_mapping(raw.get("domain"), "domain")),
            mesh=MeshConfig.from_mapping(_mapping(raw.get("mesh"), "mesh")),
            time=TimeConfig.from_mapping(_mapping(raw.get("time"), "time")),
            initial_conditions=InitialConditionConfig.from_mapping(
                _mapping(raw.get("initial_conditions"), "initial_conditions")
            ),
            solver=SolverConfig.from_mapping(_mapping(raw.get("solver"), "solver")),
            output=OutputConfig.from_mapping(_mapping(raw.get("output"), "output")),
        )

    def resolved_mapping(self) -> dict[str, Any]:
        result = _serialize(asdict(self))
        result["derived"] = {"time_steps": self.time.steps}
        return result


def _serialize(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, Mapping):
        return {str(k): _serialize(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [_serialize(v) for v in value]
    return value


def load_run_config(path: str | Path) -> RunConfig:
    source = Path(path)
    with source.open("r", encoding="utf-8") as stream:
        raw = yaml.safe_load(stream)
    return RunConfig.from_mapping(_mapping(raw, str(source)))


def write_resolved_config(config: RunConfig, path: str | Path) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8") as stream:
        yaml.safe_dump(config.resolved_mapping(), stream, sort_keys=False)
