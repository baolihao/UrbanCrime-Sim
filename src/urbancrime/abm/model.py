"""Discrete stochastic burglary model with mobile police agents.

The update follows equations (2.1)--(2.13) of the police-extension paper.  Both
criminals and police are integer populations.  Continuous arrays exposed by
the public API are their nondimensional continuum scalings, which makes ABM
output directly comparable with the PDE implementation.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray
import yaml

from urbancrime.coefficients import analytic_profile
from urbancrime.config import CoefficientSpec
from urbancrime.provenance import environment_metadata, git_revision, write_checksums

from .config import ABMConfig


FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]


@dataclass(frozen=True)
class ABMStep:
    """Stochastic counts and signal generated during one time interval."""

    events: IntArray
    generated: IntArray
    police_staying: IntArray
    realized_crime_rate: FloatArray


@dataclass(frozen=True)
class ABMResult:
    output_directory: Path
    times: FloatArray
    attractiveness: FloatArray
    criminal_density: FloatArray
    events: IntArray
    police: FloatArray
    information: FloatArray
    expected_crime_rate: FloatArray
    realized_crime_rate: FloatArray
    criminal_counts: IntArray
    police_counts: IntArray
    metric_times: FloatArray
    mean_attractiveness: FloatArray
    mean_criminal_density: FloatArray
    mean_police: FloatArray
    mean_information: FloatArray
    mean_expected_crime_rate: FloatArray
    mean_realized_crime_rate: FloatArray
    cumulative_events: IntArray


def _grid_field(spec: CoefficientSpec, x: FloatArray, y: FloatArray) -> FloatArray:
    shape = x.shape
    if spec.kind == "constant":
        return np.full(shape, float(spec.value), dtype=float)
    if spec.kind == "analytic":
        values = analytic_profile(spec)(np.vstack([x.ravel(), y.ravel()]))
        return np.asarray(values, dtype=float).reshape(shape)
    if spec.path is None or spec.path.suffix != ".npy":
        raise ValueError("ABM file fields use .npy arrays")
    values = np.asarray(np.load(spec.path), dtype=float)
    if values.shape != shape:
        raise ValueError(f"ABM field {spec.path} has shape {values.shape}, expected {shape}")
    return values


def _integer_population(expected: FloatArray, rng: np.random.Generator) -> IntArray:
    """Randomly round a density while preserving one deterministic total.

    The total is the ceiling used by the MATLAB reference.  Fractional agents
    are placed without replacement, so a uniform sub-unit density becomes a
    uniformly sampled set of occupied sites rather than independent Poisson
    noise in the initial condition.
    """

    if np.any(expected < 0) or not np.all(np.isfinite(expected)):
        raise ValueError("initial population densities must be finite and nonnegative")
    result = np.floor(expected).astype(np.int64)
    remainder = int(np.ceil(float(expected.sum()))) - int(result.sum())
    if remainder == 0:
        return result
    fractions = (expected - result).ravel()
    positive = np.flatnonzero(fractions > 0)
    if remainder > positive.size:
        raise RuntimeError("population rounding invariant failed")
    probabilities = fractions[positive] / fractions[positive].sum()
    chosen = rng.choice(positive, size=remainder, replace=False, p=probabilities)
    result.ravel()[chosen] += 1
    return result


class BurglaryABM:
    """Stateful realization of the discrete ABM in nondimensional variables."""

    def __init__(self, config: ABMConfig) -> None:
        self.config = config
        grid = config.grid
        coordinates_x = (np.arange(grid.nx) + 0.5) * grid.spacing
        coordinates_y = (np.arange(grid.ny) + 0.5) * grid.spacing
        self.x, self.y = np.meshgrid(coordinates_x, coordinates_y)
        self.rng = np.random.default_rng(config.initial.seed)

        self.ast = _grid_field(config.model.ast, self.x, self.y)
        self.generation_ratio = _grid_field(config.model.generation_ratio, self.x, self.y)
        self.B = _grid_field(config.initial.dynamic_attractiveness, self.x, self.y)
        initial_rho = _grid_field(config.initial.criminal_density, self.x, self.y)
        initial_pi = _grid_field(config.initial.police_density, self.x, self.y)
        if not np.all(np.isfinite(self.ast)) or not np.all(np.isfinite(self.B)):
            raise ValueError("ABM attractiveness fields must be finite")
        if np.min(self.ast) <= 0 or np.min(self.B) < 0:
            raise ValueError("ABM requires ast > 0 and initial B >= 0")
        if not np.all(np.isfinite(self.generation_ratio)) or np.min(self.generation_ratio) < 0:
            raise ValueError("ABM requires generation_ratio >= 0")

        scale = config.scaling
        h2 = grid.spacing**2
        self.criminal_density_per_agent = 4.0 * scale.theta * config.time.step / (
            scale.omega * h2
        )
        self.police_density_per_agent = 4.0 * scale.beta * scale.omega / (
            scale.diffusivity * h2
        )
        self.crime_rate_per_event = 4.0 * scale.theta / (scale.omega * h2)
        self.event_attractiveness_increment = scale.theta / scale.omega
        self.generation_rate_scale = scale.omega / scale.theta

        self.criminals = _integer_population(
            initial_rho / self.criminal_density_per_agent, self.rng
        )
        if config.policing.strategy == "none":
            if np.any(initial_pi != 0):
                raise ValueError("initial police density must be zero for strategy: none")
            self.police_agents = np.zeros_like(self.criminals)
        else:
            self.police_agents = _integer_population(
                initial_pi / self.police_density_per_agent, self.rng
            )

        initial_information = config.initial.information
        if config.policing.strategy == "none":
            self.H = np.zeros_like(self.B)
        elif initial_information is None:
            self.H = self.expected_crime_rate.copy()
        else:
            self.H = _grid_field(initial_information, self.x, self.y)
            if not np.all(np.isfinite(self.H)) or np.min(self.H) < 0:
                raise ValueError("initial information must be nonnegative")
        self.time = config.time.start
        self._step_number = 0
        self.total_events = 0
        self.total_generated = 0
        self.initial_criminals = int(self.criminals.sum())
        self.initial_police = int(self.police_agents.sum())

    @property
    def attractiveness(self) -> FloatArray:
        return self.ast + self.B

    @property
    def criminal_density(self) -> FloatArray:
        return self.criminal_density_per_agent * self.criminals

    @property
    def pi(self) -> FloatArray:
        return self.police_density_per_agent * self.police_agents

    @property
    def expected_crime_rate(self) -> FloatArray:
        return self.criminal_density * self.attractiveness * np.exp(-self.pi)

    def _neighbors(self, values: NDArray[Any]) -> tuple[NDArray[Any], ...]:
        if self.config.grid.boundary == "periodic":
            return (
                np.roll(values, -1, axis=1),
                np.roll(values, 1, axis=0),
                np.roll(values, 1, axis=1),
                np.roll(values, -1, axis=0),
            )
        east = np.empty_like(values)
        east[:, :-1], east[:, -1] = values[:, 1:], values[:, -2]
        south = np.empty_like(values)
        south[1:, :], south[0, :] = values[:-1, :], values[1, :]
        west = np.empty_like(values)
        west[:, 1:], west[:, 0] = values[:, :-1], values[:, 1]
        north = np.empty_like(values)
        north[:-1, :], north[-1, :] = values[1:, :], values[-2, :]
        return east, south, west, north

    def _move_population(self, movers: IntArray, preference: FloatArray) -> IntArray:
        neighbor_values = np.stack(self._neighbors(preference))
        totals = np.sum(neighbor_values, axis=0)
        probabilities = np.full_like(neighbor_values, 0.25, dtype=float)
        np.divide(
            neighbor_values,
            totals,
            out=probabilities,
            where=totals[None, :, :] > 0,
        )
        # A multinomial draw can be decomposed into three conditional binomial
        # draws.  NumPy vectorizes those draws over all sites, avoiding a Python
        # loop over the 10,000--40,000 cells in the paper-scale examples.
        east = self.rng.binomial(movers, probabilities[0]).astype(np.int64)
        remainder = movers - east
        denominator = 1.0 - probabilities[0]
        conditional = np.divide(
            probabilities[1], denominator, out=np.zeros_like(denominator), where=denominator > 0
        )
        south = self.rng.binomial(remainder, np.clip(conditional, 0.0, 1.0)).astype(np.int64)
        remainder -= south
        denominator -= probabilities[1]
        conditional = np.divide(
            probabilities[2], denominator, out=np.zeros_like(denominator), where=denominator > 0
        )
        west = self.rng.binomial(remainder, np.clip(conditional, 0.0, 1.0)).astype(np.int64)
        north = remainder - west

        if self.config.grid.boundary == "periodic":
            return (
                np.roll(east, 1, axis=1)
                + np.roll(south, -1, axis=0)
                + np.roll(west, -1, axis=1)
                + np.roll(north, 1, axis=0)
            )

        arrivals = np.zeros_like(movers)
        arrivals[:, 1:] += east[:, :-1]
        arrivals[:, -2] += east[:, -1]
        arrivals[:-1, :] += south[1:, :]
        arrivals[1, :] += south[0, :]
        arrivals[:, :-1] += west[:, 1:]
        arrivals[:, 1] += west[:, 0]
        arrivals[1:, :] += north[:-1, :]
        arrivals[-2, :] += north[-1, :]
        return arrivals

    def step(self) -> ABMStep:
        config = self.config
        dt = config.time.step
        A = self.attractiveness
        pi = self.pi
        criminals_before = self.criminals
        police_before = self.police_agents

        event_probability = burglary_probability_from_fields(A, pi, dt)
        events = self.rng.binomial(criminals_before, event_probability).astype(np.int64)
        criminal_arrivals = self._move_population(criminals_before - events, A)

        sigma = config.model.arrest_probability
        generation_intensity = (
            self.generation_ratio
            * self.generation_rate_scale
            * (1.0 - sigma)
            * np.exp(-pi)
            * dt
        )
        if config.model.generation_distribution == "bernoulli":
            generated = self.rng.binomial(1, -np.expm1(-generation_intensity)).astype(np.int64)
        else:
            generated = self.rng.poisson(generation_intensity).astype(np.int64)

        B_neighbors = self._neighbors(self.B)
        diffusion = sum(B_neighbors) - 4.0 * self.B
        candidate_B = (self.B + 0.25 * config.model.eta * diffusion) * (1.0 - dt)
        candidate_B += self.event_attractiveness_increment * events
        if not np.all(np.isfinite(candidate_B)) or np.min(candidate_B) < -1e-12:
            raise FloatingPointError("attractiveness update lost positivity")

        realized_signal = self.crime_rate_per_event * events
        if config.policing.strategy == "delayed":
            assert config.policing.tau is not None
            busy_probability = -np.expm1(-sigma * criminals_before * event_probability)
            police_staying = self.rng.binomial(police_before, busy_probability).astype(np.int64)
            police_arrivals = self._move_population(police_before - police_staying, self.H)
            next_police = police_staying + police_arrivals
            next_H = self.H + dt / config.policing.tau * (realized_signal - self.H)
        else:
            police_staying = np.zeros_like(police_before)
            next_police = police_before
            next_H = np.zeros_like(self.H)

        next_criminals = criminal_arrivals + generated
        if int(next_police.sum()) != int(police_before.sum()):
            raise RuntimeError("police movement failed to conserve the integer budget")
        if not np.all(np.isfinite(next_H)) or np.min(next_H) < -1e-12:
            raise FloatingPointError("information update lost positivity")

        self.B = np.maximum(candidate_B, 0.0)
        self.criminals = next_criminals
        self.police_agents = next_police
        self.H = np.maximum(next_H, 0.0)
        self._step_number += 1
        self.time = config.time.start + self._step_number * dt
        self.total_events += int(events.sum())
        self.total_generated += int(generated.sum())
        return ABMStep(
            events=events,
            generated=generated,
            police_staying=police_staying,
            realized_crime_rate=realized_signal,
        )


def burglary_probability_from_fields(
    attractiveness: NDArray[np.float64] | float,
    police: NDArray[np.float64] | float,
    dt: float,
) -> NDArray[np.float64]:
    """Exact per-agent event probability with rate ``A * exp(-pi)``."""

    if dt <= 0:
        raise ValueError("dt must be positive")
    A = np.asarray(attractiveness, dtype=float)
    pi = np.asarray(police, dtype=float)
    if np.any(A < 0) or np.any(pi < 0):
        raise ValueError("attractiveness and police density must be nonnegative")
    return -np.expm1(-A * np.exp(-pi) * dt)


def _plain(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {key: _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    return value


def run_abm(config: ABMConfig, *, output_directory: str | Path | None = None) -> ABMResult:
    model = BurglaryABM(config)
    output = Path(output_directory) if output_directory else config.output.root / config.name

    requested = config.output.snapshot_times or (config.time.start, config.time.final)
    snapshot_steps = {
        round((time - config.time.start) / config.time.step) for time in requested
    }
    snapshot_steps.update((0, config.time.steps))
    times: list[float] = []
    attractiveness: list[FloatArray] = []
    density: list[FloatArray] = []
    events: list[IntArray] = []
    police: list[FloatArray] = []
    information: list[FloatArray] = []
    expected_rate: list[FloatArray] = []
    realized_rate: list[FloatArray] = []
    criminal_counts: list[IntArray] = []
    police_counts: list[IntArray] = []

    metric_times: list[float] = []
    mean_A: list[float] = []
    mean_rho: list[float] = []
    mean_pi: list[float] = []
    mean_H: list[float] = []
    mean_expected: list[float] = []
    mean_realized: list[float] = []
    cumulative: list[int] = []

    last_events = np.zeros_like(model.criminals)
    last_realized = np.zeros_like(model.H)

    def save_snapshot() -> None:
        times.append(model.time)
        attractiveness.append(model.attractiveness.copy())
        density.append(model.criminal_density.copy())
        events.append(last_events.copy())
        police.append(model.pi.copy())
        information.append(model.H.copy())
        expected_rate.append(model.expected_crime_rate.copy())
        realized_rate.append(last_realized.copy())
        criminal_counts.append(model.criminals.copy())
        police_counts.append(model.police_agents.copy())

    def save_metrics() -> None:
        metric_times.append(model.time)
        mean_A.append(float(np.mean(model.attractiveness)))
        mean_rho.append(float(np.mean(model.criminal_density)))
        mean_pi.append(float(np.mean(model.pi)))
        mean_H.append(float(np.mean(model.H)))
        mean_expected.append(float(np.mean(model.expected_crime_rate)))
        mean_realized.append(float(np.mean(last_realized)))
        cumulative.append(model.total_events)

    save_snapshot()
    save_metrics()
    for step_number in range(1, config.time.steps + 1):
        step_result = model.step()
        last_events = step_result.events
        last_realized = step_result.realized_crime_rate
        if step_number in snapshot_steps:
            save_snapshot()
        if step_number % config.time.save_every == 0 or step_number == config.time.steps:
            save_metrics()

    output.mkdir(parents=True, exist_ok=True)
    result = ABMResult(
        output_directory=output,
        times=np.asarray(times),
        attractiveness=np.asarray(attractiveness),
        criminal_density=np.asarray(density),
        events=np.asarray(events),
        police=np.asarray(police),
        information=np.asarray(information),
        expected_crime_rate=np.asarray(expected_rate),
        realized_crime_rate=np.asarray(realized_rate),
        criminal_counts=np.asarray(criminal_counts),
        police_counts=np.asarray(police_counts),
        metric_times=np.asarray(metric_times),
        mean_attractiveness=np.asarray(mean_A),
        mean_criminal_density=np.asarray(mean_rho),
        mean_police=np.asarray(mean_pi),
        mean_information=np.asarray(mean_H),
        mean_expected_crime_rate=np.asarray(mean_expected),
        mean_realized_crime_rate=np.asarray(mean_realized),
        cumulative_events=np.asarray(cumulative, dtype=np.int64),
    )
    np.savez_compressed(
        output / "trajectory.npz",
        time=result.times,
        A=result.attractiveness,
        rho=result.criminal_density,
        events=result.events,
        pi=result.police,
        H=result.information,
        expected_S=result.expected_crime_rate,
        realized_S=result.realized_crime_rate,
        criminal_counts=result.criminal_counts,
        police_counts=result.police_counts,
    )
    np.savez_compressed(
        output / "metrics.npz",
        time=result.metric_times,
        mean_A=result.mean_attractiveness,
        mean_rho=result.mean_criminal_density,
        mean_pi=result.mean_police,
        mean_H=result.mean_information,
        mean_expected_S=result.mean_expected_crime_rate,
        mean_realized_S=result.mean_realized_crime_rate,
        cumulative_events=result.cumulative_events,
    )
    (output / "config.resolved.yaml").write_text(
        yaml.safe_dump(_plain(asdict(config)), sort_keys=False), encoding="utf-8"
    )
    summary: dict[str, Any] = {
        "initial_criminals": model.initial_criminals,
        "final_criminals": int(model.criminals.sum()),
        "total_generated": model.total_generated,
        "total_events": model.total_events,
        "initial_police_agents": model.initial_police,
        "final_police_agents": int(model.police_agents.sum()),
        "police_budget_conserved": int(model.police_agents.sum()) == model.initial_police,
        "criminal_density_per_agent": model.criminal_density_per_agent,
        "police_density_per_agent": model.police_density_per_agent,
        "crime_rate_per_event": model.crime_rate_per_event,
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    metadata = {
        "git_revision": git_revision(Path.cwd()),
        "environment": environment_metadata(("numpy", "PyYAML")),
        "model": "discrete burglary ABM with integer criminal and police agents",
        "information_signal": "realized burglary events",
        "arrest_probability": "1 - exp(-Sigma * n * p_crime)",
    }
    (output / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    write_checksums(output)
    return result
