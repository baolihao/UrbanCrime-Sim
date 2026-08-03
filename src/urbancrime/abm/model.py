"""Vectorized Python implementation of the burglary ABM with delayed police."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray
import yaml
from scipy.sparse import coo_matrix, eye
from scipy.sparse.linalg import spsolve

from urbancrime.coefficients import analytic_profile
from urbancrime.config import CoefficientSpec
from urbancrime.provenance import environment_metadata, git_revision, write_checksums

from .config import ABMConfig


@dataclass(frozen=True)
class ABMResult:
    output_directory: Path
    times: NDArray[np.float64]
    attractiveness: NDArray[np.float64]
    criminal_density: NDArray[np.float64]
    events: NDArray[np.int64]
    police: NDArray[np.float64]
    information: NDArray[np.float64]


def _grid_field(spec: CoefficientSpec, x: NDArray[np.float64], y: NDArray[np.float64]) -> NDArray[np.float64]:
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


class BurglaryABM:
    def __init__(self, config: ABMConfig) -> None:
        self.config = config
        grid = config.grid
        coordinates_x = (np.arange(grid.nx) + 0.5) * grid.spacing
        coordinates_y = (np.arange(grid.ny) + 0.5) * grid.spacing
        self.x, self.y = np.meshgrid(coordinates_x, coordinates_y)
        self.ast = _grid_field(config.model.ast, self.x, self.y)
        self.gamma = _grid_field(config.model.gamma, self.x, self.y)
        if np.min(self.ast) <= 0 or np.min(self.gamma) < 0:
            raise ValueError("ABM requires ast > 0 and gamma >= 0")
        self.rng = np.random.default_rng(config.initial.seed)
        self.B = np.full(self.ast.shape, config.initial.dynamic_attractiveness, dtype=float)
        self.criminals = self.rng.poisson(
            config.initial.mean_criminals_per_site, size=self.ast.shape
        ).astype(np.int64)
        self.pi = np.full(self.ast.shape, config.policing.initial_mean_density, dtype=float)
        if config.policing.strategy == "none":
            self.pi.fill(0.0)
        self.H = self.criminals / grid.spacing**2 * (self.ast + self.B) * np.exp(-self.pi)
        self.H = np.maximum(self.H, config.policing.information_floor)
        self.time = config.time.start

    def _neighbors(self, values: NDArray[np.float64]) -> tuple[NDArray[np.float64], ...]:
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

    def _move_criminals(self, movers: NDArray[np.int64], attractiveness: NDArray[np.float64]) -> NDArray[np.int64]:
        neighbor_values = np.stack(self._neighbors(attractiveness))
        probabilities = neighbor_values / np.sum(neighbor_values, axis=0)
        arrivals = np.zeros_like(movers)
        ny, nx = movers.shape
        for row in range(ny):
            for column in range(nx):
                count = int(movers[row, column])
                if count == 0:
                    continue
                allocations = self.rng.multinomial(count, probabilities[:, row, column])
                targets = [
                    (row, column + 1),
                    (row - 1, column),
                    (row, column - 1),
                    (row + 1, column),
                ]
                for amount, (target_row, target_column) in zip(allocations, targets, strict=True):
                    if self.config.grid.boundary == "periodic":
                        target_row %= ny
                        target_column %= nx
                    else:
                        if target_row < 0:
                            target_row = 1
                        elif target_row >= ny:
                            target_row = ny - 2
                        if target_column < 0:
                            target_column = 1
                        elif target_column >= nx:
                            target_column = nx - 2
                    arrivals[target_row, target_column] += amount
        return arrivals

    def _update_police(self, crime_signal: NDArray[np.float64]) -> None:
        police = self.config.policing
        if police.strategy != "delayed":
            return
        assert police.tau is not None
        dt = self.config.time.step
        self.H = (police.tau * self.H + dt * crime_signal) / (police.tau + dt)
        H = np.maximum(self.H, police.information_floor)
        substep = dt / police.diffusion_substeps
        spacing2 = self.config.grid.spacing**2
        for _ in range(police.diffusion_substeps):
            operator = self._police_transport_operator(H, spacing2)
            system = eye(self.pi.size, format="csr") - substep * police.diffusivity * operator
            candidate = np.asarray(spsolve(system, self.pi.ravel()), dtype=float).reshape(self.pi.shape)
            if np.min(candidate) < -1.0e-10:
                raise FloatingPointError("implicit police transport lost positivity")
            candidate[candidate < 0.0] = 0.0
            self.pi = candidate

    def _police_transport_operator(self, H: NDArray[np.float64], spacing2: float) -> Any:
        """Assemble a conservative graph operator for div(H^2 grad(pi/H^2))."""

        ny, nx = H.shape
        rows: list[int] = []
        columns: list[int] = []
        values: list[float] = []

        def add_edge(first: int, second: int) -> None:
            h_first = H.ravel()[first]
            h_second = H.ravel()[second]
            coefficient = 0.5 * (h_first**2 + h_second**2) / spacing2
            rows.extend((first, first, second, second))
            columns.extend((second, first, first, second))
            values.extend(
                (
                    coefficient / h_second**2,
                    -coefficient / h_first**2,
                    coefficient / h_first**2,
                    -coefficient / h_second**2,
                )
            )

        for row in range(ny):
            for column in range(nx - 1):
                add_edge(row * nx + column, row * nx + column + 1)
        for row in range(ny - 1):
            for column in range(nx):
                add_edge(row * nx + column, (row + 1) * nx + column)
        if self.config.grid.boundary == "periodic":
            for row in range(ny):
                add_edge(row * nx + nx - 1, row * nx)
            for column in range(nx):
                add_edge((ny - 1) * nx + column, column)
        return coo_matrix((values, (rows, columns)), shape=(nx * ny, nx * ny)).tocsr()

    def step(self) -> NDArray[np.int64]:
        config = self.config
        dt = config.time.step
        attractiveness = self.ast + self.B
        burglary_probability = burglary_probability_from_fields(attractiveness, self.pi, dt)
        events = self.rng.binomial(self.criminals, burglary_probability)
        movers = self.criminals - events
        arrivals = self._move_criminals(movers, attractiveness)
        if config.model.generation_distribution == "bernoulli":
            generation_probability = 1.0 - np.exp(-self.gamma * dt)
            generated = self.rng.binomial(1, generation_probability)
        else:
            generated = self.rng.poisson(self.gamma * dt)
        neighbors = self._neighbors(self.B)
        laplacian = (sum(neighbors) - 4.0 * self.B) / config.grid.spacing**2
        decay = 1.0 - config.model.omega * dt
        if decay < 0:
            raise ValueError("ABM requires omega * time.step <= 1")
        self.B = (
            self.B + config.model.eta * config.grid.spacing**2 / 4.0 * laplacian
        ) * decay + config.model.theta * events
        self.criminals = arrivals + generated
        crime_signal = (
            self.criminals / config.grid.spacing**2 * (self.ast + self.B) * np.exp(-self.pi)
        )
        self._update_police(crime_signal)
        self.time += dt
        return events


def burglary_probability_from_fields(
    attractiveness: NDArray[np.float64] | float,
    police: NDArray[np.float64] | float,
    dt: float,
) -> NDArray[np.float64]:
    """Per-agent event probability whose small-dt limit is ``A*exp(-pi)``."""

    if dt <= 0:
        raise ValueError("dt must be positive")
    A = np.asarray(attractiveness, dtype=float)
    pi = np.asarray(police, dtype=float)
    if np.any(A < 0) or np.any(pi < 0):
        raise ValueError("attractiveness and police density must be nonnegative")
    return 1.0 - np.exp(-A * np.exp(-pi) * dt)


def run_abm(config: ABMConfig, *, output_directory: str | Path | None = None) -> ABMResult:
    model = BurglaryABM(config)
    output = Path(output_directory) if output_directory else config.output.root / config.name
    times: list[float] = []
    attractiveness: list[np.ndarray] = []
    density: list[np.ndarray] = []
    events: list[np.ndarray] = []
    police: list[np.ndarray] = []
    information: list[np.ndarray] = []

    def save(event_values: NDArray[np.int64]) -> None:
        times.append(model.time)
        attractiveness.append((model.ast + model.B).copy())
        density.append(model.criminals.astype(float) / config.grid.spacing**2)
        events.append(event_values.copy())
        police.append(model.pi.copy())
        information.append(model.H.copy())

    save(np.zeros_like(model.criminals))
    for step in range(1, config.time.steps + 1):
        event_values = model.step()
        if step % config.time.save_every == 0:
            save(event_values)
    output.mkdir(parents=True, exist_ok=True)
    result = ABMResult(
        output_directory=output,
        times=np.asarray(times),
        attractiveness=np.asarray(attractiveness),
        criminal_density=np.asarray(density),
        events=np.asarray(events),
        police=np.asarray(police),
        information=np.asarray(information),
    )
    np.savez_compressed(
        output / "trajectory.npz",
        time=result.times,
        A=result.attractiveness,
        rho=result.criminal_density,
        events=result.events,
        pi=result.police,
        H=result.information,
    )
    resolved = asdict(config)
    resolved["output"]["root"] = str(resolved["output"]["root"])
    for name in ("ast", "gamma"):
        if resolved["model"][name]["path"] is not None:
            resolved["model"][name]["path"] = str(resolved["model"][name]["path"])
    (output / "config.resolved.yaml").write_text(
        yaml.safe_dump(resolved, sort_keys=False), encoding="utf-8"
    )
    summary: dict[str, Any] = {
        "total_events": int(result.events.sum()),
        "final_criminals": int(model.criminals.sum()),
        "initial_police_mass": float(result.police[0].sum() * config.grid.spacing**2),
        "final_police_mass": float(result.police[-1].sum() * config.grid.spacing**2),
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    metadata = {
        "git_revision": git_revision(Path.cwd()),
        "environment": environment_metadata(("numpy", "scipy", "PyYAML")),
    }
    (output / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    write_checksums(output)
    return result
