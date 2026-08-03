"""Executable DOLFINx runtime for the canonical heterogeneous-flux model."""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any, Iterable

import numpy as np
from numpy.typing import NDArray

from .analysis.metrics import crime_intensity
from .coefficients import build_coefficient, build_model_coefficients
from .config import (
    CoefficientSpec,
    DelayedPoliceConfig,
    NoPoliceConfig,
    OptimalPoliceConfig,
    PoliceConfig,
    PrescribedPoliceConfig,
    RunConfig,
    write_resolved_config,
)
from .mesh.integration import LumpedMeasure
from .mesh.polygon import create_polygon
from .mesh.rectangle import create_rectangle
from .models.burglary import burglary_residual
from .policing.delayed import delayed_police_residual, mass_lumped_h_update
from .policing.optimal import optimal_police_density
from .provenance import environment_metadata, git_revision, write_checksums
from .solvers.implicit import ImplicitSolver
from .state import SimulationState


@dataclass(frozen=True)
class PDEResult:
    output_directory: Path
    final_state: SimulationState
    times: NDArray[np.float64]
    metrics: dict[str, NDArray[np.float64]]


def _polygon_rings(path: Path) -> tuple[np.ndarray, tuple[np.ndarray, ...]]:
    if path.suffix == ".npy":
        return np.asarray(np.load(path), dtype=float), ()
    if path.suffix == ".npz":
        with np.load(path) as archive:
            exterior = np.asarray(archive["exterior"], dtype=float)
            holes = tuple(
                np.asarray(archive[name], dtype=float)
                for name in sorted(archive.files)
                if name.startswith("hole_")
            )
        return exterior, holes
    if path.suffix == ".json":
        raw = json.loads(path.read_text(encoding="utf-8"))
        return np.asarray(raw["exterior"], dtype=float), tuple(
            np.asarray(ring, dtype=float) for ring in raw.get("holes", [])
        )
    raise ValueError("polygon boundary must be .npy, .npz, or .json")


def build_domain(config: RunConfig, comm: Any = None) -> Any:
    """Create the configured rectangle or polygon mesh."""

    if config.domain.kind == "rectangle":
        if config.mesh.nx is None or config.mesh.ny is None:
            raise ValueError("rectangle mesh requires nx and ny")
        return create_rectangle(
            lower_left=config.domain.parameters["lower_left"],
            upper_right=config.domain.parameters["upper_right"],
            nx=config.mesh.nx,
            ny=config.mesh.ny,
            cell_type=config.mesh.cell_type,
            comm=comm,
        )
    boundary_file = Path(config.domain.parameters["boundary_file"])
    exterior, holes = _polygon_rings(boundary_file)
    scale_factor = float(config.domain.parameters.get("scale_factor", 1.0))
    if scale_factor <= 0:
        raise ValueError("polygon domain scale_factor must be positive")
    if scale_factor != 1.0:
        exterior = exterior * scale_factor
        holes = tuple(ring * scale_factor for ring in holes)
    if config.mesh.characteristic_length is None:
        raise ValueError("polygon mesh requires characteristic_length")
    mesh, _, _ = create_polygon(
        exterior,
        holes=holes,
        characteristic_length=config.mesh.characteristic_length,
        boundary_length=config.mesh.boundary_length,
        interior_length=config.mesh.interior_length,
        transition_width=config.mesh.transition_width,
        comm=comm,
    )
    return mesh


class PDERunner:
    """Run one configured PDE simulation and expose immutable branch states."""

    def __init__(
        self,
        config: RunConfig,
        *,
        output_directory: str | Path | None = None,
        initial_state: SimulationState | None = None,
        comm: Any = None,
        write_output: bool = True,
    ) -> None:
        import ufl
        from dolfinx import fem
        from mpi4py import MPI
        from petsc4py import PETSc

        self.config = config
        self.comm = MPI.COMM_WORLD if comm is None else comm
        self.mesh = build_domain(config, self.comm)
        self.V = fem.functionspace(
            self.mesh, ("Lagrange", config.mesh.polynomial_degree)
        )
        self.measure = LumpedMeasure.from_function_space(self.V)
        self.coefficients = build_model_coefficients(config.model, self.V, self.mesh)
        self.dt = fem.Constant(self.mesh, PETSc.ScalarType(config.time.step))
        self.zero = fem.Constant(self.mesh, PETSc.ScalarType(0.0))
        self.fields = {
            name: fem.Function(self.V, name=name)
            for name in self._active_field_names(config.policing)
        }
        self.previous = {
            name: fem.Function(self.V, name=f"{name}_previous") for name in self.fields
        }
        self._tau = (
            build_coefficient(config.policing.tau, self.V, self.mesh)
            if isinstance(config.policing, DelayedPoliceConfig)
            else None
        )
        self._external_pi = (
            fem.Function(self.V, name="pi")
            if isinstance(config.policing, (PrescribedPoliceConfig, OptimalPoliceConfig))
            else None
        )
        if self._external_pi is not None:
            self._external_pi.x.array[:] = 0.0
            self._external_pi.x.scatter_forward()
        self._prescribed_activated = False
        self.time = config.time.start if initial_state is None else initial_state.time
        self._initialize_fields()
        if initial_state is not None:
            self.restore(initial_state)
        for name in self.fields:
            self._copy(self.fields[name], self.previous[name])
        self._validate_state()
        self._mixed = None
        self._mixed_previous = None
        self._mixed_maps: dict[str, NDArray[np.int32]] = {}
        self._implicit_solver: ImplicitSolver | None = None
        self._linear_problems: dict[str, Any] = {}
        if config.solver.method == "implicit":
            if (
                isinstance(config.policing, DelayedPoliceConfig)
                and config.policing.h_discretization != "consistent"
            ):
                raise ValueError("implicit delayed solve requires consistent H discretization")
            self._build_implicit_solver()
        else:
            self._build_partitioned_solver()
        root = Path(config.output.root) / config.name
        self.output_directory = Path(output_directory) if output_directory else root
        self.write_output = write_output
        self._writers: dict[str, Any] = {}
        self._times: list[float] = []
        self._history = {"crime_average": [], "crime_maximum": [], "police_mass": []}
        self._ufl = ufl

    @staticmethod
    def _active_field_names(policing: PoliceConfig) -> tuple[str, ...]:
        return ("A", "rho", "pi", "H") if isinstance(policing, DelayedPoliceConfig) else ("A", "rho")

    def _coefficient_array(self, coefficient: Any) -> NDArray[np.float64]:
        if hasattr(coefficient, "x"):
            return np.asarray(coefficient.x.array, dtype=float)
        value = float(np.asarray(coefficient.value).reshape(-1)[0])
        return np.full(self.fields["A"].x.array.shape, value, dtype=float)

    def _noise(self, seed: int, stream: int, standard_deviation: float) -> NDArray[np.float64]:
        size_local = self.V.dofmap.index_map.size_local
        ghosts = self.V.dofmap.index_map.num_ghosts
        local = np.arange(size_local + ghosts, dtype=np.int32)
        global_indices = self.V.dofmap.index_map.local_to_global(local)
        global_size = self.V.dofmap.index_map.size_global
        rng = np.random.default_rng(np.random.SeedSequence([seed, stream]))
        global_values = rng.normal(0.0, standard_deviation, global_size)
        return global_values[global_indices]

    def _uniform(self, seed: int, stream: int) -> NDArray[np.float64]:
        size_local = self.V.dofmap.index_map.size_local
        ghosts = self.V.dofmap.index_map.num_ghosts
        local = np.arange(size_local + ghosts, dtype=np.int32)
        global_indices = self.V.dofmap.index_map.local_to_global(local)
        rng = np.random.default_rng(np.random.SeedSequence([seed, stream]))
        return rng.random(self.V.dofmap.index_map.size_global)[global_indices]

    def _initialize_fields(self) -> None:
        initial = self.config.initial_conditions
        ast = self._coefficient_array(self.coefficients.ast)
        q = self._coefficient_array(self.coefficients.q)
        noise = initial.noise or {}
        seed = int(noise.get("seed", 0))
        sparsity = float(noise.get("sparsity", 1.0))
        if not 0.0 <= sparsity <= 1.0:
            raise ValueError("initial noise sparsity must lie in [0, 1]")

        def assign(name: str, specification: dict[str, Any] | Any, stream: int) -> None:
            raw = dict(specification)
            kind = str(raw["kind"])
            field_sparsity = float(noise.get(f"{name}_sparsity", sparsity))
            if not 0.0 <= field_sparsity <= 1.0:
                raise ValueError(f"initial {name} noise sparsity must lie in [0, 1]")
            baseline: NDArray[np.float64]
            perturbation: NDArray[np.float64] | None = None
            if kind == "constant":
                baseline = np.full(self.fields[name].x.array.shape, float(raw["value"]))
            elif kind == "constant_plus_noise":
                sd = float(noise[f"{name}_standard_deviation"])
                baseline = np.full(self.fields[name].x.array.shape, float(raw["value"]))
                perturbation = self._noise(seed, stream, sd)
            elif kind == "equilibrium_plus_noise" and name == "A":
                sd = float(noise["A_standard_deviation"])
                baseline = ast + float(raw["dynamic_attractiveness"])
                perturbation = self._noise(seed, stream, sd)
            elif kind == "ast_plus_q_plus_noise" and name == "A":
                sd = float(noise["A_standard_deviation"])
                baseline = ast + q
                perturbation = self._noise(seed, stream, sd)
            elif kind == "uniform_random":
                lower, upper = float(raw["low"]), float(raw["high"])
                if upper <= lower:
                    raise ValueError("uniform initial condition requires high > low")
                baseline = lower + (upper - lower) * self._uniform(seed, stream)
            elif kind == "analytic":
                coefficient = build_coefficient(
                    CoefficientSpec.from_mapping(raw), self.V, self.mesh
                )
                baseline = self._coefficient_array(coefficient).copy()
            else:
                raise ValueError(f"unsupported initial condition for {name}: {kind}")
            values = baseline
            if perturbation is not None:
                if field_sparsity < 1.0:
                    mask = self._uniform(seed, stream + 1000) < field_sparsity
                    perturbation = np.where(mask, perturbation, 0.0)
                values = baseline + perturbation
            self.fields[name].x.array[:] = values
            self.fields[name].x.scatter_forward()

        assign("A", initial.A, 1)
        assign("rho", initial.rho, 2)
        if "pi" in self.fields:
            pi_spec = initial.pi or {"kind": "constant", "value": self.config.policing.initial_mean_density}
            assign("pi", pi_spec, 3)
            h_spec = dict(initial.H or {"kind": "instantaneous_crime_signal"})
            if h_spec["kind"] == "instantaneous_crime_signal":
                self.fields["H"].x.array[:] = crime_intensity(
                    self.fields["A"].x.array,
                    self.fields["rho"].x.array,
                    self.fields["pi"].x.array,
                )
                self.fields["H"].x.scatter_forward()
            else:
                assign("H", h_spec, 4)

    @staticmethod
    def _copy(source: Any, destination: Any) -> None:
        destination.x.array[:] = source.x.array
        destination.x.scatter_forward()

    def restore(self, state: SimulationState) -> None:
        for name, values in state.fields.items():
            if name not in self.fields:
                continue
            if values.shape != self.fields[name].x.array.shape:
                raise ValueError(f"state field {name} is incompatible with this mesh")
            self.fields[name].x.array[:] = values
            self.fields[name].x.scatter_forward()
        self.time = state.time

    def capture(self) -> SimulationState:
        return SimulationState.capture(
            self.time, **{name: field.x.array for name, field in self.fields.items()}
        )

    @property
    def pi(self) -> Any:
        if "pi" in self.fields:
            return self.fields["pi"]
        return self._external_pi if self._external_pi is not None else self.zero

    def _build_implicit_solver(self) -> None:
        import ufl
        from basix.ufl import element, mixed_element
        from dolfinx import fem

        names = tuple(self.fields)
        scalar_element = element("Lagrange", self.mesh.basix_cell(), self.config.mesh.polynomial_degree)
        W = fem.functionspace(self.mesh, mixed_element([scalar_element] * len(names)))
        self._mixed = fem.Function(W, name="solution")
        self._mixed_previous = fem.Function(W, name="solution_previous")
        components = ufl.split(self._mixed)
        old_components = ufl.split(self._mixed_previous)
        tests = ufl.TestFunctions(W)
        for index, name in enumerate(names):
            collapsed, mapping = W.sub(index).collapse()
            collapsed_coordinates = collapsed.tabulate_dof_coordinates()
            scalar_coordinates = self.V.tabulate_dof_coordinates()
            if (
                collapsed_coordinates.shape != scalar_coordinates.shape
                or not np.allclose(
                    collapsed_coordinates, scalar_coordinates, rtol=0.0, atol=1.0e-12
                )
                or len(mapping) != self.fields[name].x.array.size
            ):
                raise RuntimeError(
                    f"mixed subspace {name} is not dof-compatible with the scalar space"
                )
            self._mixed_maps[name] = np.asarray(mapping, dtype=np.int32)
        self._sync_fields_to_mixed(previous=False)
        self._sync_fields_to_mixed(previous=True)
        index = {name: position for position, name in enumerate(names)}
        pi = components[index["pi"]] if "pi" in index else self.pi
        residual = burglary_residual(
            A=components[index["A"]],
            A_previous=old_components[index["A"]],
            rho=components[index["rho"]],
            rho_previous=old_components[index["rho"]],
            test_A=tests[index["A"]],
            test_rho=tests[index["rho"]],
            dt=self.dt,
            coefficients=self.coefficients,
            pi=pi,
        )
        if "pi" in index:
            residual += delayed_police_residual(
                pi=components[index["pi"]],
                pi_previous=old_components[index["pi"]],
                H=components[index["H"]],
                H_previous=old_components[index["H"]],
                test_pi=tests[index["pi"]],
                test_H=tests[index["H"]],
                A=components[index["A"]],
                rho=components[index["rho"]],
                tau=self._tau,
                dt=self.dt,
            )
        self._implicit_solver = ImplicitSolver(
            mesh=self.mesh,
            residual=residual,
            solution=self._mixed,
            relative_tolerance=self.config.solver.relative_tolerance,
            absolute_tolerance=self.config.solver.absolute_tolerance,
            max_iterations=self.config.solver.max_iterations,
            convergence_criterion=self.config.solver.convergence_criterion,
            petsc_options=self.config.solver.petsc_options,
        )

    def _sync_fields_to_mixed(self, *, previous: bool) -> None:
        target = self._mixed_previous if previous else self._mixed
        source = self.previous if previous else self.fields
        if target is None:
            return
        for name, mapping in self._mixed_maps.items():
            target.x.array[mapping] = source[name].x.array
        target.x.scatter_forward()

    def _sync_mixed_to_fields(self) -> None:
        if self._mixed is None:
            return
        for name, mapping in self._mixed_maps.items():
            self.fields[name].x.array[:] = self._mixed.x.array[mapping]
            self.fields[name].x.scatter_forward()

    def _build_partitioned_solver(self) -> None:
        import ufl
        from dolfinx.fem.petsc import LinearProblem

        V = self.V
        options = dict(self.config.solver.petsc_options)
        dx = ufl.Measure("dx", domain=self.mesh)
        A_trial, a = ufl.TrialFunction(V), ufl.TestFunction(V)
        attenuation = ufl.exp(-self.pi)
        eta, ast, q = self.coefficients.eta, self.coefficients.ast, self.coefficients.q
        dt = self.dt
        A_form = (
            ufl.inner(A_trial, a) * dx
            + dt * ufl.inner(eta * ufl.grad(A_trial), ufl.grad(a)) * dx
            + dt * ufl.inner(A_trial - self.fields["rho"] * A_trial * attenuation, a) * dx
        )
        A_rhs = (
            ufl.inner(self.previous["A"], a) * dx
            + dt * ufl.inner(eta * ufl.grad(ast), ufl.grad(a)) * dx
            + dt * ufl.inner(ast, a) * dx
        )
        self._linear_problems["A"] = LinearProblem(
            A_form, A_rhs, u=self.fields["A"], petsc_options=options
        )
        rho_trial, r = ufl.TrialFunction(V), ufl.TestFunction(V)
        rho_form = (
            ufl.inner(rho_trial, r) * dx
            + dt * ufl.inner(ufl.grad(rho_trial), ufl.grad(r)) * dx
            - 2.0
            * dt
            * ufl.inner(
                rho_trial / self.fields["A"] * ufl.grad(self.fields["A"]), ufl.grad(r)
            )
            * dx
            + dt * ufl.inner(rho_trial * self.fields["A"] * attenuation, r) * dx
        )
        rho_rhs = (
            ufl.inner(self.previous["rho"], r) * dx
            + dt * ufl.inner(q * attenuation, r) * dx
        )
        self._linear_problems["rho"] = LinearProblem(
            rho_form, rho_rhs, u=self.fields["rho"], petsc_options=options
        )
        if isinstance(self.config.policing, DelayedPoliceConfig):
            if self.config.policing.h_discretization == "consistent":
                H_trial, h = ufl.TrialFunction(V), ufl.TestFunction(V)
                H_form = ufl.inner(H_trial, h) * dx + dt * ufl.inner(
                    H_trial / self._tau, h
                ) * dx
                H_rhs = ufl.inner(self.previous["H"], h) * dx + dt * ufl.inner(
                    self.fields["rho"]
                    * self.fields["A"]
                    * ufl.exp(-self.fields["pi"])
                    / self._tau,
                    h,
                ) * dx
                self._linear_problems["H"] = LinearProblem(
                    H_form, H_rhs, u=self.fields["H"], petsc_options=options
                )
            pi_trial, p = ufl.TrialFunction(V), ufl.TestFunction(V)
            pi_form = (
                ufl.inner(pi_trial, p) * dx
                + dt * ufl.inner(ufl.grad(pi_trial), ufl.grad(p)) * dx
                - 2.0
                * dt
                * ufl.inner(
                    pi_trial / self.fields["H"] * ufl.grad(self.fields["H"]), ufl.grad(p)
                )
                * dx
            )
            pi_rhs = ufl.inner(self.previous["pi"], p) * dx
            self._linear_problems["pi"] = LinearProblem(
                pi_form, pi_rhs, u=self.fields["pi"], petsc_options=options
            )

    def _relative_change(self, new: NDArray[np.float64], old: NDArray[np.float64]) -> float:
        from mpi4py import MPI

        owned = self.V.dofmap.index_map.size_local
        numerator = self.comm.allreduce(float(np.dot(new[:owned] - old[:owned], new[:owned] - old[:owned])), op=MPI.SUM)
        denominator = self.comm.allreduce(float(np.dot(old[:owned], old[:owned])), op=MPI.SUM)
        return float(np.sqrt(numerator) / max(np.sqrt(denominator), self.config.solver.absolute_tolerance))

    def _solve_partition_block(self, name: str) -> float:
        old = self.fields[name].x.array.copy()
        if name == "H" and isinstance(self.config.policing, DelayedPoliceConfig) and self.config.policing.h_discretization == "mass_lumped":
            self.fields["H"].x.array[:] = mass_lumped_h_update(
                self.previous["H"].x.array,
                self.fields["A"].x.array,
                self.fields["rho"].x.array,
                self.fields["pi"].x.array,
                self._coefficient_array(self._tau),
                self.config.time.step,
            )
            self.fields["H"].x.scatter_forward()
        else:
            self._linear_problems[name].solve()
            self.fields[name].x.scatter_forward()
        return self._relative_change(self.fields[name].x.array, old)

    def _partitioned_step(self) -> int:
        order = self.config.solver.block_order
        if set(order) != set(self.fields):
            raise ValueError(f"block_order {order} does not match active fields {tuple(self.fields)}")
        for iteration in range(1, self.config.solver.max_iterations + 1):
            if isinstance(self.config.policing, OptimalPoliceConfig) and self.config.policing.update_each_iteration:
                self._update_optimal_policy()
            changes = [self._solve_partition_block(name) for name in order]
            if max(changes) <= self.config.solver.relative_tolerance:
                for name in self.fields:
                    self._copy(self.fields[name], self.previous[name])
                return iteration
        raise RuntimeError(f"partitioned solve failed to converge; last changes={changes}")

    def _set_prescribed_policy(self) -> None:
        config = self.config.policing
        if not isinstance(config, PrescribedPoliceConfig) or self._external_pi is None:
            return
        if config.profile != "snapshot_tanh":
            raise ValueError(f"unknown prescribed profile: {config.profile}")
        kappa = float(config.parameters["kappa"])
        cutoff = float(config.parameters["cutoff"])
        attenuation = 0.5 * (1.0 - np.tanh(kappa * (self.fields["A"].x.array - cutoff)))
        raw = -np.log(np.clip(attenuation, config.positivity_floor, 1.0))
        total = config.budget.resolve_total(self.measure.domain_area)
        raw_mass = self.measure.integral(raw)
        if raw_mass <= config.positivity_floor:
            self._external_pi.x.array[:] = total / self.measure.domain_area
        else:
            self._external_pi.x.array[:] = raw * total / raw_mass
        self._external_pi.x.scatter_forward()
        self._prescribed_activated = True

    def _update_optimal_policy(self) -> None:
        config = self.config.policing
        if not isinstance(config, OptimalPoliceConfig) or self._external_pi is None:
            return
        owned = self.V.dofmap.index_map.size_local
        local_intensity = (self.fields["rho"].x.array * self.fields["A"].x.array)[:owned].copy()
        local_weights = self.measure.weights[:owned].copy()
        intensities = self.comm.gather(local_intensity, root=0)
        weights = self.comm.gather(local_weights, root=0)
        pieces = None
        if self.comm.rank == 0:
            assert intensities is not None and weights is not None
            sizes = [values.size for values in intensities]
            solution = optimal_police_density(
                np.concatenate(intensities),
                weights=np.concatenate(weights),
                budget=config.budget.resolve_total(self.measure.domain_area),
                positivity_floor=config.positivity_floor,
                root_tolerance=config.root_tolerance,
                max_iterations=config.root_max_iterations,
            )
            offsets = np.cumsum([0, *sizes])
            pieces = [solution[offsets[i] : offsets[i + 1]] for i in range(len(sizes))]
        local_solution = self.comm.scatter(pieces, root=0)
        self._external_pi.x.array[:owned] = local_solution
        self._external_pi.x.scatter_forward()

    def _prepare_policy(self, next_time: float) -> None:
        config = self.config.policing
        if isinstance(config, PrescribedPoliceConfig):
            if next_time >= config.activation_time and not self._prescribed_activated:
                self._set_prescribed_policy()
        elif isinstance(config, OptimalPoliceConfig):
            if next_time >= config.activation_time:
                self._update_optimal_policy()

    def step(self) -> int:
        next_time = self.time + self.config.time.step
        self._prepare_policy(next_time)
        if self.config.solver.method == "implicit":
            self._sync_fields_to_mixed(previous=False)
            assert self._implicit_solver is not None and self._mixed is not None
            iterations, _ = self._implicit_solver.solve(self._mixed)
            self._sync_mixed_to_fields()
            for name in self.fields:
                self._copy(self.fields[name], self.previous[name])
            self._sync_fields_to_mixed(previous=True)
        else:
            iterations = self._partitioned_step()
        self.time = next_time
        self._validate_state()
        return iterations

    def _validate_state(self) -> None:
        from mpi4py import MPI

        checked = dict(self.fields)
        if self._external_pi is not None:
            checked["pi"] = self._external_pi
        minima: dict[str, float] = {}
        for name, field in checked.items():
            values = field.x.array
            finite = self.comm.allreduce(bool(np.all(np.isfinite(values))), op=MPI.LAND)
            if not finite:
                raise FloatingPointError(f"field {name} contains non-finite values")
            owned = self.V.dofmap.index_map.size_local
            local_minimum = float(np.min(values[:owned])) if owned else np.inf
            minima[name] = self.comm.allreduce(local_minimum, op=MPI.MIN)
        tolerance = self.config.solver.state_tolerance
        for name in ("A", "H"):
            if name in minima and minima[name] <= tolerance:
                raise FloatingPointError(
                    f"{name} lost strict positivity: global minimum {minima[name]:.6e}"
                )
        for name in ("rho", "pi"):
            if name in minima and minima[name] < -tolerance:
                raise FloatingPointError(
                    f"{name} became negative: global minimum {minima[name]:.6e}"
                )

    def _global_values(self, field: Any) -> NDArray[np.float64] | None:
        owned = self.V.dofmap.index_map.size_local
        indices = self.V.dofmap.index_map.local_to_global(np.arange(owned, dtype=np.int32))
        gathered_indices = self.comm.gather(indices, root=0)
        gathered_values = self.comm.gather(field.x.array[:owned].copy(), root=0)
        if self.comm.rank != 0:
            return None
        result = np.empty(self.V.dofmap.index_map.size_global, dtype=float)
        assert gathered_indices is not None and gathered_values is not None
        for index, values in zip(gathered_indices, gathered_values, strict=True):
            result[index] = values
        return result

    def _global_mesh_data(
        self,
    ) -> tuple[NDArray[np.float64] | None, NDArray[np.int64] | None]:
        """Gather P1 dof coordinates and cell connectivity for portable snapshots."""

        owned = self.V.dofmap.index_map.size_local
        local = np.arange(owned, dtype=np.int32)
        global_indices = self.V.dofmap.index_map.local_to_global(local)
        coordinates = self.V.tabulate_dof_coordinates()[:owned, :2].copy()
        gathered_indices = self.comm.gather(global_indices, root=0)
        gathered_coordinates = self.comm.gather(coordinates, root=0)

        all_local = np.arange(
            owned + self.V.dofmap.index_map.num_ghosts, dtype=np.int32
        )
        local_to_global = self.V.dofmap.index_map.local_to_global(all_local)
        cell_map = self.mesh.topology.index_map(self.mesh.topology.dim)
        cell_dofs = np.asarray(
            [
                local_to_global[self.V.dofmap.cell_dofs(cell)]
                for cell in range(cell_map.size_local)
            ],
            dtype=np.int64,
        )
        gathered_cells = self.comm.gather(cell_dofs, root=0)
        if self.comm.rank != 0:
            return None, None

        points = np.empty((self.V.dofmap.index_map.size_global, 2), dtype=float)
        assert gathered_indices is not None and gathered_coordinates is not None
        for indices, values in zip(
            gathered_indices, gathered_coordinates, strict=True
        ):
            points[indices] = values
        assert gathered_cells is not None
        cells = np.vstack(gathered_cells)
        return points, cells

    def _record_metrics(self) -> None:
        from mpi4py import MPI

        pi_values = self._coefficient_array(self.pi)
        intensity = crime_intensity(
            self.fields["A"].x.array, self.fields["rho"].x.array, pi_values
        )
        self._times.append(self.time)
        self._history["crime_average"].append(self.measure.integral(intensity) / self.measure.domain_area)
        owned = self.V.dofmap.index_map.size_local
        self._history["crime_maximum"].append(
            self.comm.allreduce(float(np.max(intensity[:owned])), op=MPI.MAX)
        )
        self._history["police_mass"].append(self.measure.integral(pi_values))

    def _open_output(self) -> None:
        if not self.write_output:
            return
        from dolfinx.io import XDMFFile

        if self.comm.rank == 0:
            self.output_directory.mkdir(parents=True, exist_ok=True)
            write_resolved_config(self.config, self.output_directory / "config.resolved.yaml")
            metadata = {
                "git_revision": git_revision(Path.cwd()),
                "environment": environment_metadata(
                    ("dolfinx", "fenics-ufl", "petsc4py", "mpi4py", "numpy", "scipy")
                ),
                "domain_area": self.measure.domain_area,
                "global_dofs": self.V.dofmap.index_map.size_global,
                "mpi_size": self.comm.size,
            }
            (self.output_directory / "metadata.json").write_text(
                json.dumps(metadata, indent=2), encoding="utf-8"
            )
        self.comm.barrier()
        output_fields = list(self.config.output.fields)
        for name in output_fields:
            if name == "pi" and name not in self.fields and self._external_pi is not None:
                field = self._external_pi
            elif name in self.fields:
                field = self.fields[name]
            else:
                continue
            writer = XDMFFile(self.comm, self.output_directory / f"{name}.xdmf", "w")
            writer.write_mesh(self.mesh)
            self._writers[name] = writer
            writer.write_function(field, self.time)

    def _write_fields(self) -> None:
        for name, writer in self._writers.items():
            field = self._external_pi if name == "pi" and name not in self.fields else self.fields[name]
            writer.write_function(field, self.time)

    def _write_snapshot(self, label: str) -> None:
        values = {name: self._global_values(field) for name, field in self.fields.items()}
        if self._external_pi is not None:
            values["pi"] = self._global_values(self._external_pi)
        points, cells = self._global_mesh_data()
        if self.comm.rank == 0:
            np.savez_compressed(
                self.output_directory / f"snapshot_{label}.npz",
                time=np.asarray(self.time),
                points=points,
                cells=cells,
                **{name: value for name, value in values.items() if value is not None},
            )

    def _close_output(self) -> None:
        for writer in self._writers.values():
            writer.close()
        if not self.write_output:
            return
        self._write_snapshot("final")
        if self.comm.rank == 0:
            np.savez_compressed(
                self.output_directory / "timeseries.npz",
                time=np.asarray(self._times),
                **{name: np.asarray(values) for name, values in self._history.items()},
            )
            write_checksums(self.output_directory)

    def run(
        self,
        *,
        final_time: float | None = None,
        snapshot_times: Iterable[float] = (),
    ) -> PDEResult:
        stop = self.config.time.final if final_time is None else float(final_time)
        if stop < self.time:
            raise ValueError("final time precedes current state")
        raw_steps = (stop - self.time) / self.config.time.step
        steps = round(raw_steps)
        if abs(raw_steps - steps) > 1.0e-10:
            raise ValueError("run interval must be an integer multiple of time.step")
        snapshots = {round(float(value), 12) for value in snapshot_times}
        self._prepare_policy(self.time)
        self._open_output()
        self._record_metrics()
        if round(self.time, 12) in snapshots and self.write_output:
            self._write_snapshot(f"t{self.time:g}")
        for step_index in range(1, steps + 1):
            self.step()
            self._record_metrics()
            if self.write_output and step_index % self.config.time.save_every == 0:
                self._write_fields()
            if round(self.time, 12) in snapshots and self.write_output:
                self._write_snapshot(f"t{self.time:g}")
        self._close_output()
        return PDEResult(
            output_directory=self.output_directory,
            final_state=self.capture(),
            times=np.asarray(self._times),
            metrics={name: np.asarray(values) for name, values in self._history.items()},
        )


def run_pde(
    config: RunConfig,
    *,
    output_directory: str | Path | None = None,
    final_time: float | None = None,
    initial_state: SimulationState | None = None,
    snapshot_times: Iterable[float] = (),
    write_output: bool = True,
    comm: Any = None,
) -> PDEResult:
    runner = PDERunner(
        config,
        output_directory=output_directory,
        initial_state=initial_state,
        write_output=write_output,
        comm=comm,
    )
    return runner.run(final_time=final_time, snapshot_times=snapshot_times)
