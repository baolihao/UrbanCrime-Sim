"""Generic monolithic Newton solver for a preassembled UFL residual."""

from __future__ import annotations

from typing import Any, Mapping


class ImplicitSolver:
    def __init__(
        self,
        *,
        mesh: Any,
        residual: Any,
        solution: Any,
        relative_tolerance: float,
        absolute_tolerance: float,
        max_iterations: int,
        convergence_criterion: str,
        petsc_options: Mapping[str, str],
    ) -> None:
        from dolfinx.fem.petsc import NonlinearProblem
        from dolfinx.nls.petsc import NewtonSolver
        from petsc4py import PETSc

        if relative_tolerance <= 0 or absolute_tolerance <= 0 or max_iterations < 1:
            raise ValueError("invalid Newton solver tolerances or iteration limit")
        if convergence_criterion not in {"incremental", "residual"}:
            raise ValueError("invalid Newton convergence criterion")
        problem = NonlinearProblem(residual, solution)
        self._solver = NewtonSolver(mesh.comm, problem)
        self._solver.convergence_criterion = convergence_criterion
        self._solver.rtol = relative_tolerance
        self._solver.atol = absolute_tolerance
        self._solver.max_it = max_iterations
        self._solver.error_on_nonconvergence = True

        ksp = self._solver.krylov_solver
        prefix = ksp.getOptionsPrefix()
        options = PETSc.Options()
        for key, value in petsc_options.items():
            options[f"{prefix}{key}"] = value
        ksp.setFromOptions()

    def solve(self, solution: Any) -> tuple[int, bool]:
        result = self._solver.solve(solution)
        solution.x.scatter_forward()
        if isinstance(result, tuple):
            return int(result[0]), bool(result[1])
        return int(result), True
