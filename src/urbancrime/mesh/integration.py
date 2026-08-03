"""Mass-lumped integration for nodal policy fields on arbitrary meshes."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray


@dataclass(frozen=True)
class LumpedMeasure:
    weights: NDArray[np.float64]
    owned_size: int
    comm: Any = None

    @classmethod
    def from_function_space(cls, function_space: Any) -> "LumpedMeasure":
        import ufl
        from dolfinx import fem
        from dolfinx.fem import petsc
        from petsc4py import PETSc

        test = ufl.TestFunction(function_space)
        vector = petsc.assemble_vector(fem.form(test * ufl.dx))
        vector.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
        weights = np.asarray(vector.array, dtype=float).copy()
        owned = function_space.dofmap.index_map.size_local
        return cls(weights=weights, owned_size=owned, comm=function_space.mesh.comm)

    def local_integral(self, values: ArrayLike) -> float:
        array = np.asarray(values, dtype=float)
        if array.shape != self.weights.shape:
            raise ValueError("values and lumped weights must have the same local shape")
        return float(np.dot(self.weights[: self.owned_size], array[: self.owned_size]))

    def integral(self, values: ArrayLike) -> float:
        local = self.local_integral(values)
        if self.comm is None:
            return local
        from mpi4py import MPI

        return float(self.comm.allreduce(local, op=MPI.SUM))

    @property
    def domain_area(self) -> float:
        return self.integral(np.ones_like(self.weights))
