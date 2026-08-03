"""Structured rectangular meshes."""

from __future__ import annotations

from typing import Any, Sequence

import numpy as np


def create_rectangle(
    *,
    lower_left: Sequence[float],
    upper_right: Sequence[float],
    nx: int,
    ny: int,
    cell_type: str,
    comm: Any = None,
) -> Any:
    from dolfinx import mesh
    from mpi4py import MPI

    if nx < 1 or ny < 1:
        raise ValueError("nx and ny must be positive")
    if cell_type not in {"triangle", "quadrilateral"}:
        raise ValueError("cell_type must be triangle or quadrilateral")
    mpi_comm = MPI.COMM_WORLD if comm is None else comm
    kind = mesh.CellType.triangle if cell_type == "triangle" else mesh.CellType.quadrilateral
    return mesh.create_rectangle(
        mpi_comm,
        [np.asarray(lower_left, dtype=float), np.asarray(upper_right, dtype=float)],
        [nx, ny],
        cell_type=kind,
    )
