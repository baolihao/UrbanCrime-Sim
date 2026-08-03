"""Gmsh polygon meshes with optional holes and boundary-band sizing."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np

from .gmsh_size_fields import add_boundary_threshold_field


def create_polygon(
    exterior: np.ndarray,
    *,
    holes: Sequence[np.ndarray] = (),
    characteristic_length: float,
    boundary_length: float | None = None,
    interior_length: float | None = None,
    transition_width: float | None = None,
    comm: Any = None,
) -> tuple[Any, Any, Any]:
    import gmsh
    from dolfinx.io import gmshio
    from mpi4py import MPI

    if characteristic_length <= 0:
        raise ValueError("characteristic_length must be positive")
    mpi_comm = MPI.COMM_WORLD if comm is None else comm
    rank = 0

    gmsh.initialize()
    try:
        if mpi_comm.rank == rank:
            gmsh.model.add("urbancrime_polygon")

            def ring_loop(points: np.ndarray) -> tuple[int, list[int]]:
                coords = np.asarray(points, dtype=float)
                if coords.ndim != 2 or coords.shape[1] != 2 or len(coords) < 3:
                    raise ValueError("each polygon ring must have shape (n, 2), n >= 3")
                if np.allclose(coords[0], coords[-1]):
                    coords = coords[:-1]
                tags = [
                    gmsh.model.occ.addPoint(float(x), float(y), 0.0, characteristic_length)
                    for x, y in coords
                ]
                curves = [
                    gmsh.model.occ.addLine(tags[i], tags[(i + 1) % len(tags)])
                    for i in range(len(tags))
                ]
                return gmsh.model.occ.addCurveLoop(curves), curves

            outer_loop, outer_curves = ring_loop(exterior)
            hole_loops: list[int] = []
            all_curves = list(outer_curves)
            for hole in holes:
                loop, curves = ring_loop(hole)
                hole_loops.append(loop)
                all_curves.extend(curves)
            surface = gmsh.model.occ.addPlaneSurface([outer_loop, *hole_loops])
            gmsh.model.occ.synchronize()
            gmsh.model.addPhysicalGroup(2, [surface], tag=1)
            gmsh.model.setPhysicalName(2, 1, "domain")
            gmsh.model.addPhysicalGroup(1, all_curves, tag=2)
            gmsh.model.setPhysicalName(1, 2, "boundary")

            refinement = (boundary_length, interior_length, transition_width)
            if any(value is not None for value in refinement):
                if not all(value is not None for value in refinement):
                    raise ValueError("boundary refinement requires all three sizing parameters")
                add_boundary_threshold_field(
                    gmsh,
                    all_curves,
                    boundary_length=float(boundary_length),
                    interior_length=float(interior_length),
                    transition_width=float(transition_width),
                )
            gmsh.option.setNumber("Mesh.MeshSizeMin", boundary_length or characteristic_length)
            gmsh.option.setNumber("Mesh.MeshSizeMax", interior_length or characteristic_length)
            gmsh.model.mesh.generate(2)

        converted = gmshio.model_to_mesh(gmsh.model, mpi_comm, rank, gdim=2)
        if hasattr(converted, "mesh"):
            return converted.mesh, converted.cell_tags, converted.facet_tags
        return converted
    finally:
        gmsh.finalize()
