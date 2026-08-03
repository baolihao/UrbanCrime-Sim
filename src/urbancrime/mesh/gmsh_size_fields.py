"""Reusable Gmsh size fields."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any


def add_boundary_threshold_field(
    gmsh: Any,
    boundary_curves: Sequence[int],
    *,
    boundary_length: float,
    interior_length: float,
    transition_width: float,
) -> int:
    if boundary_length <= 0 or interior_length <= 0 or transition_width <= 0:
        raise ValueError("mesh lengths and transition width must be positive")
    if boundary_length > interior_length:
        raise ValueError("boundary_length must not exceed interior_length")

    distance = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(distance, "CurvesList", list(boundary_curves))
    threshold = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(threshold, "InField", distance)
    gmsh.model.mesh.field.setNumber(threshold, "SizeMin", boundary_length)
    gmsh.model.mesh.field.setNumber(threshold, "SizeMax", interior_length)
    gmsh.model.mesh.field.setNumber(threshold, "DistMin", 0.0)
    gmsh.model.mesh.field.setNumber(threshold, "DistMax", transition_width)
    gmsh.model.mesh.field.setAsBackgroundMesh(threshold)
    return threshold
