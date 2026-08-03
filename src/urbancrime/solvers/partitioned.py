"""Configurable fixed-point orchestration for partitioned field solves."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Sequence

import numpy as np
from numpy.typing import ArrayLike


def relative_l2_change(new: ArrayLike, old: ArrayLike, absolute_floor: float) -> float:
    if absolute_floor <= 0:
        raise ValueError("absolute_floor must be positive")
    new_values = np.asarray(new, dtype=float)
    old_values = np.asarray(old, dtype=float)
    numerator = float(np.linalg.norm(new_values - old_values))
    denominator = max(float(np.linalg.norm(old_values)), absolute_floor)
    return numerator / denominator


@dataclass(frozen=True)
class PartitionBlock:
    name: str
    solve: Callable[[], float]


@dataclass(frozen=True)
class PartitionedStepReport:
    iterations: int
    residuals: dict[str, float]


class PartitionedSolver:
    """Run ordered field blocks until every reported relative change converges.

    A two-field strategy supplies ``A`` and ``rho`` blocks.  Delayed policing
    supplies ``A``, ``rho``, ``H``, and ``pi`` blocks.  The orchestration is
    shared without pretending that inactive fields exist.
    """

    def __init__(
        self,
        blocks: Sequence[PartitionBlock],
        *,
        relative_tolerance: float,
        max_iterations: int,
    ) -> None:
        if not blocks:
            raise ValueError("at least one partition block is required")
        if relative_tolerance <= 0 or max_iterations < 1:
            raise ValueError("invalid partitioned solver controls")
        names = [block.name for block in blocks]
        if len(names) != len(set(names)):
            raise ValueError("partition block names must be unique")
        self.blocks = tuple(blocks)
        self.relative_tolerance = relative_tolerance
        self.max_iterations = max_iterations

    def solve_step(self) -> PartitionedStepReport:
        residuals: dict[str, float] = {}
        for iteration in range(1, self.max_iterations + 1):
            residuals = {block.name: float(block.solve()) for block in self.blocks}
            if not all(np.isfinite(value) for value in residuals.values()):
                raise FloatingPointError(f"non-finite partition residuals: {residuals}")
            if max(residuals.values()) <= self.relative_tolerance:
                return PartitionedStepReport(iterations=iteration, residuals=residuals)
        raise RuntimeError(
            f"partitioned solve did not converge in {self.max_iterations} iterations; "
            f"last residuals: {residuals}"
        )
