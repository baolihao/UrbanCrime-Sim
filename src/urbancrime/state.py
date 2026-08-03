"""Small immutable snapshots used to branch policy-comparison runs."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np
from numpy.typing import NDArray


@dataclass(frozen=True)
class SimulationState:
    time: float
    fields: Mapping[str, NDArray[np.float64]]

    @classmethod
    def capture(cls, time: float, **fields: NDArray[np.float64]) -> "SimulationState":
        frozen: dict[str, NDArray[np.float64]] = {}
        for name, values in fields.items():
            copied = np.asarray(values, dtype=float).copy()
            copied.setflags(write=False)
            frozen[name] = copied
        return cls(time=float(time), fields=frozen)

    def mutable_copy(self) -> dict[str, NDArray[np.float64]]:
        return {name: values.copy() for name, values in self.fields.items()}
