"""Shared interfaces for police strategies."""

from __future__ import annotations

from typing import Any, Protocol, runtime_checkable


@runtime_checkable
class PoliceStrategy(Protocol):
    name: str
    state_names: tuple[str, ...]

    def pi_expression(self, state: Any, time: float) -> Any:
        """Return police density pi, never the attenuation exp(-pi)."""


def require_positive(value: float, name: str) -> float:
    if value <= 0:
        raise ValueError(f"{name} must be positive")
    return value
