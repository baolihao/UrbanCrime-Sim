"""No-police closure."""

from __future__ import annotations

from typing import Any


class NoPolice:
    name = "none"
    state_names: tuple[str, ...] = ()

    def pi_expression(self, state: Any = None, time: float = 0.0) -> float:
        del state, time
        return 0.0
