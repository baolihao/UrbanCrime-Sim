"""Test-suite capability checks for optional native dependencies."""

from __future__ import annotations

from importlib.util import find_spec

import pytest


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    if find_spec("dolfinx") is not None and find_spec("ufl") is not None:
        return
    reason = pytest.mark.skip(reason="requires the conda FEniCSx environment")
    for item in items:
        if "fenics" in item.keywords:
            item.add_marker(reason)
