"""Reproducible multi-run studies built on the common PDE runtime."""

from __future__ import annotations

from dataclasses import replace
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import yaml

from .config import (
    DelayedPoliceConfig,
    NoPoliceConfig,
    OutputConfig,
    RunConfig,
    load_run_config,
    police_config_from_mapping,
)
from .runtime import PDEResult, run_pde


def _mapping(value: Any, name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise TypeError(f"{name} must be a mapping")
    return value


def compare_policies(path: str | Path) -> dict[str, PDEResult]:
    """Warm up once, then branch prescribed, optimal, and delayed strategies."""

    from mpi4py import MPI

    source = Path(path)
    with source.open(encoding="utf-8") as stream:
        study = _mapping(yaml.safe_load(stream), str(source))
    base_path = (source.parent / str(study["base_simulation"])).resolve()
    base = load_run_config(base_path)
    warmup = _mapping(study["warmup"], "warmup")
    activation_time = float(warmup["activation_time"])
    if str(warmup["strategy"]) != "none":
        raise ValueError("policy comparison warmup must use strategy: none")
    output = _mapping(study["output"], "output")
    study_root = Path(output["root"]) / str(output["name"])
    warm_config = replace(
        base,
        name=f"{output['name']}-warmup",
        policing=NoPoliceConfig(),
        solver=replace(base.solver, block_order=("A", "rho")),
        output=OutputConfig(root=study_root, fields=("A", "rho"), format="xdmf"),
    )
    warm_result = run_pde(
        warm_config,
        output_directory=study_root / "warmup",
        final_time=activation_time,
        snapshot_times=(activation_time,),
    )
    branch_state = warm_result.final_state
    results: dict[str, PDEResult] = {}
    for raw_strategy in study["strategies"]:
        strategy = police_config_from_mapping(_mapping(raw_strategy, "strategy"))
        strategy_name = strategy.strategy
        if hasattr(strategy, "activation_time") and float(strategy.activation_time) != activation_time:
            raise ValueError("all policy branches must activate at the warmup time")
        fields = ("A", "rho", "pi", "H") if isinstance(strategy, DelayedPoliceConfig) else ("A", "rho", "pi")
        block_order = ("A", "rho", "H", "pi") if isinstance(strategy, DelayedPoliceConfig) else ("A", "rho")
        branch_config = replace(
            base,
            name=f"{output['name']}-{strategy_name}",
            policing=strategy,
            solver=replace(base.solver, block_order=block_order),
            output=OutputConfig(root=study_root, fields=fields, format="xdmf"),
        )
        results[strategy_name] = run_pde(
            branch_config,
            output_directory=study_root / strategy_name,
            initial_state=branch_state,
            snapshot_times=tuple(float(value) for value in study["snapshot_times"]),
        )
    if MPI.COMM_WORLD.rank == 0:
        study_root.mkdir(parents=True, exist_ok=True)
        combined: dict[str, np.ndarray] = {}
        for strategy_name, result in results.items():
            combined[f"{strategy_name}_time"] = result.times
            for metric_name, values in result.metrics.items():
                combined[f"{strategy_name}_{metric_name}"] = values
        np.savez_compressed(study_root / "comparison.npz", **combined)
        (study_root / "study.resolved.json").write_text(
            json.dumps(study, indent=2), encoding="utf-8"
        )
        _plot_policy_comparison(study_root, results)
    MPI.COMM_WORLD.barrier()
    return results


def _plot_policy_comparison(root: Path, results: Mapping[str, PDEResult]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    labels = {"prescribed": "Prescribed", "optimal": "Optimal", "delayed": "Delayed"}
    for name, result in results.items():
        axes[0].plot(result.times, result.metrics["crime_average"], label=labels.get(name, name))
        axes[1].plot(result.times, result.metrics["crime_maximum"], label=labels.get(name, name))
    axes[0].set(xlabel="time", ylabel="spatial average", title="Average crime intensity")
    axes[1].set(xlabel="time", ylabel="maximum", title="Maximum crime intensity")
    for axis in axes:
        axis.grid(alpha=0.25)
        axis.legend()
    figure.savefig(root / "policy_comparison.png", dpi=200)
    plt.close(figure)
