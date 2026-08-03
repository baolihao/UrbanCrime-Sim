#!/usr/bin/env python3
"""Run homogeneous stability sweeps from a YAML study configuration."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, replace
from pathlib import Path
from typing import Any, Mapping

import matplotlib.pyplot as plt
import numpy as np
import yaml

from urbancrime.analysis.stability import (
    StabilityParameters,
    critical_delay,
    parameter_sweep,
    transversality,
)
from urbancrime.provenance import environment_metadata, git_revision, write_checksums


def values_from_spec(spec: Mapping[str, Any]) -> np.ndarray:
    scale = spec["scale"]
    if scale == "linear":
        return np.linspace(
            float(spec["start"]),
            float(spec["stop"]),
            int(spec["count"]),
            endpoint=bool(spec.get("endpoint", True)),
        )
    if scale == "log":
        return np.geomspace(float(spec["start"]), float(spec["stop"]), int(spec["count"]))
    if scale == "concatenate":
        return np.concatenate([values_from_spec(part) for part in spec["parts"]])
    raise ValueError(f"unknown sweep scale: {scale!r}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("config", help="stability study YAML file")
    args = parser.parse_args()

    config_path = Path(args.config)
    with config_path.open("r", encoding="utf-8") as stream:
        raw = yaml.safe_load(stream)

    parameters = StabilityParameters(**raw["base_parameters"])
    search = {
        "max_mode": int(raw["spectral"]["max_mode"]),
        "tau_min": float(raw["root_search"]["tau_min"]),
        "tau_max": float(raw["root_search"]["tau_max"]),
        "bracket_samples": int(raw["root_search"]["bracket_samples"]),
        "root_tolerance": float(raw["root_search"]["root_tolerance"]),
    }
    base_point = critical_delay(parameters, **search)
    if base_point is None:
        raise RuntimeError("no Hopf point found for the base parameters")
    derivative = transversality(
        base_point,
        parameters,
        relative_step=float(raw["transversality_relative_step"]),
    )

    output_spec = raw["output"]
    output_dir = Path(output_spec["root"]) / output_spec["name"]
    output_dir.mkdir(parents=True, exist_ok=True)

    sweep_results: dict[str, dict[str, np.ndarray]] = {}
    for sweep in raw["sweeps"]:
        name = str(sweep["parameter"])
        values = values_from_spec(sweep["values"])
        fixed = {key: float(value) for key, value in sweep.get("fixed", {}).items()}
        local_base = replace(parameters, **fixed)
        delays = parameter_sweep(local_base, name, values, **search)
        sweep_results[name] = {"values": values, "critical_delays": delays}

    formats = set(output_spec["formats"])
    if "npz" in formats:
        arrays: dict[str, np.ndarray] = {}
        for name, result in sweep_results.items():
            arrays[f"{name}_values"] = result["values"]
            arrays[f"{name}_critical_delays"] = result["critical_delays"]
        np.savez(output_dir / "phase_diagrams.npz", **arrays)

    metadata = {
        "schema_version": raw["schema_version"],
        "name": raw["name"],
        "base_parameters": asdict(parameters),
        "critical_point": asdict(base_point),
        "transversality": derivative,
        "search": search,
        "git_revision": git_revision(Path.cwd()),
        "environment": environment_metadata(("numpy", "scipy", "matplotlib", "PyYAML")),
    }
    with (output_dir / "config.resolved.yaml").open("w", encoding="utf-8") as stream:
        yaml.safe_dump(raw, stream, sort_keys=False)
    if "json" in formats:
        with (output_dir / "metadata.json").open("w", encoding="utf-8") as stream:
            json.dump(metadata, stream, indent=2)
            stream.write("\n")

    if "png" in formats:
        figure, axes = plt.subplots(1, len(sweep_results), figsize=(5.4 * len(sweep_results), 5))
        axes = np.atleast_1d(axes)
        labels = {"eta": r"$\eta$", "q_bar": r"$\bar q$", "pi_bar": r"$\bar\pi$"}
        for axis, (name, result) in zip(axes, sweep_results.items()):
            values = result["values"]
            delays = result["critical_delays"]
            finite = np.isfinite(delays)
            axis.semilogy(values[finite], delays[finite], color="black", linewidth=2)
            axis.fill_between(values[finite], search["tau_min"], delays[finite], alpha=0.3)
            axis.fill_between(values[finite], delays[finite], search["tau_max"], alpha=0.25)
            axis.set_xlabel(labels.get(name, name))
            axis.set_ylabel(r"$\tau_c^*$")
            axis.grid(alpha=0.25, which="both")
        figure.tight_layout()
        figure.savefig(output_dir / "phase_diagrams.png", dpi=200, bbox_inches="tight")

    write_checksums(output_dir)

    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
