#!/usr/bin/env python3
"""Create compact time-series or field-panel figures from a completed run."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np


def _triangles(cells: np.ndarray, points: np.ndarray) -> np.ndarray:
    connectivity = np.asarray(cells, dtype=np.int64)
    if connectivity.ndim != 2:
        raise ValueError("snapshot cells must be a two-dimensional array")
    if connectivity.shape[1] == 3:
        return connectivity
    if connectivity.shape[1] == 4:
        coordinates = points[connectivity]
        centers = coordinates.mean(axis=1, keepdims=True)
        angles = np.arctan2(
            coordinates[:, :, 1] - centers[:, :, 1],
            coordinates[:, :, 0] - centers[:, :, 0],
        )
        ordered = np.take_along_axis(connectivity, np.argsort(angles, axis=1), axis=1)
        return np.vstack([ordered[:, [0, 1, 2]], ordered[:, [0, 2, 3]]])
    raise ValueError("only triangular and quadrilateral P1 snapshots are supported")


def make_timeseries(run_directory: Path, output: Path) -> None:
    with np.load(run_directory / "timeseries.npz") as data:
        time = data["time"]
        figure, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
        axes[0].plot(time, data["crime_average"])
        axes[1].plot(time, data["crime_maximum"])
    axes[0].set(xlabel="time", ylabel="spatial average", title="Average crime intensity")
    axes[1].set(xlabel="time", ylabel="maximum", title="Maximum crime intensity")
    for axis in axes:
        axis.grid(alpha=0.25)
    figure.savefig(output, dpi=200)
    plt.close(figure)


def make_field_panels(
    run_directory: Path,
    labels: list[str],
    fields: list[str],
    output: Path,
) -> None:
    snapshots: list[dict[str, np.ndarray]] = []
    for label in labels:
        source = run_directory / f"snapshot_{label}.npz"
        with np.load(source) as data:
            missing = {"time", "points", "cells", *fields}.difference(data.files)
            if missing:
                raise KeyError(f"{source} is missing arrays: {sorted(missing)}")
            snapshots.append({name: np.asarray(data[name]) for name in data.files})

    figure, axes = plt.subplots(
        len(fields),
        len(snapshots),
        figsize=(3.1 * len(snapshots), 2.8 * len(fields)),
        squeeze=False,
        constrained_layout=True,
    )
    for row, field in enumerate(fields):
        minimum = min(float(snapshot[field].min()) for snapshot in snapshots)
        maximum = max(float(snapshot[field].max()) for snapshot in snapshots)
        if np.isclose(minimum, maximum):
            maximum = minimum + 1.0e-12
        levels = np.linspace(minimum, maximum, 80)
        for column, snapshot in enumerate(snapshots):
            points = snapshot["points"]
            triangulation = mtri.Triangulation(
                points[:, 0], points[:, 1], _triangles(snapshot["cells"], points)
            )
            artist = axes[row, column].tricontourf(
                triangulation,
                snapshot[field].reshape(-1),
                levels=levels,
                cmap="turbo",
                extend="both",
            )
            axes[row, column].set_aspect("equal")
            axes[row, column].set_xticks([])
            axes[row, column].set_yticks([])
            if row == 0:
                axes[row, column].set_title(f"t = {float(snapshot['time']):g}")
            if column == 0:
                axes[row, column].set_ylabel(field, rotation=0, labelpad=18)
        figure.colorbar(artist, ax=axes[row, :].tolist(), shrink=0.78, pad=0.01)
    figure.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(figure)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_directory", type=Path)
    parser.add_argument("--snapshots", nargs="+", help="snapshot labels such as t0 t2 t200")
    parser.add_argument("--fields", nargs="+", default=["A", "rho"])
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    if args.snapshots:
        output = args.output or args.run_directory / "field_panels.png"
        make_field_panels(args.run_directory, args.snapshots, args.fields, output)
    else:
        output = args.output or args.run_directory / "crime_timeseries.png"
        make_timeseries(args.run_directory, output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
