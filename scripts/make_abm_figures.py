#!/usr/bin/env python3
"""Plot field panels or mean time series from a completed Python ABM run."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def make_panels(run: Path, requested_times: list[float], fields: list[str], output: Path) -> None:
    with np.load(run / "trajectory.npz") as data:
        available = np.asarray(data["time"], dtype=float)
        missing = set(fields).difference(data.files)
        if missing:
            raise KeyError(f"trajectory is missing fields: {sorted(missing)}")
        indices: list[int] = []
        for value in requested_times:
            matches = np.flatnonzero(np.isclose(available, value, rtol=0, atol=1e-9))
            if matches.size == 0:
                raise ValueError(
                    f"t={value:g} was not saved; available times are {available.tolist()}"
                )
            indices.append(int(matches[0]))
        arrays = {name: np.asarray(data[name])[indices] for name in fields}

    single_time = len(indices) == 1
    rows = 1 if single_time else len(fields)
    columns = len(fields) if single_time else len(indices)
    figure, axes = plt.subplots(
        rows,
        columns,
        figsize=(3 * columns, 2.7 * rows),
        squeeze=False,
        constrained_layout=True,
    )
    for row, field in enumerate(fields):
        minimum = float(np.min(arrays[field]))
        maximum = float(np.max(arrays[field]))
        if np.isclose(minimum, maximum):
            maximum = minimum + 1e-12
        for column, (time, values) in enumerate(zip(requested_times, arrays[field], strict=True)):
            axis = axes[0, row] if single_time else axes[row, column]
            artist = axis.imshow(
                values, origin="lower", cmap="turbo", vmin=minimum, vmax=maximum
            )
            axis.set_xticks([])
            axis.set_yticks([])
            if single_time:
                axis.set_title(field)
            elif row == 0:
                axis.set_title(f"t = {time:g}")
            if not single_time and column == 0:
                axis.set_ylabel(field, rotation=0, labelpad=22)
        colorbar_axes = [axes[0, row]] if single_time else axes[row, :].tolist()
        figure.colorbar(artist, ax=colorbar_axes, shrink=0.78, pad=0.01)
    if single_time:
        figure.suptitle(f"t = {requested_times[0]:g}")
    figure.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(figure)


METRIC_LABELS = {
    "mean_A": "A",
    "mean_rho": "rho",
    "mean_H": "H",
    "mean_pi": "pi",
    "mean_expected_S": "expected S",
    "mean_realized_S": "realized S",
}

METRIC_MATH_LABELS = {
    "mean_A": "A",
    "mean_rho": r"\rho",
    "mean_H": "H",
    "mean_pi": r"\pi",
    "mean_expected_S": "S",
    "mean_realized_S": r"S_{\mathrm{realized}}",
}


def make_timeseries(run: Path, metrics: list[str], output: Path, *, overlay: bool) -> None:
    with np.load(run / "metrics.npz") as data:
        time = data["time"]
        missing = set(metrics).difference(data.files)
        if missing:
            raise KeyError(f"metrics file is missing fields: {sorted(missing)}")
        if overlay:
            figure, axis = plt.subplots(figsize=(9, 4.8), constrained_layout=True)
            styles = {
                "mean_A": {"color": "tab:blue", "linestyle": "-"},
                "mean_rho": {"color": "tab:green", "linestyle": "-"},
                "mean_H": {"color": "tab:red", "linestyle": "--"},
                "mean_pi": {"color": "0.35", "linestyle": "--"},
            }
            for key in metrics:
                label = rf"$\langle {METRIC_MATH_LABELS.get(key, key)} \rangle$"
                axis.plot(time, data[key], label=label, **styles.get(key, {}))
            axis.set(xlabel="time", ylabel="spatial average")
            axis.grid(alpha=0.25)
            axis.legend(ncol=2)
            figure.savefig(output, dpi=220)
            plt.close(figure)
            return
        columns = min(2, len(metrics))
        rows = (len(metrics) + columns - 1) // columns
        figure, axes = plt.subplots(
            rows,
            columns,
            figsize=(4.5 * columns, 3.5 * rows),
            squeeze=False,
            constrained_layout=True,
        )
        for axis, key in zip(axes.ravel(), metrics, strict=False):
            axis.plot(time, data[key])
            axis.set(xlabel="time", ylabel=METRIC_LABELS.get(key, key))
            axis.grid(alpha=0.25)
        for axis in axes.ravel()[len(metrics) :]:
            axis.set_visible(False)
    figure.savefig(output, dpi=220)
    plt.close(figure)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_directory", type=Path)
    parser.add_argument("--times", nargs="+", type=float)
    parser.add_argument("--fields", nargs="+", default=["A", "rho"])
    parser.add_argument(
        "--metrics",
        nargs="+",
        default=["mean_A", "mean_rho", "mean_H", "mean_pi"],
        help="metrics.npz arrays to plot when --times is omitted",
    )
    parser.add_argument(
        "--overlay",
        action="store_true",
        help="draw selected metrics on one axis, as in extension-paper Figure 6",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.times:
        make_panels(
            args.run_directory,
            args.times,
            args.fields,
            args.output or args.run_directory / "field_panels.png",
        )
    else:
        make_timeseries(
            args.run_directory,
            args.metrics,
            args.output or args.run_directory / "mean_timeseries.png",
            overlay=args.overlay,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
