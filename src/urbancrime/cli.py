"""Console entry points kept deliberately thin."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

from .config import load_run_config


def pde_main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run the heterogeneous-flux PDE model")
    parser.add_argument("config", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--final-time", type=float)
    parser.add_argument("--snapshot", type=float, action="append", default=[])
    args = parser.parse_args(argv)
    from .runtime import run_pde

    result = run_pde(
        load_run_config(args.config),
        output_directory=args.output,
        final_time=args.final_time,
        snapshot_times=args.snapshot,
    )
    if result.final_state.fields and result.times.size:
        print(f"completed t={result.times[-1]:g}; output={result.output_directory}")
    return 0


def policy_main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Compare policing strategies from one warmup")
    parser.add_argument("config", type=Path)
    args = parser.parse_args(argv)
    from .studies import compare_policies

    compare_policies(args.config)
    return 0


def abm_main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run the Python burglary ABM")
    parser.add_argument("config", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    from .abm import load_abm_config, run_abm

    result = run_abm(load_abm_config(args.config), output_directory=args.output)
    print(f"completed t={result.times[-1]:g}; output={result.output_directory}")
    return 0
