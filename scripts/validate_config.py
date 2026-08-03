#!/usr/bin/env python3
"""Validate and print a fully resolved simulation configuration."""

from __future__ import annotations

import argparse
import json

from urbancrime.config import load_run_config


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("config", help="simulation YAML file")
    args = parser.parse_args()
    config = load_run_config(args.config)
    print(json.dumps(config.resolved_mapping(), indent=2))


if __name__ == "__main__":
    main()
