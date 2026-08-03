#!/usr/bin/env python3
"""Extract a normalized exterior boundary and provenance metadata."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import geopandas as gpd
import numpy as np


def checksum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source")
    parser.add_argument("output_prefix")
    parser.add_argument("--target-long-side", type=float)
    args = parser.parse_args()

    source = Path(args.source)
    prefix = Path(args.output_prefix)
    frame = gpd.read_file(source)
    geometry = frame.geometry.union_all()
    if geometry.geom_type == "MultiPolygon":
        geometry = max(geometry.geoms, key=lambda item: item.area)
    if geometry.geom_type != "Polygon":
        raise TypeError(f"expected polygonal boundary, got {geometry.geom_type}")

    boundary = np.asarray(geometry.exterior.coords, dtype=float)
    boundary -= boundary.min(axis=0)
    original_long_side = float(np.ptp(boundary, axis=0).max())
    scale = 1.0
    if args.target_long_side is not None:
        if args.target_long_side <= 0 or original_long_side <= 0:
            raise ValueError("target and original long-side lengths must be positive")
        scale = args.target_long_side / original_long_side
        boundary *= scale

    prefix.parent.mkdir(parents=True, exist_ok=True)
    np.save(prefix.with_suffix(".npy"), boundary)
    sidecars = sorted(source.parent.glob(f"{source.stem}.*"))
    metadata = {
        "source": str(source),
        "source_files_sha256": {path.name: checksum(path) for path in sidecars},
        "source_crs": None if frame.crs is None else str(frame.crs),
        "source_feature_count": len(frame),
        "selected_geometry": "largest polygon exterior",
        "hole_count": len(geometry.interiors),
        "original_long_side": original_long_side,
        "target_long_side": args.target_long_side,
        "scale": scale,
    }
    with prefix.with_suffix(".json").open("w", encoding="utf-8") as stream:
        json.dump(metadata, stream, indent=2)
        stream.write("\n")


if __name__ == "__main__":
    main()
