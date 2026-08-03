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
    parser.add_argument(
        "--source-crs",
        help="CRS to record when the source lacks embedded CRS metadata (for example EPSG:3435)",
    )
    parser.add_argument(
        "--provenance",
        type=Path,
        help="JSON object containing dataset-level provenance to preserve in the output metadata",
    )
    args = parser.parse_args()

    source = Path(args.source)
    prefix = Path(args.output_prefix)
    frame = gpd.read_file(source)
    embedded_crs = None if frame.crs is None else str(frame.crs)
    if embedded_crs is not None and args.source_crs is not None:
        requested_crs = str(args.source_crs)
        if frame.crs != requested_crs:
            raise ValueError(
                f"source embeds CRS {embedded_crs}, which conflicts with "
                f"--source-crs {requested_crs}"
            )
    source_crs = embedded_crs or args.source_crs

    provenance = {}
    if args.provenance is not None:
        provenance = json.loads(args.provenance.read_text(encoding="utf-8"))
        if not isinstance(provenance, dict):
            raise TypeError("provenance JSON must contain an object")

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
        **provenance,
        "source": str(source),
        "source_files_sha256": {path.name: checksum(path) for path in sidecars},
        "source_crs": source_crs,
        "source_crs_embedded": embedded_crs is not None,
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
