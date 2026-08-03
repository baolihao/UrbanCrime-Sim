# Boundary data

`raw/` contains source boundary files.  `processed/` contains derived arrays and
metadata created by `build_boundaries.py`.

Before publication, each raw dataset must have metadata recording its source,
license, CRS, acquisition date, processing steps, and checksum.  The current
Chicago shapefile has no `.prj` sidecar, so its CRS must be established before it
is treated as a reproducible geographic input.
Only the Chicago boundary used by the published and delayed-policing examples is
tracked.  Derived coordinates are normalized to a long side of 5; paper-scale
runs apply an explicit `domain.scale_factor` in YAML.
