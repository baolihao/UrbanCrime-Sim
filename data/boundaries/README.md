# Boundary data

`raw/` contains source boundary files.  `processed/` contains derived arrays and
metadata created by `build_boundaries.py`.

Each raw dataset must have metadata recording its source, terms, CRS,
acquisition date, processing steps, and checksums.  The repository software
license does not apply to third-party geographic inputs; see [NOTICE.md](NOTICE.md).

Only the Chicago boundary used by the published and delayed-policing examples is
tracked.  Its `.prj` sidecar was lost before it entered this repository.  EPSG:3435
is restored as source metadata from the official City GIS layer and the matching
coordinate range; the raw bytes are left unchanged.  Derived coordinates are
normalized to a long side of 5, and paper-scale runs apply an explicit
`domain.scale_factor` in YAML.

Regenerate the processed array and metadata with:

```bash
python data/boundaries/build_boundaries.py \
  data/boundaries/raw/chicago/chicago_simplified.shp \
  data/boundaries/processed/chicago \
  --target-long-side 5 \
  --source-crs EPSG:3435 \
  --provenance data/boundaries/raw/chicago/chicago_simplified.provenance.json
```

Dataset-level provenance is stored beside the raw geometry; the builder combines
it with freshly calculated processing fields and file checksums.
