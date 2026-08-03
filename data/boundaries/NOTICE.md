# Boundary-data notice

## Chicago city boundary

The Chicago files are derived from the City of Chicago Data Portal dataset
[Boundaries - City](https://data.cityofchicago.org/Facilities-Geographic-Boundaries/Boundaries-City/qqq8-j68g)
(dataset ID `qqq8-j68g`), provided by the City of Chicago.

The dataset page does not grant the repository's BSD-3-Clause software license;
it directs users to the City's terms of use.  Users and redistributors of these
files remain responsible for following those terms, including the attribution
and disclaimer requirements for modified or derivative applications.  The City
does not warrant the accuracy, timeliness, or completeness of this modified
geometry, and use of it is at the user's own risk.  See the
[City data terms and disclaimer](https://www.chicago.gov/city/en/narr/foia/data_disclaimer.html).

### Provenance

- Provider: City of Chicago.
- Source dataset: `Boundaries - City` (`qqq8-j68g`).
- Source CRS: EPSG:3435, NAD83 / Illinois East (ftUS).
- CRS evidence: the coordinate range of the historical local file matches the
  City boundary layer in the City's
  [official GIS service](https://gisapps.cityofchicago.org/arcgis/rest/services/CachedMaps/CartoCache/MapServer),
  which identifies spatial reference 3435 and feet as its units.
- Acquisition date: unknown; this is a historical download.
- Local raw input: `raw/chicago/chicago_simplified.shp` and its sidecars.  The
  file is already simplified to one closed exterior ring with 33 coordinates
  and has no `.prj` file.  The method and date of that earlier simplification
  cannot be reconstructed from the files currently available.
- Repository processing: select the largest polygon exterior, discard holes and
  attributes, translate the coordinate minima to zero, and scale the longest
  side to 5 model units.  The generated array and machine-readable metadata are
  under `processed/`.

The processed coordinates are dimensionless simulation coordinates and must not
be treated as a survey-grade or authoritative city boundary.
