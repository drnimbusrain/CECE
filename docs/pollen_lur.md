# Pollen Sense station-to-grid land-use regression

`tools/prepare_pollen_lur.py` is a provider-neutral workflow for fitting a simple
station-to-grid pollen model from the documented CSV contract below and a
preprocessed, co-registered predictor NetCDF. It writes a nonnegative
`base_pollen_potential` field in **grains m-2 s-1**. This is an independently
calibrated base emission-potential product, not a retrieval or inversion of
atmospheric transport, deposition, or source strength.

## Reproducible synthetic example

The fixture generator creates freely generated data only; it uses no network,
proprietary service, credentials, or satellite data. The predictor NetCDF is
marked with synthetic provenance, and the paired station CSV contains only
synthetic rows, including
taxon-specific synthetic concentrations for mugwort/Artemisia, chenopod,
ragweed, grass, alder, ash, birch, cottonwood, cypress, elm, hazel, juniper,
maple, oak, olive, pine, plane, nettle, and total:

```sh
python tools/prepare_pollen_lur.py make-synthetic-demo --output-dir pollen-demo
python tools/prepare_pollen_lur.py fit \
  --stations pollen-demo/synthetic_stations.csv \
  --predictors pollen-demo/synthetic_predictors.nc \
  --output pollen-demo/ragweed_potential.nc \
  --metrics-out pollen-demo/ragweed_metrics.json \
  --taxon ragweed \
  --calibration-factor-m-s 0.001
```

The generator also writes `base_pollen_potential_{taxon}_2025.nc` for every
distinct canonical taxon and `base_pollen_potential_all_taxa_2025.nc`, whose
variables are `{taxon}_base`. These are arbitrary, explicitly synthetic
potentials, not fitted or observed emissions. To place the example mugwort,
total, and all-taxa files at the Earth example paths, run:

```sh
python tools/prepare_pollen_lur.py make-synthetic-demo --output-dir data/pollen --year 2025
```

These fixtures cover a **regional** demonstration grid, not the global UFS
domain. Do not interpret their limited coverage as global data or fill outside
coverage with zero. Prepare/regrid an appropriate global product separately.
All maps use `area_basis="per_grid_cell_area"` and a singleton `time` at January
1 of `--year` (default 2025). `fit --no-output-time` omits that optional time
dimension. The nominal year is a stream coordinate, not a claim of year-specific
observations or a seasonal emission model.

The calibration value above is illustrative and is not a recommended real-world
coefficient. “Artemisia” is normalized to “mugwort”; the other supported taxon
labels are case-insensitive. A CSV containing multiple taxa must be run once per
selected `--taxon`. Output dimensions follow the predictor dimensions, including
structured `lat`/`lon`, curvilinear `tile`/`y`/`x`, or unstructured `ncol`
grids.

## Input contract

### Station CSV

Required columns:

| Column | Contract |
| --- | --- |
| `site_id` | Non-empty station/group identifier; repeated observations at one site share this ID. |
| `latitude`, `longitude` | Finite decimal degrees, within [-90, 90] and [-180, 180]. Rows for a site must have one location. |
| `taxon` | One supported taxon; one model is fitted per canonical taxon. |
| `pollen_count` | Nonnegative pollen concentration in **grains m-3** (not a count per filter or a raw count). |
| `timestamp` | Optional; retained by neither the fit nor its group split. Repeated rows are grouped by `site_id`. |

The file contract is provider-neutral; no Pollen Sense private API or
undocumented export format is assumed. Use a positive, explicit
`--calibration-factor-m-s` coefficient to multiply concentration by m s-1:

```text
base_pollen_potential = pollen_count × calibration_factor_m_s
```

This dimensional conversion produces grains m-2 s-1. The factor is an external
calibration assumption that must be justified for the intended application; it
does not estimate/invert meteorological transport or deposition. The workflow
does not accept another concentration unit and silently convert it.

Use comparable observation windows, sampling durations, and taxon-specific
seasons across stations, with predictors representing a compatible period.
Aggregate concentrations consistently before fitting when records have
different sampling frequencies; repeated rows otherwise give more frequently
sampled sites more weight. Temperature, humidity, precipitation, wind,
boundary-layer mixing, transport, and deposition confound station
concentrations. A scalar calibration alone does **not** remove this
meteorological/transport variability or isolate land-use effects. The current
spatial regression does not normalize weather or model seasonality: evaluate
these limitations separately and avoid interpreting associations as causal
emission estimates.

### Predictor NetCDF

Supply spatial predictors after QA, georeferencing, reprojection/aggregation,
and co-registration to one common grid. Latitude and longitude must be
coordinates/variables named `latitude`/`longitude` or `lat`/`lon`, in decimal
degrees, and may be 1-D rectilinear coordinates or co-located 2-D / unstructured
coordinates. Every selected feature must use the same spatial dimensions and
shape as every other selected feature. Geographic coordinates must be finite.
Predictor values must be finite on every active land cell. A missing active
cell is an error, not a zero-filled value.

An optional co-registered `land_mask` variable must contain only **1 = active
land** and **0 = explicitly masked ocean**. Without a mask, every grid cell is
active. Missing predictors are allowed only at explicitly masked ocean cells;
the model never fits or predicts those cells. Output is zero there by explicit
mask policy, not by imputing missing vegetation. A station whose nearest cell
is masked is rejected (it is not moved to the nearest land cell). Mask policy
and active/masked cell counts are recorded in provenance. Invalid mask values,
an entirely masked domain, or missing active-cell predictors are errors.

By default all variables named `land_cover_*`, `ndvi`, and `canopy_height` are
selected. `--features name1,name2,...` selects an explicit subset. Feature
units are checked for the known feature families:

* `land_cover_*`: fraction/dimensionless (`1`, `fraction`, etc.) or percentage
  (`%`, `percent`); values must lie between 0–1 or 0–100 respectively.
* `ndvi`: dimensionless units (`1` or equivalent).
* `canopy_height`: metres (`m`, `meter`, or `metre`).
* Other explicitly selected predictor variables must declare a non-empty
  `units` attribute. Standardization handles differing scales, not unit
  conversion or inconsistent definitions.

Station points are sampled at the nearest predictor cell using geographic
distance. Stations outside the coordinate extent are rejected, not
extrapolated. Predictor domain and resolution should be appropriate for the
station network; nearest-cell sampling does not correct for representativeness,
urban point-vs-pixel mismatch, or sub-grid land cover.

For land-cover data, independently perform product-specific quality filtering,
decode scale/offset and fill values, aggregate class counts to fractional
cover, and reproject/resample to the target grid using a method appropriate for
categorical data. Compute or select NDVI and canopy height for the same
time/domain and grid. Explicitly mask ocean and other out-of-scope cells before
writing the predictors. Do not encode ocean, cloud, invalid QA, or absent
coverage as zero-valued vegetation predictors. Mark ocean with `land_mask=0`
and retain missing ocean vegetation values if appropriate; missing vegetation
on `land_mask=1` still fails. Do not mask missing land coverage as “ocean” to
hide gaps. The tool is not a satellite QA,
reprojection, aggregation, or predictor-generation utility.

## Fit, held-out evaluation, and output

`fit` uses a scikit-learn `StandardScaler` + Ridge pipeline. A deterministic
`GroupShuffleSplit` holds out complete `site_id` groups for evaluation, then a
fresh pipeline is fitted on all observations for the grid prediction.
Evaluation MAE, RMSE, and (when defined) R² are reported in grains m-2 s-1;
provenance records the split method/seed, held-out and training station/sample
counts, zero overlap assertion, selected features, model, calibration factor,
target units, and nonnegative clipping. The metrics are held-out-site metrics,
not training-fit scores. The final grid model is trained on all supplied
stations after evaluation.

The output NetCDF contains `base_pollen_potential`, preserves the predictor
grid dimensions/coordinates (with the optional singleton `time`), and stores JSON model/evaluation provenance in
the global `provenance_json` attribute. Negative model predictions are clipped
to zero. That physical nonnegativity constraint does not make extrapolated
predictions reliable: Ridge is not an extrapolation model, and values outside
the feature range represented by stations can be unsupported. Inspect feature
coverage/ranges and spatial station coverage before using the output.

## Optional Earthaccess search and download

The optional commands use Earthaccess with user-selected collection, version,
time range, and geographic bounding box; they make no assumptions about a
particular mission/product schema. Install/use the repository's `pollen`
dependency extra and configure Earthdata authentication through the standard
Earthaccess login/authentication flow in the user's environment. No credentials
are stored by this script. Search does not download:

```sh
python tools/prepare_pollen_lur.py earthaccess-search \
  --collection YOUR_COLLECTION_SHORT_NAME \
  --version YOUR_VERSION \
  --start 2024-01-01 --end 2024-12-31 \
  --bbox=-125,24,-66,50
```

Download only when explicitly requested:

```sh
python tools/prepare_pollen_lur.py earthaccess-download \
  --collection YOUR_COLLECTION_SHORT_NAME \
  --version YOUR_VERSION \
  --start 2024-01-01 --end 2024-12-31 \
  --bbox=-125,24,-66,50 \
  --output-dir raw-granules
```

Raw granules require separate product-specific QA, geolocation/CRS handling,
ocean/invalid-cell masking, reprojection and aggregation before they can serve
as the co-registered NetCDF predictors described above.

## Tests

Run the focused offline tests (they generate and remove local synthetic files;
no Earthaccess search or network request is made):

```sh
python -m pytest tests/test_prepare_pollen_lur.py --basetemp=.pollen-lur-pytest
```

Tests use pytest's per-test `tmp_path`; choose a distinct project-local
`--basetemp` for each concurrent invocation rather than sharing that example
directory.
