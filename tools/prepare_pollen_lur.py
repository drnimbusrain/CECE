#!/usr/bin/env python3
"""Prepare a station-to-grid pollen land-use regression product."""

from __future__ import annotations

import argparse
import json
import math
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import xarray as xr
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GroupShuffleSplit
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

CONCENTRATION_UNITS = "grains m-3"
FLUX_UNITS = "grains m-2 s-1"
REQUIRED_STATION_COLUMNS = {
    "site_id",
    "latitude",
    "longitude",
    "taxon",
    "pollen_count",
}
TAXON_ALIASES = {
    "mugwort": "mugwort",
    "artemisia": "mugwort",
    "chenopod": "chenopod",
    "ragweed": "ragweed",
    "grass": "grass",
    "alder": "alder",
    "ash": "ash",
    "birch": "birch",
    "cottonwood": "cottonwood",
    "cypress": "cypress",
    "elm": "elm",
    "hazel": "hazel",
    "juniper": "juniper",
    "maple": "maple",
    "oak": "oak",
    "olive": "olive",
    "pine": "pine",
    "plane": "plane",
    "nettle": "nettle",
    "total": "total",
}
SYNTHETIC_TAXA = tuple(TAXON_ALIASES)


class WorkflowError(ValueError):
    """Raised when an input violates the workflow contract."""


def canonical_taxon(value: str) -> str:
    taxon = str(value).strip().lower()
    try:
        return TAXON_ALIASES[taxon]
    except KeyError as exc:
        raise WorkflowError(
            f"Unsupported taxon {value!r}; supported labels: "
            + ", ".join(SYNTHETIC_TAXA)
        ) from exc


def _validate_coordinates(latitude: np.ndarray, longitude: np.ndarray, label: str) -> None:
    if not np.isfinite(latitude).all() or not np.isfinite(longitude).all():
        raise WorkflowError(f"{label} latitude/longitude coordinates must be finite")
    if ((latitude < -90) | (latitude > 90)).any():
        raise WorkflowError(f"{label} latitude must be within [-90, 90]")
    if ((longitude < -180) | (longitude > 180)).any():
        raise WorkflowError(f"{label} longitude must be within [-180, 180]")


def _get_coordinate(ds: xr.Dataset, names: Sequence[str], label: str) -> xr.DataArray:
    for name in names:
        if name in ds.coords:
            return ds.coords[name]
        if name in ds:
            return ds[name]
    raise WorkflowError(
        f"{label} must provide latitude/longitude coordinates named "
        f"{names[0]!r} or {names[1]!r}"
    )


def _validate_coordinate_units(latitude: xr.DataArray, longitude: xr.DataArray) -> None:
    lat_units = str(latitude.attrs.get("units", "")).strip().lower().replace(" ", "").replace("_", "")
    lon_units = str(longitude.attrs.get("units", "")).strip().lower().replace(" ", "").replace("_", "")
    if lat_units and lat_units not in {
        "degree",
        "degrees",
        "degreenorth",
        "degreesnorth",
        "degreen",
        "degreesn",
    }:
        raise WorkflowError(f"Latitude coordinate units must be degrees, not {latitude.attrs['units']!r}")
    if lon_units and lon_units not in {
        "degree",
        "degrees",
        "degreeeast",
        "degreeseast",
        "degreee",
        "degreese",
    }:
        raise WorkflowError(f"Longitude coordinate units must be degrees, not {longitude.attrs['units']!r}")


def _units(var: xr.DataArray, name: str) -> str:
    units = var.attrs.get("units")
    if not isinstance(units, str) or not units.strip():
        raise WorkflowError(f"Predictor {name!r} must declare a non-empty 'units' attribute")
    return units.strip().lower().replace(" ", "")


def _validate_feature_units(name: str, units: str) -> None:
    dimensionless = {"1", "fraction", "dimensionless", "unitless"}
    fraction_units = dimensionless | {"%", "percent", "percentage"}
    if name.startswith("land_cover_") and units not in fraction_units:
        raise WorkflowError(
            f"Land-cover predictor {name!r} must use fraction/dimensionless or percent units"
        )
    if name == "ndvi" and units not in dimensionless:
        raise WorkflowError("NDVI predictor 'ndvi' must use dimensionless units")
    if name == "canopy_height" and units not in {
        "m",
        "meter",
        "meters",
        "metre",
        "metres",
    }:
        raise WorkflowError("Canopy-height predictor 'canopy_height' must use metre units")


def _select_features(ds: xr.Dataset, requested: str | None) -> list[str]:
    if requested:
        features = [item.strip() for item in requested.split(",") if item.strip()]
        if len(set(features)) != len(features):
            raise WorkflowError("Feature list contains duplicate variable names")
    else:
        features = sorted(
            name
            for name in ds.data_vars
            if name.startswith("land_cover_") or name in {"ndvi", "canopy_height"}
        )
    if not features:
        raise WorkflowError(
            "No predictor features selected; provide --features or land_cover_*, ndvi, "
            "and/or canopy_height variables"
        )
    missing = [name for name in features if name not in ds.data_vars]
    if missing:
        raise WorkflowError(f"Predictor variables not found: {', '.join(missing)}")
    return features


def _active_cells(ds: xr.Dataset, reference: xr.DataArray) -> np.ndarray:
    if "land_mask" not in ds:
        return np.ones(reference.shape, dtype=bool)
    mask = ds["land_mask"]
    if mask.dims != reference.dims or mask.shape != reference.shape:
        raise WorkflowError("land_mask must be co-registered with the predictor spatial dimensions")
    values = np.asarray(mask.values, dtype=float)
    if not np.isfinite(values).all() or not np.isin(values, [0, 1]).all():
        raise WorkflowError("land_mask must contain only finite 0 (ocean) or 1 (active land)")
    if not values.any():
        raise WorkflowError("land_mask must contain at least one active land cell")
    return values.astype(bool)


def _read_predictors(path: Path, feature_arg: str | None) -> tuple[xr.Dataset, list[str], tuple[str, ...]]:
    try:
        with xr.open_dataset(path) as source:
            ds = source.load()
    except Exception as exc:
        raise WorkflowError(f"Could not open predictor NetCDF {path}: {exc}") from exc
    features = _select_features(ds, feature_arg)
    reference = ds[features[0]]
    dims = tuple(reference.dims)
    if not dims:
        ds.close()
        raise WorkflowError("Predictor variables must have at least one spatial dimension")
    if "time" in dims:
        raise WorkflowError("Predictor features must be spatial only; output time is set with --year")
    active = _active_cells(ds, reference)
    for name in features:
        var = ds[name]
        if tuple(var.dims) != dims:
            ds.close()
            raise WorkflowError(
                f"Predictor {name!r} dimensions {var.dims} do not match {features[0]!r} "
                f"dimensions {dims}; predictors must be co-registered"
            )
        units = _units(var, name)
        _validate_feature_units(name, units)
        values = np.asarray(var.values, dtype=float)[active]
        if not np.isfinite(values).all():
            count = int((~np.isfinite(values)).sum())
            ds.close()
            raise WorkflowError(
                f"Predictor {name!r} has {count} missing or non-finite grid cells; "
                "provide complete predictor coverage or explicitly subset the domain"
            )
        if name.startswith("land_cover_"):
            maximum = 100.0 if units in {"%", "percent", "percentage"} else 1.0
            if ((values < 0) | (values > maximum + 1e-6)).any():
                ds.close()
                raise WorkflowError(
                    f"Land-cover fractions in {name!r} must be between 0 and {maximum:g}"
                )
    latitude = _get_coordinate(ds, ("latitude", "lat"), "Predictor NetCDF")
    longitude = _get_coordinate(ds, ("longitude", "lon"), "Predictor NetCDF")
    _validate_coordinate_units(latitude, longitude)
    try:
        lat_b, lon_b = xr.broadcast(latitude, longitude)
        lat_b = lat_b.broadcast_like(reference).transpose(*dims)
        lon_b = lon_b.broadcast_like(reference).transpose(*dims)
    except (ValueError, TypeError) as exc:
        ds.close()
        raise WorkflowError(
            "Predictor latitude/longitude coordinates cannot be co-registered with "
            f"feature dimensions {dims}: {exc}"
        ) from exc
    if set(lat_b.dims) != set(dims) or set(lon_b.dims) != set(dims):
        ds.close()
        raise WorkflowError(
            "Predictor latitude/longitude coordinates must cover the same spatial dimensions "
            f"as the features ({dims})"
        )
    lat_values = np.asarray(lat_b.values, dtype=float)
    lon_values = np.asarray(lon_b.values, dtype=float)
    _validate_coordinates(lat_values, lon_values, "Predictor")
    return ds, features, dims


def _sample_stations(
    stations: pd.DataFrame,
    ds: xr.Dataset,
    features: Sequence[str],
    dims: Sequence[str],
) -> np.ndarray:
    lat_da = _get_coordinate(ds, ("latitude", "lat"), "Predictor NetCDF")
    lon_da = _get_coordinate(ds, ("longitude", "lon"), "Predictor NetCDF")
    lat_b, lon_b = xr.broadcast(lat_da, lon_da)
    reference = ds[features[0]]
    lat_grid = np.asarray(lat_b.broadcast_like(reference).transpose(*dims).values, dtype=float).ravel()
    lon_grid = np.asarray(lon_b.broadcast_like(reference).transpose(*dims).values, dtype=float).ravel()
    station_lat = stations["latitude"].to_numpy(dtype=float)
    station_lon = stations["longitude"].to_numpy(dtype=float)
    _validate_coordinates(station_lat, station_lon, "Station")
    if (
        (station_lat < lat_grid.min()).any()
        or (station_lat > lat_grid.max()).any()
        or (station_lon < lon_grid.min()).any()
        or (station_lon > lon_grid.max()).any()
    ):
        raise WorkflowError(
            "A station lies outside the predictor grid's latitude/longitude extent; "
            "no extrapolation or out-of-domain nearest-cell sampling is performed"
        )

    feature_matrix = np.column_stack(
        [np.asarray(ds[name].values, dtype=float).ravel() for name in features]
    )
    sampled_indices = np.empty(len(stations), dtype=int)
    for row, (lat, lon) in enumerate(zip(station_lat, station_lon)):
        lat1 = np.radians(lat)
        lat2 = np.radians(lat_grid)
        delta_lat = lat2 - lat1
        delta_lon = np.radians(lon_grid - lon)
        delta_lon = (delta_lon + np.pi) % (2 * np.pi) - np.pi
        haversine = (
            np.sin(delta_lat / 2) ** 2
            + np.cos(lat1) * np.cos(lat2) * np.sin(delta_lon / 2) ** 2
        )
        sampled_indices[row] = int(np.argmin(haversine))
    active = _active_cells(ds, reference).ravel()
    if not active[sampled_indices].all():
        raise WorkflowError("A station samples an explicitly masked ocean cell in land_mask")
    return feature_matrix[sampled_indices]


def _validate_calibration(calibration_factor_m_s: float, concentration_units: str) -> None:
    if not math.isfinite(calibration_factor_m_s) or calibration_factor_m_s <= 0:
        raise WorkflowError("Calibration factor must be a finite, explicitly positive value")
    if concentration_units.strip().lower() != CONCENTRATION_UNITS:
        raise WorkflowError(
            f"Station pollen_count must have units {CONCENTRATION_UNITS!r}; "
            "no implicit unit conversion is performed"
        )


def load_training_data(
    stations_path: Path,
    predictors_path: Path,
    taxon: str | None,
    features_arg: str | None,
    calibration_factor_m_s: float,
    concentration_units: str = CONCENTRATION_UNITS,
) -> tuple[pd.DataFrame, np.ndarray, np.ndarray, xr.Dataset, list[str], tuple[str, ...]]:
    _validate_calibration(calibration_factor_m_s, concentration_units)
    try:
        stations = pd.read_csv(stations_path)
    except Exception as exc:
        raise WorkflowError(f"Could not read station CSV {stations_path}: {exc}") from exc
    missing_columns = sorted(REQUIRED_STATION_COLUMNS - set(stations.columns))
    if missing_columns:
        raise WorkflowError(f"Station CSV missing required columns: {', '.join(missing_columns)}")
    if stations.empty:
        raise WorkflowError("Station CSV contains no observations")
    if stations["site_id"].isna().any() or (stations["site_id"].astype(str).str.strip() == "").any():
        raise WorkflowError("Every station row must have a non-empty site_id")
    for column in ("latitude", "longitude", "pollen_count"):
        stations[column] = pd.to_numeric(stations[column], errors="coerce")
        if not np.isfinite(stations[column].to_numpy(dtype=float)).all():
            raise WorkflowError(f"Station column {column!r} must contain only finite numeric values")
    _validate_coordinates(
        stations["latitude"].to_numpy(dtype=float),
        stations["longitude"].to_numpy(dtype=float),
        "Station",
    )
    if (stations["pollen_count"] < 0).any():
        raise WorkflowError("Station pollen_count concentration must be nonnegative")
    normalized_taxa = stations["taxon"].map(canonical_taxon)
    stations["taxon"] = normalized_taxa
    if taxon is not None:
        selected_taxon = canonical_taxon(taxon)
        stations = stations.loc[stations["taxon"] == selected_taxon].copy()
        if stations.empty:
            raise WorkflowError(f"No station rows found for selected taxon {selected_taxon!r}")
    else:
        unique_taxa = sorted(stations["taxon"].unique())
        if len(unique_taxa) != 1:
            raise WorkflowError(
                "Station CSV contains multiple taxa; select one with --taxon to prevent "
                "mixing taxon-specific responses"
            )
    stations["site_id"] = stations["site_id"].astype(str)
    locations_per_site = stations.groupby("site_id")[["latitude", "longitude"]].nunique()
    if (locations_per_site > 1).any().any():
        raise WorkflowError("Rows for one site_id must use a consistent latitude/longitude")
    ds, features, dims = _read_predictors(predictors_path, features_arg)
    try:
        sampled = _sample_stations(stations, ds, features, dims)
        target = stations["pollen_count"].to_numpy(dtype=float) * calibration_factor_m_s
        if not np.isfinite(target).all():
            raise WorkflowError("Calibrated target contains non-finite values")
        return stations.reset_index(drop=True), sampled, target, ds, features, dims
    except Exception:
        ds.close()
        raise


def train_evaluate_predict(
    stations: pd.DataFrame,
    sampled_features: np.ndarray,
    target: np.ndarray,
    grid_features: np.ndarray,
    features: Sequence[str],
    taxon: str,
    calibration_factor_m_s: float,
    holdout_fraction: float = 0.25,
    random_seed: int = 42,
    ridge_alpha: float = 1.0,
) -> tuple[np.ndarray, dict[str, Any]]:
    groups = stations["site_id"].to_numpy()
    unique_sites = np.unique(groups)
    if len(unique_sites) < 2:
        raise WorkflowError("Held-out-by-station evaluation requires at least two distinct sites")
    if not 0 < holdout_fraction < 1:
        raise WorkflowError("Holdout fraction must be strictly between 0 and 1")
    if not math.isfinite(ridge_alpha) or ridge_alpha < 0:
        raise WorkflowError("Ridge alpha must be a finite nonnegative value")
    splitter = GroupShuffleSplit(
        n_splits=1, test_size=holdout_fraction, random_state=random_seed
    )
    train_index, test_index = next(splitter.split(sampled_features, target, groups))
    train_sites = set(groups[train_index])
    test_sites = set(groups[test_index])
    if not train_sites or not test_sites or train_sites & test_sites:
        raise WorkflowError("Could not create a non-overlapping held-out station split")

    model = make_pipeline(StandardScaler(), Ridge(alpha=ridge_alpha))
    model.fit(sampled_features[train_index], target[train_index])
    heldout_prediction = np.maximum(model.predict(sampled_features[test_index]), 0.0)
    metrics = {
        "mae": float(mean_absolute_error(target[test_index], heldout_prediction)),
        "rmse": float(np.sqrt(mean_squared_error(target[test_index], heldout_prediction))),
        "r2": float(r2_score(target[test_index], heldout_prediction))
        if len(test_index) > 1
        else None,
    }
    final_model = make_pipeline(StandardScaler(), Ridge(alpha=ridge_alpha))
    final_model.fit(sampled_features, target)
    prediction = np.maximum(final_model.predict(grid_features), 0.0)
    provenance: dict[str, Any] = {
        "workflow": "station-to-grid land-use regression",
        "taxon": taxon,
        "features": list(features),
        "model": "sklearn StandardScaler + Ridge",
        "ridge_alpha": float(ridge_alpha),
        "calibration_factor_m_s": float(calibration_factor_m_s),
        "input_concentration_units": CONCENTRATION_UNITS,
        "target_units": FLUX_UNITS,
        "evaluation": {
            "method": "GroupShuffleSplit by site_id",
            "holdout_fraction": float(holdout_fraction),
            "random_seed": int(random_seed),
            "training_rows": len(train_index),
            "heldout_rows": len(test_index),
            "training_stations": len(train_sites),
            "heldout_stations": len(test_sites),
            "station_overlap": False,
            "metrics_units": FLUX_UNITS,
            **metrics,
        },
        "final_fit_rows": len(target),
        "prediction_clipping": "negative fitted fluxes clipped to zero",
    }
    return prediction, provenance


def _write_output(
    path: Path,
    predictor_ds: xr.Dataset,
    features: Sequence[str],
    dims: Sequence[str],
    prediction: np.ndarray,
    provenance: dict[str, Any],
    year: int | None = 2025,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    shape = tuple(predictor_ds[features[0]].shape)
    if prediction.size != int(np.prod(shape)):
        raise WorkflowError("Prediction size does not match predictor grid shape")
    coordinates = {
        name: coord
        for name, coord in predictor_ds.coords.items()
        if set(coord.dims).issubset(set(dims))
    }
    for names in (("latitude", "lat"), ("longitude", "lon")):
        for name in names:
            if name in predictor_ds and name not in coordinates:
                coord = predictor_ds[name]
                if set(coord.dims).issubset(set(dims)):
                    coordinates[name] = coord
                break
    out = xr.Dataset(
        {
            "base_pollen_potential": xr.DataArray(
                prediction.reshape(shape),
                dims=tuple(dims),
                attrs={
                    "units": FLUX_UNITS,
                    "area_basis": "per_grid_cell_area",
                    "long_name": "calibrated base pollen emission potential",
                    "nonnegative": "true",
                    "calibration_note": (
                        "Concentration multiplied by an explicit positive calibration factor; "
                        "not an inversion of atmospheric transport or deposition."
                    ),
                },
            )
        },
        coords=coordinates,
        attrs={
            "title": "Station-to-grid pollen land-use regression base potential",
            "provenance_json": json.dumps(provenance, sort_keys=True),
            "workflow": "prepare_pollen_lur.py",
            "area_basis": "per_grid_cell_area",
        },
    )
    if provenance.get("synthetic_input") or provenance.get("synthetic_output"):
        out.attrs["synthetic_data"] = "true"
    if year is not None:
        if not 1 <= year <= 9999:
            raise WorkflowError("Output year must lie within [1, 9999]")
        out = out.expand_dims(time=[np.datetime64(f"{year:04d}-01-01", "s")])
        out["time"].attrs["long_name"] = "nominal year of base emission potential"
    try:
        out.to_netcdf(path)
    finally:
        out.close()


def _synthetic_predictor_arrays(rows: int = 6, columns: int = 7) -> dict[str, np.ndarray]:
    yy, xx = np.meshgrid(np.linspace(0, 1, rows), np.linspace(0, 1, columns), indexing="ij")
    tree = 0.1 + 0.72 * (0.65 * yy + 0.35 * (1 - xx))
    grass = 0.05 + 0.55 * (0.35 * yy + 0.65 * xx)
    urban = np.clip(1.0 - tree - grass, 0.03, 0.8)
    return {
        "land_cover_tree": tree,
        "land_cover_grass": grass,
        "land_cover_urban": urban,
        "ndvi": 0.15 + 0.75 * (0.55 * yy + 0.45 * xx),
        "canopy_height": 2.0 + 23.0 * (0.7 * yy + 0.3 * (1 - xx)),
    }


def generate_synthetic_demo(
    output_dir: Path, seed: int = 17, year: int = 2025
) -> tuple[Path, Path]:
    """Generate openly reproducible, explicitly synthetic predictor and station files."""
    output_dir.mkdir(parents=True, exist_ok=True)
    rows, columns = 6, 7
    latitudes = np.linspace(35.0, 40.0, rows)
    longitudes = np.linspace(-120.0, -112.0, columns)
    arrays = _synthetic_predictor_arrays(rows, columns)
    data_vars = {
        name: xr.DataArray(
            values,
            dims=("lat", "lon"),
            attrs={
                "units": "m" if name == "canopy_height" else "1",
                "synthetic": "true",
            },
        )
        for name, values in arrays.items()
    }
    predictor_path = output_dir / "synthetic_predictors.nc"
    xr.Dataset(
        data_vars,
        coords={
            "latitude": ("lat", latitudes, {"units": "degrees_north"}),
            "longitude": ("lon", longitudes, {"units": "degrees_east"}),
        },
        attrs={
            "synthetic_data": "true",
            "description": "Freely generated synthetic predictors; no satellite source data.",
            "generator_seed": int(seed),
        },
    ).to_netcdf(predictor_path)

    rng = np.random.default_rng(seed)
    calibration_factor = 1e-3
    labels = list(SYNTHETIC_TAXA)
    canonical_labels = list(dict.fromkeys(TAXON_ALIASES.values()))
    feature_stack = np.column_stack([arrays[name].ravel() for name in arrays])
    records: list[dict[str, Any]] = []
    for site_number, grid_index in enumerate(
        np.linspace(0, rows * columns - 1, 18, dtype=int)
    ):
        row, column = divmod(int(grid_index), columns)
        values = feature_stack[grid_index]
        site_id = f"synthetic-site-{site_number:02d}"
        taxon_noise: dict[str, float] = {}
        for label in labels:
            canonical = canonical_taxon(label)
            canonical_index = canonical_labels.index(canonical)
            profile = (canonical_index + 1) / len(canonical_labels)
            weights = np.array(
                [
                    0.2 + 0.8 * profile,
                    0.9 - 0.55 * profile,
                    0.15 + (canonical_index % 5) * 0.08,
                    0.35 + (canonical_index % 7) * 0.06,
                    0.1 + (canonical_index % 4) * 0.04,
                ]
            )
            flux = (0.04 + profile * 0.3) + float(np.dot(values, weights)) * (0.2 + profile)
            if canonical not in taxon_noise:
                taxon_noise[canonical] = float(rng.normal(0.0, 0.012))
            flux *= 1.0 + taxon_noise[canonical]
            records.append(
                {
                    "site_id": site_id,
                    "latitude": float(latitudes[row]),
                    "longitude": float(longitudes[column]),
                    "taxon": label,
                    "pollen_count": max(flux, 0.001) / calibration_factor,
                }
            )
    station_path = output_dir / "synthetic_stations.csv"
    pd.DataFrame.from_records(records).to_csv(station_path, index=False)
    features = list(arrays)
    with xr.open_dataset(predictor_path) as source:
        predictor_ds = source.load()
    combined_fields = {}
    for taxon_index, taxon in enumerate(canonical_labels):
        # Arbitrary synthetic maps exercise the stream contract, not fitted emissions.
        potential = (0.02 + 0.01 * taxon_index) * (
            1 + arrays["ndvi"] + arrays["canopy_height"] / 25
        )
        single_path = output_dir / f"base_pollen_potential_{taxon}_{year}.nc"
        provenance = {
            "synthetic_input": True,
            "synthetic_output": True,
            "taxon": taxon,
            "generator_seed": seed,
            "description": "Freely generated regional demonstration; arbitrary potentials, not observations or a fitted model",
            "area_basis": "per_grid_cell_area",
        }
        _write_output(single_path, predictor_ds, features, ("lat", "lon"), potential.ravel(), provenance, year)
        with xr.open_dataset(single_path) as source:
            combined_fields[f"{taxon}_base"] = source["base_pollen_potential"].load()
    with xr.Dataset(
        combined_fields,
        attrs={
            "synthetic_data": "true",
            "description": "Freely generated regional synthetic potentials; not a global product",
            "area_basis": "per_grid_cell_area",
            "generator_seed": seed,
        },
    ) as combined:
        combined.to_netcdf(output_dir / f"base_pollen_potential_all_taxa_{year}.nc")
    predictor_ds.close()
    return predictor_path, station_path


def _parse_bbox(value: str) -> tuple[float, float, float, float]:
    try:
        bbox = tuple(float(part.strip()) for part in value.split(","))
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "bbox must be west,south,east,north"
        ) from exc
    if len(bbox) != 4:
        raise argparse.ArgumentTypeError("bbox must be west,south,east,north")
    west, south, east, north = bbox
    try:
        _validate_coordinates(np.array([south, north]), np.array([west, east]), "Bounding box")
    except WorkflowError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc
    if south > north:
        raise argparse.ArgumentTypeError("bbox south must not exceed north")
    return bbox  # type: ignore[return-value]


def _earthaccess_search(args: argparse.Namespace, download: bool) -> int:
    if args.count <= 0:
        raise WorkflowError("Earthaccess search count must be positive")
    try:
        import earthaccess
    except ImportError as exc:
        raise WorkflowError(
            "earthaccess is required for this command; install the repository's pollen extra"
        ) from exc
    earthaccess.login()
    results = earthaccess.search_data(
        short_name=args.collection,
        version=args.version,
        temporal=(args.start, args.end),
        bounding_box=args.bbox,
        count=args.count,
    )
    print(f"Found {len(results)} granules for {args.collection} version {args.version}")
    if download:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        downloaded = earthaccess.download(results, local_path=str(args.output_dir))
        print(f"Downloaded {len(downloaded)} granule files to {args.output_dir}")
    else:
        for result in results:
            print(result)
    return 0


def _add_earthaccess_args(parser: argparse.ArgumentParser, download: bool) -> None:
    parser.add_argument("--collection", required=True, help="Earthdata short_name/collection ID")
    parser.add_argument("--version", required=True, help="Exact collection version")
    parser.add_argument("--start", required=True, help="Start date/time (Earthaccess temporal format)")
    parser.add_argument("--end", required=True, help="End date/time")
    parser.add_argument(
        "--bbox", required=True, type=_parse_bbox, help="west,south,east,north in degrees"
    )
    parser.add_argument("--count", type=int, default=100, help="Maximum granules to search")
    if download:
        parser.add_argument("--output-dir", required=True, type=Path, help="Granule download folder")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Fit a station-to-grid pollen LUR and write base_pollen_potential. "
            "Station pollen_count is concentration in grains m-3."
        )
    )
    subparsers = parser.add_subparsers(dest="command")
    workflow = subparsers.add_parser("fit", help="Fit/evaluate the station-to-grid model")
    workflow.add_argument("--stations", type=Path, required=True, help="Station schema CSV")
    workflow.add_argument("--predictors", type=Path, required=True, help="Co-registered predictor NetCDF")
    workflow.add_argument("--output", type=Path, required=True, help="Output NetCDF path")
    workflow.add_argument("--taxon", help="One taxon to model (artemisia is a mugwort synonym)")
    workflow.add_argument("--features", help="Comma-separated predictor variable names")
    workflow.add_argument(
        "--calibration-factor-m-s",
        type=float,
        required=True,
        help="Explicit positive conversion: grains m-3 concentration × m s-1 = grains m-2 s-1",
    )
    workflow.add_argument(
        "--concentration-units",
        default=CONCENTRATION_UNITS,
        help=f"Declared pollen_count units (must be {CONCENTRATION_UNITS})",
    )
    workflow.add_argument("--holdout-fraction", type=float, default=0.25)
    workflow.add_argument("--random-seed", type=int, default=42)
    workflow.add_argument("--ridge-alpha", type=float, default=1.0)
    workflow.add_argument("--metrics-out", type=Path, help="Optional held-out metrics/provenance JSON")
    workflow.add_argument("--year", type=int, default=2025, help="Nominal singleton output year (default 2025)")
    workflow.add_argument("--no-output-time", action="store_true", help="Omit the optional singleton time dimension")

    demo = subparsers.add_parser(
        "make-synthetic-demo", help="Generate free synthetic station and predictor fixtures"
    )
    demo.add_argument("--output-dir", type=Path, required=True)
    demo.add_argument("--seed", type=int, default=17)
    demo.add_argument("--year", type=int, default=2025)

    search = subparsers.add_parser("earthaccess-search", help="Search Earthaccess granules only")
    _add_earthaccess_args(search, download=False)
    download = subparsers.add_parser(
        "earthaccess-download", help="Search then explicitly download Earthaccess granules"
    )
    _add_earthaccess_args(download, download=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command == "make-synthetic-demo":
        predictors, stations = generate_synthetic_demo(args.output_dir, args.seed, args.year)
        print(f"Synthetic predictors: {predictors}")
        print(f"Synthetic stations: {stations}")
        print(f"Synthetic taxon labels: {', '.join(SYNTHETIC_TAXA)}")
        return 0
    if args.command in {"earthaccess-search", "earthaccess-download"}:
        return _earthaccess_search(args, args.command == "earthaccess-download")
    if args.command != "fit":
        parser.print_help()
        return 2
    try:
        stations, sampled, target, ds, features, dims = load_training_data(
            args.stations,
            args.predictors,
            args.taxon,
            args.features,
            args.calibration_factor_m_s,
            args.concentration_units,
        )
        try:
            grid_features = np.column_stack(
                [np.asarray(ds[name].values, dtype=float).ravel() for name in features]
            )
            active = _active_cells(ds, ds[features[0]]).ravel()
            selected_taxon = str(stations["taxon"].iloc[0])
            predictions, provenance = train_evaluate_predict(
                stations,
                sampled,
                target,
                grid_features[active],
                features,
                selected_taxon,
                args.calibration_factor_m_s,
                args.holdout_fraction,
                args.random_seed,
                args.ridge_alpha,
            )
            full_prediction = np.zeros(active.size, dtype=float)
            full_prediction[active] = predictions
            provenance["area_basis"] = "per_grid_cell_area"
            provenance["land_mask"] = {
                "provided": "land_mask" in ds,
                "active_cells": int(active.sum()),
                "masked_ocean_cells": int((~active).sum()),
                "policy": "predict active cells only; zero only explicitly masked ocean cells",
            }
            if str(ds.attrs.get("synthetic_data", "")).lower() == "true":
                provenance["synthetic_input"] = True
                provenance["synthetic_input_description"] = str(
                    ds.attrs.get("description", "Generator-marked synthetic predictors")
                )
            _write_output(
                args.output, ds, features, dims, full_prediction, provenance,
                None if args.no_output_time else args.year,
            )
        finally:
            ds.close()
        if args.metrics_out:
            args.metrics_out.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
        print(
            f"Wrote {args.output} with base_pollen_potential "
            f"({FLUX_UNITS}); features={','.join(features)}"
        )
        print(json.dumps(provenance["evaluation"], sort_keys=True))
        return 0
    except (WorkflowError, OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
