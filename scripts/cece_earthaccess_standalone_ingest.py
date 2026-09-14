#!/usr/bin/env python3
"""Fetch one timestep of CECE earthaccess streams for the standalone driver.

This helper is invoked by the native C++ standalone driver when a YAML config
contains ``cece_data.streams`` entries with ``source: earthaccess``. It opens
NASA Earthdata granules through earthaccess/xarray, applies configured field
transforms, interpolates to the CECE target grid, and writes raw float64 arrays
that the C++ driver injects into the CECE import state.
"""

from __future__ import annotations

import argparse
import re
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable, List, Optional, Tuple

import numpy as np
import yaml


def _load_coords(path: Path) -> np.ndarray:
    return np.asarray([float(line.strip()) for line in path.read_text().splitlines() if line.strip()], dtype=np.float64)


def _bounding_box_from_grid(grid: dict) -> Optional[Tuple[float, float, float, float]]:
    required = ("lon_min", "lon_max", "lat_min", "lat_max")
    if not all(key in grid for key in required):
        return None
    return (float(grid["lon_min"]), float(grid["lat_min"]), float(grid["lon_max"]), float(grid["lat_max"]))


def _earthaccess_streams(config: dict) -> List[dict]:
    grid = config.get("driver", {}).get("grid", {}) or {}
    streams = []
    for stream in config.get("cece_data", {}).get("streams", []) or []:
        if stream.get("source") != "earthaccess":
            continue
        item = dict(stream)
        if item.get("bounding_box") is None:
            item["bounding_box"] = _bounding_box_from_grid(grid)
        streams.append(item)
    return streams


def _mapping_target_and_transform(mapping: Any) -> Tuple[str, Optional[str]]:
    if isinstance(mapping, str):
        return mapping, None
    if isinstance(mapping, dict):
        target = mapping.get("model") or mapping.get("field") or mapping.get("name")
        if not target:
            raise ValueError("earthaccess variable mapping dict must contain 'model'")
        return str(target), mapping.get("transform")
    raise TypeError("earthaccess variable mapping must be a string or mapping dict")


def _apply_transform(values: np.ndarray, transform: Optional[str]) -> np.ndarray:
    if transform is None or transform == "none":
        return values
    if transform == "cos_degrees":
        return np.clip(np.cos(np.deg2rad(values)), 0.0, 1.0)
    if transform == "cos_radians":
        return np.clip(np.cos(values), 0.0, 1.0)
    raise ValueError(f"Unsupported earthaccess variable transform: {transform}")


def _validate_field(field_name: str, values: np.ndarray) -> None:
    if not np.all(np.isfinite(values)):
        raise ValueError(f"earthaccess field {field_name!r} contains non-finite values")
    if field_name == "solar_cosine" and (float(np.nanmin(values)) < -1.0e-12 or float(np.nanmax(values)) > 1.0 + 1.0e-12):
        raise ValueError(
            "earthaccess field 'solar_cosine' must be in [0, 1]; "
            "use transform: cos_degrees or transform: cos_radians for solar zenith angle inputs"
        )


def _coord_name(data_array: Any, candidates: Iterable[str]) -> Optional[str]:
    names = set(data_array.coords) | set(data_array.dims)
    for candidate in candidates:
        if candidate in names:
            return candidate
    lowered = {name.lower(): name for name in names}
    for candidate in candidates:
        match = lowered.get(candidate.lower())
        if match is not None:
            return match
    return None


def _normalize_longitudes(lon: np.ndarray) -> np.ndarray:
    return ((lon + 180.0) % 360.0) - 180.0


def _select_time(data_array: Any, timestamp: np.datetime64) -> Any:
    if "time" in data_array.coords or "time" in data_array.dims:
        return data_array.sel(time=timestamp, method="nearest")
    return data_array


def _interp_to_target(data_array: Any, target_lons: np.ndarray, target_lats: np.ndarray) -> np.ndarray:
    lon_name = _coord_name(data_array, ("lon", "longitude", "LON", "grid_lon", "grid_lont", "x"))
    lat_name = _coord_name(data_array, ("lat", "latitude", "LAT", "grid_lat", "grid_latt", "y"))

    data_array = data_array.squeeze(drop=True)
    if lon_name is None or lat_name is None:
        values = np.asarray(data_array.values, dtype=np.float64).squeeze()
        if values.shape == (target_lats.size, target_lons.size):
            return values
        if values.shape == (target_lons.size, target_lats.size):
            return values.T
        raise ValueError(f"Cannot identify latitude/longitude coordinates for variable {data_array.name!r}")

    if lon_name in data_array.coords:
        lon_values = np.asarray(data_array[lon_name].values, dtype=np.float64)
        if lon_values.ndim == 1:
            data_array = data_array.assign_coords({lon_name: _normalize_longitudes(lon_values)})
            data_array = data_array.sortby(lon_name)

    if lat_name in data_array.coords:
        lat_values = np.asarray(data_array[lat_name].values, dtype=np.float64)
        if lat_values.ndim == 1:
            data_array = data_array.sortby(lat_name)

    interpolated = data_array.sel({lon_name: target_lons, lat_name: target_lats}, method="nearest")
    values = np.asarray(interpolated.values, dtype=np.float64).squeeze()

    expected = (target_lats.size, target_lons.size)
    if values.shape == expected:
        return values
    transposed = (target_lons.size, target_lats.size)
    if values.shape == transposed:
        return values.T
    raise ValueError(f"Interpolated field {data_array.name!r} has shape {values.shape}, expected {expected}")


def _safe_name(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name)


def _open_dataset(stream: dict) -> Any:
    import earthaccess
    import xarray as xr

    earthaccess.login(strategy="all")

    granules = earthaccess.search_data(
        short_name=stream["short_name"],
        temporal=(stream["temporal_start"], stream["temporal_end"]),
        bounding_box=stream.get("bounding_box"),
        version=stream.get("version"),
        cloud_hosted=stream.get("cloud_hosted", True),
        daac=stream.get("daac"),
        count=-1,
    )
    if not granules:
        raise RuntimeError(f"No earthaccess granules found for stream {stream.get('name', '<unnamed>')!r}")

    open_kwargs = {}
    if stream.get("block_size") is not None:
        open_kwargs["block_size"] = stream["block_size"]
    if stream.get("cache_type") is not None:
        open_kwargs["cache_type"] = stream["cache_type"]

    try:
        file_objs = earthaccess.open(granules, **open_kwargs) if open_kwargs else earthaccess.open(granules)
    except TypeError:
        file_objs = earthaccess.open(granules)

    return xr.open_mfdataset(file_objs, engine="h5netcdf", combine="by_coords")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--time", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--lon-file", required=True, type=Path)
    parser.add_argument("--lat-file", required=True, type=Path)
    args = parser.parse_args()

    config = yaml.safe_load(args.config.read_text())
    if not isinstance(config, dict):
        raise ValueError("CECE config must be a YAML mapping")

    streams = _earthaccess_streams(config)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    target_lons = _load_coords(args.lon_file)
    target_lats = _load_coords(args.lat_file)
    timestamp = np.datetime64(datetime.fromisoformat(args.time.replace("Z", "+00:00")).replace(tzinfo=None))

    manifest_lines = []
    for stream in streams:
        dataset = _open_dataset(stream)
        try:
            for nasa_var, mapping in (stream.get("variables") or {}).items():
                field_name, transform = _mapping_target_and_transform(mapping)
                if nasa_var not in dataset:
                    raise KeyError(f"Variable {nasa_var!r} not found in earthaccess stream {stream.get('name', '<unnamed>')!r}")
                selected = _select_time(dataset[nasa_var], timestamp)
                values = _interp_to_target(selected, target_lons, target_lats)
                values = _apply_transform(values, transform)
                values = np.asarray(values, dtype=np.float64)
                _validate_field(field_name, values)

                flat = np.ascontiguousarray(values.reshape(-1))
                binary_path = args.output_dir / f"{_safe_name(field_name)}.f64"
                flat.tofile(binary_path)
                manifest_lines.append(f"{field_name} {target_lons.size} {target_lats.size} 1 {binary_path.name} {float(np.min(flat)):.17g} {float(np.max(flat)):.17g}")
        finally:
            dataset.close()

    (args.output_dir / "manifest.txt").write_text("\n".join(manifest_lines) + ("\n" if manifest_lines else ""))
    print(f"Wrote {len(manifest_lines)} earthaccess field(s) to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
