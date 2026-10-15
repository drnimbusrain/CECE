#!/usr/bin/env python3
"""Generate free synthetic base potentials and hourly meteorology for CECE."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import xarray as xr


def build_dataset() -> xr.Dataset:
    lon = np.linspace(115.05, 117.95, 30)
    lat = np.linspace(39.05, 41.45, 25)
    hours = np.arange(25)
    lon_grid, lat_grid = np.meshgrid(lon, lat)
    peak = np.exp(-(((lon_grid - 116.5) / 0.9) ** 2 + ((lat_grid - 40.2) / 0.7) ** 2))
    diurnal = np.sin(2.0 * np.pi * (hours[:, None, None] - 8.0) / 24.0)
    shape = (len(hours), len(lat), len(lon))
    fields = {
        "artemisia_base": (2.0 + 8.0 * peak, "grains m-2 s-1"),
        "chenopod_base": (1.0 + 5.0 * (1.0 - peak), "grains m-2 s-1"),
        "total_base": (10.0 + 25.0 * peak, "grains m-2 s-1"),
        "temperature_2m": (293.15 + 7.0 * diurnal, "K"),
        "wind_speed_10m": (2.0 + 3.0 * np.maximum(diurnal, 0.0), "m s-1"),
        "relative_humidity_2m": (68.0 - 18.0 * diurnal, "%"),
        "precipitation_rate": (
            np.where(
                (hours[:, None, None] >= 15) & (hours[:, None, None] < 18), 1.2, 0.0
            ),
            "mm h-1",
        ),
    }
    dataset = xr.Dataset(
        {
            name: xr.DataArray(
                np.broadcast_to(values, shape).copy(),
                dims=("time", "lat", "lon"),
                attrs={"units": units},
            )
            for name, (values, units) in fields.items()
        },
        coords={
            "time": np.datetime64("2020-08-25T00:00:00")
            + hours.astype("timedelta64[h]"),
            "lat": lat,
            "lon": lon,
        },
        attrs={
            "Conventions": "CF-1.8",
            "title": "Synthetic inline pollen demonstration",
            "source": "CECE arbitrary demonstration potentials; not Pollen Sense measurements or a calibrated LUR product",
            "synthetic": "true",
        },
    )
    dataset.lat.attrs.update(units="degrees_north", standard_name="latitude")
    dataset.lon.attrs.update(units="degrees_east", standard_name="longitude")
    dataset.time.attrs.update(standard_name="time")
    for name in ("artemisia_base", "chenopod_base", "total_base"):
        dataset[name].attrs.update(area_basis="per_grid_cell_area", synthetic="true")
    return dataset


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, default=Path("data/pollen/pollen_inline_2020.nc")
    )
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    build_dataset().to_netcdf(args.output)
    print(args.output)


if __name__ == "__main__":
    main()
