import importlib.util
import json
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import numpy as np
import pandas as pd
import pytest
import xarray as xr

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "tools" / "prepare_pollen_lur.py"
SPEC = importlib.util.spec_from_file_location("prepare_pollen_lur", SCRIPT)
LUR = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(LUR)


@pytest.fixture
def synthetic_inputs(tmp_path):
    predictors, stations = LUR.generate_synthetic_demo(tmp_path / "demo")
    return predictors, stations


def test_synthetic_fixtures_cover_taxa_and_are_marked_synthetic(synthetic_inputs):
    predictors, stations = synthetic_inputs
    data = pd.read_csv(stations)
    assert set(data["taxon"]) == set(LUR.SYNTHETIC_TAXA)
    assert {"mugwort", "artemisia", "total"} <= set(data["taxon"])
    synonym_rows = data.loc[data["site_id"] == "synthetic-site-00"].set_index("taxon")
    assert synonym_rows.loc["mugwort", "pollen_count"] == synonym_rows.loc[
        "artemisia", "pollen_count"
    ]
    with xr.open_dataset(predictors) as ds:
        assert ds.attrs["synthetic_data"] == "true"
        assert ds["ndvi"].attrs["synthetic"] == "true"


def test_taxon_alias_and_group_holdout_are_station_disjoint(synthetic_inputs):
    predictors, station_path = synthetic_inputs
    stations, sampled, target, ds, features, _ = LUR.load_training_data(
        station_path, predictors, "artemisia", None, 1e-3
    )
    try:
        grid = np.column_stack([ds[name].values.ravel() for name in features])
        _, provenance = LUR.train_evaluate_predict(
            stations, sampled, target, grid, features, "mugwort", 1e-3, random_seed=8
        )
    finally:
        ds.close()
    evaluation = provenance["evaluation"]
    assert provenance["taxon"] == "mugwort"
    assert evaluation["method"] == "GroupShuffleSplit by site_id"
    assert evaluation["station_overlap"] is False
    assert evaluation["training_stations"] + evaluation["heldout_stations"] == 18
    assert evaluation["rmse"] >= 0
    assert evaluation["metrics_units"] == LUR.FLUX_UNITS


def test_cli_end_to_end_writes_nonnegative_flux_and_provenance(synthetic_inputs, tmp_path):
    predictors, stations = synthetic_inputs
    output = tmp_path / "ragweed.nc"
    metrics = tmp_path / "metrics.json"
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "fit",
            "--stations",
            str(stations),
            "--predictors",
            str(predictors),
            "--output",
            str(output),
            "--metrics-out",
            str(metrics),
            "--taxon",
            "ragweed",
            "--calibration-factor-m-s",
            "0.001",
        ],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    with xr.open_dataset(output) as ds:
        assert ds["base_pollen_potential"].dims == ("time", "lat", "lon")
        assert ds.sizes["time"] == 1
        assert ds.time.dt.year.item() == 2025
        assert ds["base_pollen_potential"].attrs["area_basis"] == "per_grid_cell_area"
        assert ds["base_pollen_potential"].attrs["units"] == LUR.FLUX_UNITS
        assert np.isfinite(ds["base_pollen_potential"].values).all()
        assert (ds["base_pollen_potential"].values >= 0).all()
        provenance = json.loads(ds.attrs["provenance_json"])
        assert provenance["taxon"] == "ragweed"
        assert provenance["synthetic_input"] is True
        assert provenance["evaluation"]["station_overlap"] is False
    assert json.loads(metrics.read_text())["calibration_factor_m_s"] == 0.001


def test_missing_predictor_coverage_fails_explicitly(synthetic_inputs, tmp_path):
    predictors, stations = synthetic_inputs
    with xr.open_dataset(predictors) as source:
        ds = source.load()
    ds["ndvi"].values[1, 2] = np.nan
    broken = tmp_path / "missing_predictor.nc"
    ds.to_netcdf(broken)
    ds.close()
    with pytest.raises(LUR.WorkflowError, match="missing or non-finite grid cells"):
        LUR.load_training_data(stations, broken, "ragweed", None, 1e-3)


def test_out_of_domain_station_fails_instead_of_nearest_extrapolation(synthetic_inputs, tmp_path):
    predictors, station_path = synthetic_inputs
    data = pd.read_csv(station_path)
    one = data.loc[data["taxon"] == "ragweed"].iloc[:3].copy()
    one.loc[one.index[0], "latitude"] = 50.0
    csv = tmp_path / "outside.csv"
    one.to_csv(csv, index=False)
    with pytest.raises(LUR.WorkflowError, match="outside the predictor grid"):
        LUR.load_training_data(csv, predictors, "ragweed", None, 1e-3)


def test_invalid_units_features_and_calibration_are_rejected(synthetic_inputs, tmp_path):
    predictors, stations = synthetic_inputs
    with xr.open_dataset(predictors) as source:
        ds = source.load()
    ds["canopy_height"].attrs["units"] = "km"
    bad_units = tmp_path / "bad_units.nc"
    ds.to_netcdf(bad_units)
    ds.close()
    with pytest.raises(LUR.WorkflowError, match="must use metre units"):
        LUR.load_training_data(
            stations,
            bad_units,
            "oak",
            "canopy_height",
            1e-3,
        )
    with xr.open_dataset(predictors) as source:
        ds = source.load()
    ds["latitude"].attrs["units"] = "km"
    bad_coordinates = tmp_path / "bad_coordinate_units.nc"
    ds.to_netcdf(bad_coordinates)
    ds.close()
    with pytest.raises(LUR.WorkflowError, match="Latitude coordinate units"):
        LUR.load_training_data(stations, bad_coordinates, "oak", None, 1e-3)
    with pytest.raises(LUR.WorkflowError, match="not found"):
        LUR.load_training_data(stations, predictors, "oak", "does_not_exist", 1e-3)
    with pytest.raises(LUR.WorkflowError, match="explicitly positive"):
        LUR.load_training_data(stations, predictors, "oak", None, 0)


def test_unstructured_ncol_predictor_grid_preserves_output_dimensions(tmp_path):
    ncol = 12
    n = np.arange(ncol, dtype=float)
    predictor_path = tmp_path / "ncol.nc"
    stations_path = tmp_path / "ncol.csv"
    xr.Dataset(
        {
            "land_cover_tree": ("ncol", 0.1 + n / 20, {"units": "1"}),
            "ndvi": ("ncol", 0.2 + n / 15, {"units": "1"}),
            "canopy_height": ("ncol", 1 + n, {"units": "m"}),
        },
        coords={
            "ncol": np.arange(ncol),
            "latitude": ("ncol", 30 + n * 0.2),
            "longitude": ("ncol", -110 + n * 0.3),
        },
    ).to_netcdf(predictor_path)
    pd.DataFrame(
        {
            "site_id": [f"s{i}" for i in range(ncol)],
            "latitude": 30 + n * 0.2,
            "longitude": -110 + n * 0.3,
            "taxon": ["pine"] * ncol,
            "pollen_count": 100 + n * 10,
        }
    ).to_csv(stations_path, index=False)
    output_path = tmp_path / "ncol_output.nc"
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "fit",
            "--stations",
            str(stations_path),
            "--predictors",
            str(predictor_path),
            "--output",
            str(output_path),
            "--taxon",
            "pine",
            "--calibration-factor-m-s",
            "0.001",
        ],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    with xr.open_dataset(output_path) as result_ds:
        assert result_ds["base_pollen_potential"].dims == ("time", "ncol")
        assert result_ds.sizes["ncol"] == ncol
        assert np.all(result_ds["base_pollen_potential"].values >= 0)


def test_curvilinear_tile_y_x_grid_preserves_output_dimensions(tmp_path):
    tile, rows, columns = np.indices((2, 3, 4))
    latitude = 32 + rows * 0.4 + tile * 0.1
    longitude = -105 + columns * 0.4 + tile * 0.2
    predictor_path = tmp_path / "tile_grid.nc"
    stations_path = tmp_path / "tile_stations.csv"
    xr.Dataset(
        {
            "land_cover_tree": (
                ("tile", "y", "x"),
                0.1 + (rows + columns + tile) / 15,
                {"units": "1"},
            ),
            "ndvi": (("tile", "y", "x"), 0.2 + rows / 10 + columns / 20, {"units": "1"}),
            "canopy_height": (
                ("tile", "y", "x"),
                2 + rows * 2 + columns + tile,
                {"units": "m"},
            ),
        },
        coords={
            "latitude": (("tile", "y", "x"), latitude),
            "longitude": (("tile", "y", "x"), longitude),
        },
    ).to_netcdf(predictor_path)
    pd.DataFrame(
        {
            "site_id": [f"site-{i}" for i in range(latitude.size)],
            "latitude": latitude.ravel(),
            "longitude": longitude.ravel(),
            "taxon": ["birch"] * latitude.size,
            "pollen_count": 100 + np.arange(latitude.size),
        }
    ).to_csv(stations_path, index=False)
    output_path = tmp_path / "tile_output.nc"
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "fit",
            "--stations",
            str(stations_path),
            "--predictors",
            str(predictor_path),
            "--output",
            str(output_path),
            "--taxon",
            "birch",
            "--calibration-factor-m-s",
            "0.001",
        ],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    with xr.open_dataset(output_path) as result_ds:
        field = result_ds["base_pollen_potential"]
        assert field.dims == ("time", "tile", "y", "x")
        assert field.shape == (1, 2, 3, 4)
        assert np.all(field.values >= 0)


def test_station_contract_rejects_wrong_concentration_units(synthetic_inputs):
    predictors, stations = synthetic_inputs
    with pytest.raises(LUR.WorkflowError, match="no implicit unit conversion"):
        LUR.load_training_data(
            stations, predictors, "oak", None, 1e-3, concentration_units="grains/l"
        )


def test_synthetic_cli_emits_all_distinct_taxa_stream_maps(tmp_path):
    output_dir = tmp_path / "pollen"
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "make-synthetic-demo", "--output-dir", str(output_dir)],
        cwd=ROOT, text=True, capture_output=True, check=False,
    )
    assert result.returncode == 0, result.stderr
    with xr.open_dataset(output_dir / "base_pollen_potential_all_taxa_2025.nc") as combined:
        assert set(combined.data_vars) == {
            f"{taxon}_base" for taxon in set(LUR.TAXON_ALIASES.values())
        }
        assert combined.attrs["synthetic_data"] == "true"
        assert combined.sizes["time"] == 1
        assert combined.time.dt.year.item() == 2025
        assert not np.array_equal(combined["mugwort_base"], combined["total_base"])
        for taxon in ("mugwort", "total"):
            with xr.open_dataset(output_dir / f"base_pollen_potential_{taxon}_2025.nc") as single:
                assert single["base_pollen_potential"].attrs["area_basis"] == "per_grid_cell_area"
                assert single["base_pollen_potential"].attrs["units"] == LUR.FLUX_UNITS
                np.testing.assert_equal(single["base_pollen_potential"].values, combined[f"{taxon}_base"].values)
                assert json.loads(single.attrs["provenance_json"])["synthetic_output"] is True


def test_explicit_mask_allows_ocean_missing_predictors_and_rejects_ocean_stations(
    synthetic_inputs, tmp_path
):
    predictors, station_path = synthetic_inputs
    with xr.open_dataset(predictors) as source:
        ds = source.load()
    mask = np.ones(ds["ndvi"].shape, dtype=int)
    mask[0, 0] = 0
    ds["land_mask"] = (("lat", "lon"), mask)
    for name in LUR._select_features(ds, None):
        ds[name].values[0, 0] = np.nan
    masked = tmp_path / "masked.nc"
    ds.to_netcdf(masked)
    ds.close()
    with pytest.raises(LUR.WorkflowError, match="masked ocean cell"):
        LUR.load_training_data(station_path, masked, "grass", None, 0.001)
    stations = pd.read_csv(station_path)
    stations = stations.loc[stations["site_id"] != "synthetic-site-00"]
    land_stations = tmp_path / "land_stations.csv"
    stations.to_csv(land_stations, index=False)
    output = tmp_path / "masked_output.nc"
    assert LUR.main([
        "fit", "--stations", str(land_stations), "--predictors", str(masked),
        "--output", str(output), "--taxon", "grass", "--calibration-factor-m-s", "0.001",
        "--year", "2026",
    ]) == 0
    with xr.open_dataset(output) as result:
        assert result.time.dt.year.item() == 2026
        values = result["base_pollen_potential"].values[0]
        assert values[0, 0] == 0
        assert np.isfinite(values).all()
        assert (values[mask == 1] > 0).all()
        provenance = json.loads(result.attrs["provenance_json"])
        assert provenance["land_mask"]["masked_ocean_cells"] == 1
        assert provenance["land_mask"]["provided"] is True
    with xr.open_dataset(masked) as source:
        invalid = source.load()
    invalid["land_mask"].values[0, 1] = 2
    invalid_path = tmp_path / "invalid_mask.nc"
    invalid.to_netcdf(invalid_path)
    invalid.close()
    with pytest.raises(LUR.WorkflowError, match="only finite 0"):
        LUR.load_training_data(land_stations, invalid_path, "grass", None, 0.001)
    with xr.open_dataset(masked) as source:
        missing_land = source.load()
    missing_land["ndvi"].values[0, 1] = np.nan
    missing_land_path = tmp_path / "missing_active_land.nc"
    missing_land.to_netcdf(missing_land_path)
    missing_land.close()
    with pytest.raises(LUR.WorkflowError, match="missing or non-finite"):
        LUR.load_training_data(land_stations, missing_land_path, "grass", None, 0.001)


def test_predictor_validation_closes_source_on_failure(synthetic_inputs, monkeypatch):
    predictors, _ = synthetic_inputs
    opened = []
    real_open = LUR.xr.open_dataset

    def tracked_open(*args, **kwargs):
        ds = real_open(*args, **kwargs)
        original_close = ds._close
        state = {"closed": False}

        def close():
            state["closed"] = True
            original_close()

        ds.set_close(close)
        opened.append(state)
        return ds

    monkeypatch.setattr(LUR.xr, "open_dataset", tracked_open)
    with pytest.raises(LUR.WorkflowError, match="not found"):
        LUR._read_predictors(predictors, "missing_feature")
    assert opened and all(state["closed"] for state in opened)


def test_time_can_be_omitted(synthetic_inputs, tmp_path):
    predictors, stations = synthetic_inputs
    output = tmp_path / "no_time.nc"
    assert LUR.main([
        "fit", "--stations", str(stations), "--predictors", str(predictors),
        "--output", str(output), "--taxon", "oak", "--calibration-factor-m-s", "0.001",
        "--no-output-time",
    ]) == 0
    with xr.open_dataset(output) as ds:
        assert ds["base_pollen_potential"].dims == ("lat", "lon")


@pytest.mark.parametrize(
    "bbox",
    ["bad", "1,2,3", "1,2,3,4,5", "1,2,nan,4", "-181,0,10,20", "0,-91,10,20", "0,30,10,20"],
)
def test_earthaccess_rejects_invalid_bbox(bbox):
    parser = LUR.build_parser()
    with pytest.raises(SystemExit) as exc:
        parser.parse_args([
            "earthaccess-search", "--collection", "EXAMPLE", "--version", "001",
            "--start", "2025-01-01", "--end", "2025-12-31", f"--bbox={bbox}",
        ])
    assert exc.value.code == 2


@pytest.mark.parametrize("download", [False, True])
def test_earthaccess_search_and_download_are_explicit_and_mocked(
    download, tmp_path, monkeypatch, capsys
):
    results = ["synthetic-granule-1", "synthetic-granule-2"]
    fake = SimpleNamespace(
        login=Mock(),
        search_data=Mock(return_value=results),
        download=Mock(return_value=["synthetic-file-1", "synthetic-file-2"]),
    )
    monkeypatch.setitem(sys.modules, "earthaccess", fake)
    args = [
        "earthaccess-download" if download else "earthaccess-search",
        "--collection", "USER_COLLECTION", "--version", "007",
        "--start", "2025-01-01", "--end", "2025-02-01",
        "--bbox=-120,30,-110,40", "--count", "2",
    ]
    destination = tmp_path / "granules"
    if download:
        args.extend(["--output-dir", str(destination)])
    assert LUR.main(args) == 0
    fake.login.assert_called_once_with()
    fake.search_data.assert_called_once_with(
        short_name="USER_COLLECTION", version="007",
        temporal=("2025-01-01", "2025-02-01"),
        bounding_box=(-120.0, 30.0, -110.0, 40.0), count=2,
    )
    if download:
        fake.download.assert_called_once_with(results, local_path=str(destination))
        assert destination.is_dir()
    else:
        fake.download.assert_not_called()
        assert not destination.exists()
    assert "Found 2 granules" in capsys.readouterr().out
