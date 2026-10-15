"""Contract checks for the runnable free inline pollen examples."""

from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
import xarray as xr
import yaml

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "tools" / "generate_pollen_example_stream.py"
SPEC = importlib.util.spec_from_file_location("pollen_example", SCRIPT)
GENERATOR = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(GENERATOR)
CONFIG_SPEC = importlib.util.spec_from_file_location(
    "pollen_config_schema", ROOT / "src/python/config.py"
)
CONFIG = importlib.util.module_from_spec(CONFIG_SPEC)
sys.modules[CONFIG_SPEC.name] = CONFIG
CONFIG_SPEC.loader.exec_module(CONFIG)
REQUIRED = {
    "base_pollen_potential",
    "temperature",
    "wind_speed",
    "relative_humidity",
    "precipitation_rate",
}


def test_generator_cli_writes_exact_example_contract(tmp_path):
    output = tmp_path / "pollen.nc"
    subprocess.run(
        [sys.executable, str(SCRIPT), "--output", str(output)],
        check=True,
        capture_output=True,
        text=True,
    )
    config = yaml.safe_load((ROOT / "examples/cece_config_pollen.yaml").read_text())
    stream = config["cece_data"]["streams"][0]
    with xr.open_dataset(output) as dataset:
        assert dict(dataset.sizes) == {"time": 25, "lat": 25, "lon": 30}
        assert dataset.attrs["synthetic"] == "true"
        for mapping in stream["variables"]:
            assert mapping["file"] in dataset
            assert np.isfinite(dataset[mapping["file"]]).all()
        assert dataset.time.values[0] == np.datetime64(config["driver"]["start_time"])
        assert dataset.time.values[-1] == np.datetime64(config["driver"]["end_time"])
        for name in ("artemisia_base", "chenopod_base", "total_base"):
            assert dataset[name].attrs["units"] == "grains m-2 s-1"
            assert dataset[name].attrs["area_basis"] == "per_grid_cell_area"
            assert (dataset[name] >= 0).all()
            assert np.array_equal(dataset[name][0], dataset[name][-1])
        assert dataset.precipitation_rate.attrs["units"] == "mm h-1"
        assert float(dataset.precipitation_rate.isel(time=16).min()) >= 1.0
        assert float(dataset.precipitation_rate.isel(time=12).max()) == 0.0
        assert float(dataset.temperature_2m.max()) > float(dataset.temperature_2m.min())


@pytest.mark.parametrize(
    "name",
    [
        "cece_config_pollen.yaml",
        "cece_config_earthaccess_pollen.yaml",
        "cece_config_earthaccess_pollen_total.yaml",
        "cece_config_earthaccess_pollen_all_taxa.yaml",
    ],
)
def test_all_pollen_examples_resolve_new_fields_at_every_dt(name):
    text = (ROOT / "examples" / name).read_text()
    CONFIG.CeceConfig.from_yaml(text)
    config = yaml.safe_load(text)
    available = set()
    for stream in config["cece_data"]["streams"]:
        variables = stream["variables"]
        if isinstance(variables, list):
            available.update(variable["model"] for variable in variables)
        else:
            available.update(
                variable["model"] if isinstance(variable, dict) else variable
                for variable in variables.values()
            )
        available.update(stream.get("derived_variables", {}))
    names = set()
    for scheme in config["physics_schemes"]:
        assert scheme["name"] not in names
        names.add(scheme["name"])
        assert scheme.get("refresh_interval_seconds", 0) in (
            0,
            config["driver"]["timestep_seconds"],
        )
        options = scheme.get("options", {})
        inputs = options.get("input_mapping", {})
        assert set(inputs) <= REQUIRED
        assert {inputs.get(key, key) for key in REQUIRED} <= available
        outputs = options["output_mapping"]
        assert {"pollen_number_emissions", "pollen_mass_emissions"} <= outputs.keys()
        for field in outputs.values():
            assert config["species"][field] == []
        assert "phenology_method" not in options
    assert not {"pollen_mugwort", "pollen_artemisia"} <= names


def test_earthaccess_examples_use_true_u10_and_rain_rate():
    for config_path in (ROOT / "examples").glob("cece_config_earthaccess_pollen*.yaml"):
        config = yaml.safe_load(config_path.read_text())
        streams = config["cece_data"]["streams"]
        state = next(
            stream for stream in streams if stream.get("short_name") == "M2T1NXSLV"
        )
        assert state["derived_variables"]["wind_speed"] == {
            "method": "wind_speed",
            "u": "U10M",
            "v": "V10M",
        }
        rain = next(
            stream for stream in streams if stream.get("short_name") == "M2T1NXFLX"
        )
        assert rain["variables"]["PRECTOTCORR"] == {
            "model": "precipitation_rate",
            "scale": 3600.0,
        }
        assert "USTAR" not in rain["variables"]
