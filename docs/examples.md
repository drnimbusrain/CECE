# Using the Provided Examples

CECE includes several example configurations that demonstrate common emission stacking scenarios, modeled after examples in the HEMCO guide and showcasing advanced CECE features.

## Example Scenarios

The `examples/` directory contains several YAML configuration files:

-   `cece_config_ex1.yaml`: Basic single CO species with data stream ingestion
-   `cece_config_ex2.yaml`: Overlaying a regional European CO inventory on a global background
-   `cece_config_ex3.yaml`: Simple testing configuration with minimal grid
-   `cece_config_ex4.yaml`: Using the GFED4 extension for biomass burning
-   `cece_config_ex5.yaml`: Multi-species (CO and NO) emissions with multi-timestep execution
-   `cece_config_ex6.yaml`: Handling non-separated inventories
-   `cece_config_advanced.yaml`: **NEW** - Comprehensive example demonstrating advanced Stacking Engine features
-   `cece_config_earthaccess.yaml`: Cloud-native NASA Earthdata streaming with `earthaccess`

### Advanced Example Highlights

The `cece_config_advanced.yaml` example showcases sophisticated emission processing capabilities:

- **Hierarchical Layer Processing**: Multiple priority levels within categories
- **Temporal Scaling**: Diurnal, weekly, and seasonal emission cycles
- **Vertical Distribution**: Multiple algorithms (PBL, HEIGHT, PRESSURE) for different source types
- **Environmental Dependencies**: Temperature, PAR, and LAI-dependent scaling
- **Geographical Masking**: Land/ocean/vegetation/regional masks
- **Physics Scheme Integration**: Active MEGAN, sea salt, and dust schemes
- **Multi-Source Integration**: Data streams from multiple emission inventories

For complete technical details about how these features work, see the [Stacking Engine Documentation](stacking_engine.md).

---

## Configuration Features by Example

| Example | Grid Size | Species | Key Features |
|---------|-----------|---------|---------------|
| ex1 | 4×4 | CO | Basic data stream integration, simple stacking |
| ex2 | (varies) | CO | Regional override with hierarchy |
| ex3 | 2×2 | CO | Minimal test configuration |
| ex4 | (varies) | Multiple | Biomass burning with GFED4 |
| ex5 | 4×4 | CO, NO | Multi-species, multi-timestep execution |
| ex6 | (varies) | Multiple | Non-separated inventory handling |
| **advanced** | **144×91** | **CO, NOx, Isoprene** | **All advanced features demonstrated** |
| earthaccess | HEMCO 4×5 | Isoprene, soil NO, dust | NASA Earthdata cloud streams via `earthaccess` |

---

## Installing Cloud Extras for Earthdata Tests

The Earthaccess example reads NASA Earthdata cloud-hosted granules directly instead of staging local NetCDF files. Install CECE's optional `cloud` dependencies before running that workflow.

For tests against this source checkout or development branch, install CECE in editable mode from the repository root:

```bash
cd /path/to/CECE
python -m pip install --upgrade pip setuptools wheel
python -m pip install -e '.[cloud,test]'
```

The quoted `'.[cloud,test]'` argument installs the local package plus optional dependency groups declared in `pyproject.toml`:

-   `cloud`: `earthaccess`, `xarray`, `h5netcdf`, `fsspec`, `s3fs`, and `dask`
-   `test`: `pytest`

The quotes are important because many shells treat square brackets as glob characters. Quoting ensures `pip` receives the extras expression unchanged.

If you are installing a published CECE package instead of the local checkout, use:

```bash
python -m pip install 'cece-tools[cloud]'
```

For an isolated local test environment, create and activate a virtual environment first:

```bash
cd /path/to/CECE
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip setuptools wheel
python -m pip install -e '.[cloud,test]'
```

Verify that the cloud stack imports in the same Python environment that will run the tests:

```bash
python - <<'PY'
import earthaccess
import xarray
import h5netcdf
import fsspec
import s3fs
import dask

print("earthaccess:", getattr(earthaccess, "__version__", "installed"))
print("xarray:", xarray.__version__)
print("h5netcdf:", h5netcdf.__version__)
print("fsspec:", fsspec.__version__)
print("s3fs:", s3fs.__version__)
print("dask:", dask.__version__)
PY
```

### Earthdata Credentials

Live Earthdata tests require NASA Earthdata Login credentials. Create an Earthdata account at <https://urs.earthdata.nasa.gov/> if needed, then make the credentials available to `earthaccess`.

The recommended approach for local testing is `~/.netrc`:

```bash
cat > ~/.netrc <<'EOF'
machine urs.earthdata.nasa.gov
	login YOUR_EARTHDATA_USERNAME
	password YOUR_EARTHDATA_PASSWORD_OR_TOKEN
EOF

chmod 600 ~/.netrc
```

Alternatively, export credentials in the current shell:

```bash
export EARTHDATA_USERNAME='your-username'
export EARTHDATA_TOKEN='your-token-or-password'
```

Do not commit credentials, tokens, `.netrc` files, or shell history snippets containing secrets to the repository.

Smoke-test authentication and a small CMR search before running the full workflow:

```bash
python - <<'PY'
import earthaccess

auth = earthaccess.login(strategy="all")
print("Authenticated:", auth.authenticated)

granules = earthaccess.search_data(
		short_name="SPL4SMGP",
		temporal=("2022-07-01", "2022-07-02"),
		count=1,
		cloud_hosted=True,
)
print("Granules found:", len(granules))
PY
```

Run the fast mocked and fixture-backed tests with:

```bash
pytest tests/test_earthaccess_stream_bdsnp_megan3.py -v
```

Run the live Earthdata tests with:

```bash
pytest tests/test_earthaccess_stream_bdsnp_megan3.py -v -m live_earthdata
```

The native standalone driver also uses the cloud extras when a config contains `source: earthaccess` streams. It invokes `scripts/cece_earthaccess_standalone_ingest.py` at each timestep, then injects the returned arrays into the C++ import state before physics runs. Run the driver from the CECE repository root, or set `CECE_EARTHACCESS_HELPER` to the helper path. Set `CECE_PYTHON` if the `earthaccess` environment is not the default `python3`:

```bash
export CECE_PYTHON=/path/to/venv/bin/python
export CECE_EARTHACCESS_HELPER=/path/to/CECE/scripts/cece_earthaccess_standalone_ingest.py
./build/bin/cece_nuopc_driver examples/cece_config_earthaccess_megan3.yaml
```

Earthaccess variable mappings can also apply simple transforms before fields are injected into CECE. For example, MEGAN3 expects `solar_cosine` in the range `[0, 1]`, while some NASA products expose solar zenith angle in degrees. Use `transform: cos_degrees` to convert degrees to daylight cosine during injection:

```yaml
variables:
  solar_zenith_angle:
    model: solar_cosine
    transform: cos_degrees
```

Use `transform: cos_radians` for radian inputs. If `solar_cosine` is mapped without a transform, the bridge validates that the incoming values are already in `[0, 1]`.

If the live test fails to authenticate, confirm that `~/.netrc` is mode `600` and that `python -c 'import earthaccess; print(earthaccess.login(strategy="all").authenticated)'` returns `True` in the active environment.

## Setting Up Examples

To run these examples, you need the associated NetCDF data files. CECE provides a script to automate the setup process.

### 1. Run the Setup Script
```bash
./scripts/setup_hemco_examples.sh
```
This script will:

-   Create the `examples/` directory if it doesn't exist.
-   Copy or generate the example configuration files.
-   Create a `scripts/data_download/` directory with shell scripts to download the required data from S3.

### 2. Download Data
Choose an example to run and execute its data download script:
```bash
# Example: Download data for Example 1
./scripts/data_download/download_ex1.sh
```
This will download the necessary NetCDF files into the `data/` directory.

### 3. Run the Example
You can use the standalone NUOPC driver to run any of the example configurations:
```bash
# Example: Run Example 1
./build/bin/cece_nuopc_driver --config examples/cece_config_ex1.yaml
```
The driver will perform the simulation steps and produce diagnostic output as configured in the YAML file.

---

## Visualizing Example Plans

To better understand the stacking hierarchy of an example, use the visualization utility:
```bash
python scripts/visualize_stack.py --config examples/cece_config_ex2.yaml --output ex2_stack.png
```
This will generate a graph showing how the different layers (global background and regional override) are prioritized and combined.
