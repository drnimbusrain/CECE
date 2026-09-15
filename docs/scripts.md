# CECE Scripts and Utilities

CECE provides several Python scripts to facilitate data management, configuration migration, and visualization of the emission stacking process.

## Data Management

### `download_hemco_data.py`
Downloads required emission inventories from the public GEOS-Chem S3 bucket.
```bash
python scripts/download_hemco_data.py --config cece_config.yaml --dest data/
```

### `verify_hemco_data.py`
Validates the integrity of downloaded NetCDF files and ensures all required variables are present.
```bash
python scripts/verify_hemco_data.py --config cece_config.yaml --data-dir data/
```

### `setup_hemco_examples.sh`
Automates the creation of example CECE configuration files and generates download scripts for the associated data.
```bash
./scripts/setup_hemco_examples.sh
```

### `stage_earthaccess_streams.py`
Stages NASA Earthdata-backed `source: earthaccess` streams into CECE's native per-step cache format. Run it on a login or data-transfer node with outbound HTTPS access, then point compute jobs at the staged directory with `CECE_EARTHACCESS_STAGE_DIR`.

The staging helper downloads protected granules into each staged step's `granules/` cache before opening them with `xarray+h5netcdf`, so `source: earthaccess` streams should use NetCDF4/HDF5-readable collections. HDF-EOS `.hdf` granules are detected during preflight and should be pre-converted to NetCDF before compute-node runs.
```bash
python scripts/stage_earthaccess_streams.py \
	--config examples/cece_config_earthaccess_megan3.yaml \
	--stage-dir /scratch/$USER/cece_earthaccess_stage \
	--overwrite

python scripts/stage_earthaccess_streams.py \
	--config examples/cece_config_earthaccess_megan3.yaml \
	--stage-dir /scratch/$USER/cece_earthaccess_stage \
	--preflight-only \
	--check-download-access \
	--auth-strategy netrc

export CECE_EARTHACCESS_STAGE_DIR=/scratch/$USER/cece_earthaccess_stage
```

For `--auth-strategy netrc`, store the Earthdata username and account password in `~/.netrc` and run `unset EARTHDATA_TOKEN` before preflight. `EARTHDATA_TOKEN` is only for a current bearer token; an expired token causes CMR `401 Unauthorized: Token does not exist` responses.

### `ursa_earthaccess_staged_run.slurm`
Example Ursa workflow that prepares the EarthAccess Python environment and stages remote streams on a login node, then submits a Slurm job that consumes the staged cache without network access from compute nodes.
```bash
bash scripts/ursa_earthaccess_staged_run.slurm
```

---

## Configuration Migration

### `hemco_to_cece.py`
Converts legacy HEMCO `.rc` configuration files to the CECE YAML format. It handles:

- Recursive includes (`>>>include`)
- `$ROOT` token replacement
- Mapping scale factors and masks to CECE layers
- Parsing grid and diagnostic definitions from auxiliary files

```bash
python scripts/hemco_to_cece.py HEMCO_Config.rc -o cece_config.yaml
```

---

## Visualization

### `visualize_stack.py`
Generates a visual representation (graph) of the emission stacking hierarchy defined in an CECE configuration file. This is useful for verifying that layers, masks, and scale factors are correctly prioritized.
```bash
python scripts/visualize_stack.py --config cece_config.yaml --output stacking_plan.png
```

### `visualize_optimized_stack.py`
Similar to `visualize_stack.py`, but specifically visualizes the fused kernel plan used by the optimized CECE engine.
```bash
python scripts/visualize_optimized_stack.py --config cece_config.yaml --output optimized_plan.png
```
