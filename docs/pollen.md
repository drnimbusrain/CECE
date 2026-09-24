# Pollen Emissions

## Overview

The `pollen` scheme implements equations 1-9 from Li et al. (2025) as a native Kokkos surface-emission parameterization. It converts taxon-specific vegetation, RF annual production, phenology, and hourly meteorology into inert number and mass fluxes suitable for coarse-particle CTM tracers.

Reference: Li, J. et al. (2025), *Construction and application of a pollen emissions model based on phenology and random forests*, Atmos. Chem. Phys., 25, 3583-3602, https://doi.org/10.5194/acp-25-3583-2025.

Registration names are `pollen`, `pollen_artemisia`, `pollen_chenopod`, and `pollen_total`. The aliases allow CECE to schedule several independently mapped taxa.

## Modeling Workflow

1. Train the offline RF against station/year pollen production and meteorological predictors. `tools/train_pollen_rf.py` writes the gridded `annual_pollen_production` field used by CECE.
2. Supply gridded vegetation fraction, RF production, and start/end DOY. Start/end dates may be generated with the paper's `Rs1`, `Rs2`, or `Rssig` autumn forcing and cumulative threshold.
3. CECE computes the Gaussian daily pollen pool and applies hourly wind, precipitation, RH, and temperature release controls.
4. Map each number/mass pair to inert, non-reactive CTM tracers. Diameter and density outputs preserve settling and coarse-mode metadata needed by CMAQ-like consumers.

## Equations

For taxon $i$, the daily potential is

$$E_i(t)=f_i P_{annual,i}\exp\left[-\frac{(t-\mu)^2}{2\delta^2}\right],\quad
\mu=\frac{sDOY+eDOY}{2},\quad \delta=\frac{eDOY-sDOY}{a}.$$

The default is $a=4$. The three autumn forcing options are

$$R_{s1}=(T_{base}-T_i)^x\left(\frac{L_i}{L_{base}}\right)^y,$$

$$R_{s2}=(T_{base}-T_i)^x\left(1-\frac{L_i}{L_{base}}\right)^y,$$

when $T_i<T_{base}$ and $L_i<L_{base}$, and zero otherwise, and

$$R_{ssig}=\frac{1}{1+\exp(aT_iL_i-b)}.$$

The autumn trigger occurs on the first $t_n$ satisfying $\sum_{t=t_0}^{t_n}R_s(t)\ge Y$. CECE can enforce this with `use_autumn_trigger` and the `autumn_accumulated_forcing` input.

Hourly mobilization follows

$$E_{pollen,i}=E_i f_w f_r f_h,$$

where $f_w=1.5-\exp[-(u_{10}+u_{conv})/5]$. Rain and RH factors are one below their low threshold, zero above their high threshold, and decrease linearly between them. Precipitation is the accumulation in millimeters over one physics execution interval; configure `precipitation_low_mm_interval` and `precipitation_high_mm_interval` when that interval differs from the paper's hourly application. CECE additionally applies a configurable logistic temperature dehiscence factor; set `temperature_slope: 0` for a constant factor or use a small threshold to effectively disable temperature gating.

The daily number potential is divided by 86400 to produce `grains m-2 s-1`. Mass flux assumes spherical grains:

$$m_{grain}=\frac{\pi}{6}d^3\rho.$$

## Fields

Required inputs are `day_of_year`, `vegetation_fraction`, `annual_pollen_production`, `wind_speed`, `convective_velocity`, `precipitation`, `relative_humidity`, `temperature`, and `sunshine_hours`. Optional spatial inputs `season_start_doy` and `season_end_doy` override configured constants. `autumn_accumulated_forcing` is required only when the autumn trigger is enabled.

Required outputs are `pollen_number_emissions` and `pollen_mass_emissions`. Optional outputs are `pollen_diameter`, `pollen_density`, and `pollen_phenology_forcing`. All emissions are written only to the surface layer.

Use `input_mapping` and `output_mapping` to bind taxon-specific names. See `examples/cece_config_pollen.yaml` for Artemisia, chenopod, and total-pollen mappings and CTM metadata.

## RF Mapping Tool

Install the optional dependencies and generate a grid:

```bash
python -m pip install -e '.[pollen]'
python tools/train_pollen_rf.py observations.csv meteorology_grid.nc artemisia_rf.nc \
  --taxon artemisia --model-output artemisia_rf.joblib
```

The training table and predictor NetCDF must contain the selected feature names. Defaults include temperature, wind, precipitation, RH, sunshine, altitude, and pressure. The tool uses a 4:1 train/test split and cross-validated RF hyperparameter search, matching the paper's workflow.

### Global 2025 Climatology

`tools/prepare_global_pollen_rf.py` prepares a global annual predictor grid from MERRA-2 and joins pollen observations to that grid. NASA Earthdata authentication is required. Configure Earthdata credentials outside CECE using Earthaccess (`~/.netrc`, environment-based login, or its interactive login); do not put credentials in YAML or command history.

```bash
python -m pip install -e '.[pollen]'

python tools/prepare_global_pollen_rf.py download-merra2 \
  --year 2025 \
  --output-dir data/pollen/merra2/2025

python tools/prepare_global_pollen_rf.py aggregate-merra2 \
  --year 2025 \
  --surface-glob 'data/pollen/merra2/2025/MERRA2_*tavg1_2d_slv_Nx*.nc4' \
  --radiation-glob 'data/pollen/merra2/2025/MERRA2_*tavg1_2d_rad_Nx*.nc4' \
  --constant-glob 'data/pollen/merra2/2025/MERRA2_*const_2d_asm_Nx*.nc4' \
  --output data/pollen/merra2_annual_predictors_2025.nc
```

Obtain an authorized 2025 historical pollen-count export from a provider such as Ambee. Normalize its columns to `site_id`, `latitude`, `longitude`, `timestamp`, `taxon`, and `pollen_count`, or pass the corresponding `--*-column` options. Ambee API keys must remain in the provider client or environment and must never be committed.

For station-based API retrieval, prepare `pollen_sites.csv` with `site_id`, `latitude`, and `longitude`, set `AMBEE_API_KEY` directly in the shell, and specify the response paths documented for your Ambee subscription. The public endpoint documentation does not guarantee one universal species-count JSON layout, so `--count-path` is deliberately explicit.

```bash
python tools/prepare_global_pollen_rf.py download-ambee \
  pollen_sites.csv \
  data/pollen/ambee_pollen_history_2025.csv \
  --year 2025 \
  --taxon mugwort \
  --records-path data \
  --timestamp-path updatedAt \
  --count-path YOUR_SUBSCRIPTION_COUNT_PATH
```

Use Ambee's licensed bulk historical product for dense global mapping. Calling a 500 m global grid through the point API would require an impractical number of billable requests.

Airborne pollen concentration is affected by transport and removal and is not identical to source production. Therefore, preparation requires a positive `--concentration-to-production` factor calibrated against source measurements or an inverse transport model. This prevents concentration or pollen-index values from being silently labeled as `grains m-2 yr-1`.

The preparer infers each site's reporting cadence, integrates concentration in concentration-days, and defaults to retaining only sites with at least 75% annual temporal coverage. Adjust `--minimum-coverage-fraction` only when the observation product has a documented seasonal sampling design.

```bash
python tools/prepare_global_pollen_rf.py prepare-training \
  ambee_pollen_history_2025.csv \
  data/pollen/merra2_annual_predictors_2025.nc \
  data/pollen/mugwort_training_2025.csv \
  --year 2025 \
  --taxon mugwort \
  --concentration-to-production CALIBRATED_FACTOR

python tools/train_pollen_rf.py \
  data/pollen/mugwort_training_2025.csv \
  data/pollen/merra2_annual_predictors_2025.nc \
  data/pollen/annual_pollen_production_mugwort_2025.nc \
  --taxon mugwort \
  --year 2025 \
  --training-source 'Ambee historical pollen export, 2025' \
  --meteorology-source 'NASA MERRA-2 M2T1NXSLV, M2T1NXRAD, and M2C0NXASM' \
  --model-output data/pollen/mugwort_rf_2025.joblib
```

Repeat the preparation and training stages for each supported taxon. The Google Maps Pollen API is not suitable for reconstructing calendar year 2025: it provides a rolling forecast of up to five days and a Universal Pollen Index rather than a historical concentration archive. Its values must not be used as annual production observations.

Map each generated climatology directly into the corresponding CECE pollen scheme. For example:

```yaml
cece_data:
  streams:
    - name: MUGWORT_RF_CLIMATOLOGY_2025
      file: "data/pollen/annual_pollen_production_mugwort_2025.nc"
      yearFirst: 2025
      yearLast: 2025
      yearAlign: 2025
      taxmode: extend
      tintalgo: nearest
      mapalgo: bilinear
      variables:
        - file: annual_pollen_production
          model: RF_MUGWORT_PANNUAL

physics_schemes:
  - name: pollen_artemisia
    options:
      input_mapping:
        annual_pollen_production: RF_MUGWORT_PANNUAL
```

The climatology provides the annual pool only. Reanalysis, analysis, or forecast meteorology must separately populate the online fields mapped to `temperature`, `wind_speed`, `convective_velocity`, `precipitation`, and `relative_humidity` at each CECE physics interval.
The RF writer stores that pool on a singleton January 1 time coordinate; use `taxmode: extend` and `tintalgo: nearest` so CECE holds it constant throughout the target year.

### Earthaccess-Driven CECE Example

`examples/cece_config_earthaccess_pollen.yaml` combines the generated 2025 mugwort RF climatology with MERRA-2 meteorology streamed by the native standalone driver's Earthaccess helper. `M2T1NXSLV` supplies temperature, pressure, and specific humidity; CECE derives RH, fractional DOY, and astronomical daylight duration. `M2T1NXFLX` supplies wind, friction velocity, and corrected precipitation, with precipitation converted from `mm s-1` to accumulation over the one-hour physics interval.

The example maps MERRA-2 `USTAR` to `convective_velocity` as an explicitly documented turbulence proxy. For production coupling, map a diagnosed convective velocity scale from the driving forecast model when available.

```bash
python -m pip install -e '.[cloud,pollen]'
export EARTHDATA_TOKEN=<NASA-Earthdata-token>
./build/cece_standalone_driver examples/cece_config_earthaccess_pollen.yaml
```

The example expects `data/pollen/annual_pollen_production_mugwort_2025.nc`, generated by the global 2025 workflow above. Earthdata credentials can alternatively be supplied through `~/.netrc`; see `docs/examples.md` for live and staged Earthaccess operation.

To run the bundled configuration without observational training data, generate its explicitly synthetic demonstration stream:

```bash
python tools/generate_pollen_example_stream.py
```

This writes `data/pollen/pollen_rf_phenology_2020.nc`, the path used by `examples/cece_config_pollen.yaml`. Replace this demonstration file with RF, phenology, land-cover, and meteorological products for scientific simulations.