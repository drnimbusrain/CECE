# Inline Pollen Emissions

## Model

The native Kokkos `pollen` scheme follows CECE's inline MEGAN/FENGSHA pattern:
read a taxon-specific base potential and current host meteorology, evaluate
release controls, and write surface number and mass fluxes. With the scheme's
`refresh_interval_seconds` omitted (or zero), CECE executes it at every base
model timestep. Do not configure a longer refresh interval for online release.

For each surface cell and execution:

$$F = B\, f_T(T)\, f_{RH}(RH)\, f_U(U_{10})\, f_P(P).$$

Here $B$ is `base_pollen_potential` in **grains m-2 s-1 per total grid-cell
area**, already incorporating land-cover/taxon abundance. It is a reference
potential flux, not an annual production, observed concentration, or finite
reservoir. Do not multiply it by vegetation fraction a second time.
CECE writes flux, not an amount per timestep: the host integrates $F\,dt$.
Unchanged meteorology therefore gives identical flux at different timesteps.

### Release functions

Temperature is a thermal proxy for diurnal flowering/dehiscence, not an
astronomical clock or a prognostic seasonal phenology model:

$$f_T = [1+\exp(-s(T_C-T_0))]^{-1},\qquad T_C=T_K-273.15.$$

The selected wind trigger estimates friction velocity from 10 m wind using
the **neutral log law**, zero displacement height, and von Karman constant 0.4:

```math
u_{\ast}=\frac{0.4\,U_{10}}{\ln(10/z_0)},\qquad
f_U=\begin{cases}
0,&u_{\ast}\le u_{\ast,0},\\
1-\exp[-(u_{\ast}-u_{\ast,0})/u_{\ast,s}],&u_{\ast}>u_{\ast,0}.
\end{cases}
```

Do not compare U10 directly to a friction-velocity threshold, map a surface
layer wind at another height as U10, or substitute USTAR as convective velocity.
This approximation does not include atmospheric stability or canopy displacement;
choose and validate roughness and release thresholds for each taxon/domain.

Both RH and precipitation use the same bounded linear suppression:

$$f(x;l,h)=\max\left(0,\min\left(1,\frac{h-x}{h-l}\right)\right),\qquad l<h.$$

The factor is 1 at or below the low threshold $l$, decreases linearly as
$(h-x)/(h-l)$ between the thresholds, and is 0 at or above the high threshold
$h$. For RH, use $x=RH$ with `rh_low_percent` and `rh_high_percent`; for
precipitation, use $x=P$ with `precipitation_low_mm_h` and
`precipitation_high_mm_h`.

Precipitation is a **rate in mm h-1**, independent of `dt`. Convert a host rate
in kg m-2 s-1 (equivalently mm s-1) by multiplying by 3600. If the host
provides interval accumulation, divide by the interval duration in hours first.
Never multiply a rate by the CECE timestep before supplying this field.

All factors are in [0,1]. Base potential is an upper bound, not an annual
emission budget. A static map and temperature alone do **not** represent
dormancy, seasonal depletion, or taxon-specific flowering calendars; for
seasonal applications supply externally calibrated time-varying base maps.
These empirical defaults are demonstration parameters, not a validated global
pollen parameterization.

## Inputs, outputs, and options

| Required logical input | Units | Meaning |
| --- | --- | --- |
| `base_pollen_potential` | grains m-2 s-1 | Nonnegative grid-cell reference potential |
| `temperature` | K | Positive surface/2 m air temperature |
| `relative_humidity` | % | RH in [0,100], not fraction |
| `wind_speed` | m s-1 | Nonnegative 10 m wind magnitude |
| `precipitation_rate` | mm h-1 | Nonnegative precipitation rate |

Required outputs: `pollen_number_emissions` (grains m-2 s-1) and
`pollen_mass_emissions` (kg m-2 s-1). Optional outputs: `pollen_diameter` (m)
and `pollen_density` (kg m-3). All writes are at surface index zero; upper
layers are untouched. Mass assumes spherical grains:
$F_{mass}=F(\pi/6)d^3\rho$. Inert/coarse-mode tracer transport, settling,
deposition and atmospheric processing belong to the host CTM.

| Option | Default | Constraint |
| --- | --- | --- |
| `roughness_length_m` | 0.1 | 0 < z0 < 10 m |
| `friction_velocity_threshold_m_s` | 0.2 | >= 0 |
| `friction_velocity_scale_m_s` | 0.3 | > 0 |
| `temperature_threshold_c` | 10 | Finite |
| `temperature_slope` | 0.35 | > 0, per degree C |
| `rh_low_percent`, `rh_high_percent` | 50, 90 | 0 <= low < high <= 100 |
| `precipitation_low_mm_h`, `precipitation_high_mm_h` | 0, 1 | 0 <= low < high |
| `particle_diameter_um` | 30 | > 0 |
| `particle_density_kg_m3` | 1000 | > 0 |

Missing fields, mismatched horizontal shapes, nonfinite/out-of-range surface
inputs and invalid options raise errors rather than silently producing zero.
The runtime fields carry no enforced units metadata; unit conversion is the
coupler/stream configuration's responsibility.

Use `input_mapping` / `output_mapping` for independent taxon fields. Empty
species layer lists advertise physics-owned outputs without stack overwrites.
Each taxon instance needs a unique registered name. Supported aliases are:
`pollen` (custom), `pollen_mugwort`, `pollen_artemisia` (mugwort synonym),
`pollen_chenopod`, `pollen_ragweed`, `pollen_grass`, `pollen_alder`,
`pollen_ash`, `pollen_birch`, `pollen_cottonwood`, `pollen_cypress`,
`pollen_elm`, `pollen_hazel`, `pollen_juniper`, `pollen_maple`, `pollen_oak`,
`pollen_olive`, `pollen_pine`, `pollen_plane`, `pollen_nettle`, `pollen_total`.
Do not schedule mugwort and Artemisia for the same source. Total pollen is an
aggregate diagnostic; never add it to component taxa in a mass/number budget.
Provider availability and species groupings must be checked, not inferred from
the alias list.

## Pollen Sense and land-use regression

[Station-to-grid LUR preparation](pollen_lur.md) documents licensed CSV
ingestion, land-cover fractions, NDVI, canopy height, Earthaccess acquisition,
held-out-station evaluation, explicit calibration, and free synthetic fixtures.
Station concentrations do not determine emissions without transport/deposition
information. The included regression is a calibrated empirical proxy, **not an
emissions inversion**; raw particle counts first require sampling-volume
normalization to concentration. No Pollen Sense private API or proprietary
observations are bundled or fetched automatically. Confirm account-specific
rights for data export, storage, training, and redistribution.

Apply the fitted spatial relationship across co-registered land predictors,
then supply the static NetCDF to CECE. Spatial infill is performed in
preprocessing, **not by contacting a proprietary service in the model loop**.
Runtime meteorological scaling remains online. Curvilinear cubed-sphere and
unstructured predictor dimensions are supported by the LUR tool; preserve
latitude/longitude and mesh ordering when binding a native host grid, or use
CECE's existing stream regridding for a rectilinear source map. A regional
station network does not by itself validate global extrapolation.

## Examples and validation

Free, offline, 25-hour demonstration:

```sh
python -m pip install -e '.[pollen,test]'
python tools/generate_pollen_example_stream.py
./build/cece_standalone_driver examples/cece_config_pollen.yaml
python -m pytest tests/test_pollen_inline_examples.py tests/test_prepare_pollen_lur.py
```

The generated `data/pollen/pollen_inline_2020.nc` contains arbitrary but
plausible-order reference potentials, time-varying T/U10/RH, and a rain pulse.
It requires neither Pollen Sense rights nor Earthdata access and is marked
synthetic. This validates plumbing, not scientific skill.

Earthdata examples use live MERRA-2 `U10M`/`V10M` magnitude and derived RH from
`M2T1NXSLV`, plus rate-converted `PRECTOTCORR` from `M2T1NXFLX`:

* [Mugwort](../examples/cece_config_earthaccess_pollen.yaml)
* [Aggregate total only](../examples/cece_config_earthaccess_pollen_total.yaml)
* [All distinct registered taxa](../examples/cece_config_earthaccess_pollen_all_taxa.yaml)

Prepare their base maps first using the LUR guide and configure Earthdata
authentication as in [the examples guide](examples.md). No live service or
scientific validation is implied by the offline tests. C++ physics regression
tests in `tests/test_physics.cpp` cover the exact release equation, surface-only
mass conversion, online weather changes, suppression, threshold crossings,
base scaling, field mappings and invalid input failures.

## Migration from the Li/Ambee annual RF approach

This intentionally replaces the native Li et al. (2025) equations, Gaussian
daily pool and autumn forcing. Old options/mappings are rejected with migration
errors. Refit/calibrate a base **flux** product; simply renaming annual
production or dividing it by 86400 is not a defensible unit/science conversion.
Remove annual production, vegetation fraction, DOY/season, sunshine, convective
velocity, autumn-forcing fields and `pollen_phenology_forcing` output. Replace
interval precipitation with `precipitation_rate` and use the options above.
The old RF/Ambee preparation scripts and their tests remain legacy utilities
only; their annual products are not inputs to this inline scheme.
