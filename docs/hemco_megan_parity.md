# HEMCO 3.12.1 MEGAN Isoprene Parity Test

This page describes how to run the CECE HEMCO 3.12.1 MEGAN isoprene parity test, which
validates the CECE `megan` and `megan3` biogenic emission schemes against a source-pinned
reference derived from HEMCO 3.12.1 `hcox_megan_mod.F90`.

---

## Overview

The parity implementation provides:

1. **`hemco_megan_stateless.hpp`** — a frozen header-only scalar reference layer containing
   the exact HEMCO 3.12.1 isoprene emission equations and constants.  All constants are
   verified against the pinned HEMCO source (tag `3.12.1`).

2. **`megan_method: hemco_3_12_1`** — a runtime option for the CECE `megan` C++ scheme that
   routes through the frozen stateless reference instead of the native CECE path.

3. **Unit tests** (`tests/test_hemco_megan_runtime.cpp`) that pin every constant,
   verify boundary/branch behavior, and compare a 15-case oracle CSV to double-precision
   tolerance (1×10⁻¹²).

4. **Example YAML configs** and a **comparison script** for running all three schemes
   against each other on the global HEMCO 4°×5° grid.

---

## Science: what the `hemco_3_12_1` mode pins

The mode is called *stateless* because it does not maintain the 15-day temperature average
or the 24-hr PAR average that a coupled HEMCO run would accumulate.  Instead it uses
climatological defaults and accepts externally supplied overrides.

### Constants frozen from HEMCO 3.12.1

| Constant | Value | Notes |
|---|---|---|
| `NORM_FAC` | `1/1.0101081` | Normalisation factor |
| `LDF` (ISOP) | `0.9996` | Light-dependent fraction for isoprene |
| `β` | `0.13 K⁻¹` | Light-independent temperature sensitivity |
| `T_std` | `303.0 K` | Standard temperature |
| `CT1` | `95 kJ mol⁻¹` | Activation energy |
| `CT2` | `200 kJ mol⁻¹` | De-activation energy |
| `CEO` | `2.0` | Empirical scaling at T_opt |
| `T_opt_c1` | `313.0 K` | T_opt intercept |
| `T_opt_c2` | `0.6 K K⁻¹` | T_opt slope |
| `e_opt_coeff` | `0.08 K⁻¹` | e_opt T_avg dependence |
| **`PTOA_C1`** | **`2650 µmol m⁻² s⁻¹`** | **Differs from CECE native (3000)** |
| **`PTOA_C2`** | **`130 µmol m⁻² s⁻¹`** | **Differs from CECE native (99)** |
| **`DOY_offset`** | **`18`** | **Differs from CECE native (10)** |
| `PAR_AVG` | `400 µmol m⁻² s⁻¹` | Direct µmol convention (HEMCO), not W m⁻² × 4.766 |
| `T_avg_15` | `297 K` (default) | Climatological stateless default |
| `γ_CO₂ c₁` | `8.9406` | Possell et al. (2005) |
| `γ_CO₂ c₂` | `0.0024 ppm⁻¹` | Possell et al. (2005) |
| Leaf-age weights | `0.05, 0.60, 1.00, 0.90` | new, growing, mature, old |

### Key differences from CECE native MEGAN defaults

The three most significant formula differences are:

1. **PTOA** — HEMCO uses `PTOA = 2650 + 130·cos(2π·(DOY−18)/365)` while CECE native
   uses `3000 + 99·cos(2π·(DOY−10)/365)`.  This produces a systematically different
   top-of-atmosphere PAR estimate and hence a different `φ` in the PCEEA.

2. **PAR_AVG convention** — HEMCO's 24-hr running average is already in µmol m⁻² s⁻¹,
   so `bbb = 1 + 0.0005 × (PAR_AVG − 400) = 1.0` exactly at the climatological default.
   CECE native starts from a `par_avg = 400 W m⁻²` which is then multiplied by 4.766,
   producing `pac_daily = 1906 µmol m⁻² s⁻¹` and `bbb ≈ 1.75`.

3. **CO₂ reference** — HEMCO 3.12.1 offline runs use 390 ppm; CECE native defaults
   to 400 ppm.

---

## Running the parity test

### 1 — Build and run the unit tests

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc) test_hemco_megan_runtime
ctest -R HEMCO3121Megan -V
```

All constant contract tests and the 15 oracle regression tests must pass.

### 2 — Run HEMCO 3.12.1 parity simulation

```bash
./cece examples/cece_config_hemco_megan_parity.yaml
```

Update the `cece_data` file paths to point to your MERRA-2 4°×5° fields.

### 3 — Run CECE native MEGAN comparison

```bash
./cece examples/cece_config_megan_hemco_comparison.yaml
```

This produces a single NetCDF with both `isoprene_native` and `isoprene_hemco3121` fields.

### 4 — Run CECE MEGAN3 comparison

```bash
./cece examples/cece_config_megan3_hemco_comparison.yaml
```

### 5 — Generate comparison plots

```bash
pip install netCDF4 matplotlib numpy
python scripts/compare_megan_hemco_parity.py \
    --hemco   cece_hemco_megan_parity_4x5.nc \
    --native  cece_megan_comparison.nc \
    --megan3  cece_megan3_hemco_comparison.nc \
    --outdir  plots/
```

This writes six figures and a `summary_stats.txt` to the `plots/` directory.

---

## Configuration reference

```yaml
physics_schemes:
  - name: megan
    language: cpp
    options:
      megan_method: hemco_3_12_1    # Select HEMCO 3.12.1 stateless path
      aef: 1.0e-9                   # AEF [kg m⁻² s⁻¹]
      hemco_co2_ppm: 390.0          # CO₂ concentration [ppm]  (default: 390)
      hemco_par_avg_umol: 400.0     # 24-hr PAR average [µmol m⁻² s⁻¹] (default: 400)
      hemco_t_avg_15_k: 297.0       # 15-day T average [K]  (default: 297)
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `megan_method` | string | `"native"` | `"native"` or `"hemco_3_12_1"` |
| `aef` | double | `1.0e-9` | Isoprene AEF [kg m⁻² s⁻¹] |
| `hemco_co2_ppm` | double | `390.0` | Ambient CO₂ [ppm] |
| `hemco_par_avg_umol` | double | `400.0` | 24-hr PAR average [µmol m⁻² s⁻¹] |
| `hemco_t_avg_15_k` | double | `297.0` | 15-day temperature average [K] |

The `co2_concentration` key is accepted as a synonym for `hemco_co2_ppm` in
`hemco_3_12_1` mode.

---

## Parity scope and limits

| Feature | Supported in `hemco_3_12_1` mode |
|---|---|
| Exact γ_LAI, γ_T_LI, γ_T_LD, γ_PAR, γ_age, γ_CO₂ | ✅ |
| Frozen HEMCO 3.12.1 constants | ✅ |
| HEMCO PTOA formula (2650+130·cos, offset 18) | ✅ |
| HEMCO PAR_AVG µmol convention | ✅ |
| Gridded AEF from MEGAN2.1 netCDF | via `aef` field or gridded input |
| State-evolving 15-day T average | ❌ — fixed at `hemco_t_avg_15_k` |
| State-evolving 24-hr PAR average | ❌ — fixed at `hemco_par_avg_umol` |
| Non-isoprene species | ❌ — use `megan3` scheme |
| γ_SM for ALD2/ETOH | ❌ — isoprene γ_SM = 1 always |

---

## Regenerating the oracle CSV

```bash
python scripts/generate_hemco_megan_oracle.py \
    > tests/data/hemco_megan/hemco_3_12_1_megan_reference.csv
```

Two independent executions must produce byte-identical output before the CSV
may be committed.  Any change to constants triggers a full re-audit of
`hemco_megan_stateless.hpp`.
