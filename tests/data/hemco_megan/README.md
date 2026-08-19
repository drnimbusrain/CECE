# HEMCO 3.12.1 MEGAN isoprene scalar reference vectors

`hemco_3_12_1_megan_reference.csv` is an independently generated oracle for the
stateless scalar MEGAN emission activity factor from HEMCO 3.12.1.  It is not
CECE output and does not represent a full gridded HEMCO run.

## Provenance

- Repository: `https://github.com/geoschem/HEMCO`
- Release: `3.12.1`
- Science source: `src/Extensions/hcox_megan_mod.F90`
- Generation method: analytic evaluation of the frozen scalar equations using
  Python 3 with IEEE-754 double precision; results are identical across two
  independent runs.
- Generator script: `scripts/generate_hemco_megan_oracle.py`

## Key constants pinned from HEMCO 3.12.1

| Constant | HEMCO 3.12.1 | CECE `megan` default | Description |
|---|---|---|---|
| `PTOA_C1` | 2650 µmol m⁻² s⁻¹ | 3000 | Top-of-atmosphere PAR baseline |
| `PTOA_C2` | 130 µmol m⁻² s⁻¹ | 99 | TOA PAR seasonal amplitude |
| `DOY_offset` | 18 | 10 | Phase offset in PTOA cosine |
| `PAR_AVG` | 400 µmol m⁻² s⁻¹ (direct) | 400 W m⁻² × 4.766 | 24-hr average PAR convention |
| `CO₂` | pinned (per-case) | configurable | Ambient CO₂ for γ_CO₂ |
| `T_avg_15` | 297 K (stateless default) | 297 K | 15-day temperature average |

## Coverage

The fixture contains 15 cases covering:

- standard daytime conditions (T near T_std, mid-day PAR);
- nighttime (solar cosine ≤ 0 → gamma_PAR = 0);
- zero LAI gate (→ zero emission);
- cold temperature (T ≪ T_std → near-zero emission);
- hot temperature (T ≫ T_std → suppressed light-dependent branch);
- low and high CO₂ (γ_CO₂ sensitivity);
- growing and senescing leaf-age states (γ_age non-unity);
- low and high LAI (γ_LAI shape);
- zero solar cosine (LDF light-independent only);
- high-PAR / high-suncos (γ_PAR > 1);
- T_opt (γ_T_LD at maximum); and
- warm mid-day with CO₂ = 400 ppm.

## Tolerance

The C++ implementation is expected to reproduce each oracle value to within
a relative tolerance of 1 × 10⁻¹² (double-precision rounding budget).
