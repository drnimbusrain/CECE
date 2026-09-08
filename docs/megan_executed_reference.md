# Executed HEMCO MEGAN reference and CECE comparisons

The independent reference is the HEMCO 3.12.1 executable at commit
`07da3c29fd85abc3824cb6288578b0b68c2395a3`. The case is a controlled synthetic,
one-hour cold start on 20 June 2021, 12:00-13:00 UTC, on the ordered global
72 x 46 GEOS grid. It is not a real-meteorology or state-evolution validation.
Numerical results must be reported from the actual paired runtime artifact;
the checked-in scalar regressions alone do not execute HEMCO.

## Matching the inputs that enter the emission calculation

- HEMCO's NetCDF reader stores input fields in REAL(sp), even when the input
  file stores doubles. Promote those rounded effective values for CECE.
- The constant effective AEF is also a REAL(sp) value promoted to double:
  `float(1e-9)` kg isoprene m-2 s-1, not the unrounded double `1e-9`.
- With no restart and LAI normalization disabled, current effective LAI and
  initial previous-day LAI are identical. The arbitrary `LAI_PREV` in the
  synthetic source fixture is not what this HEMCO cold start uses.
- Cold-start histories are T_DAVG = REAL(sp)(288.15 K), PARDR_DAVG = 30 W/m2,
  and PARDF_DAVG = 48 W/m2. The LAI interval is one day and UTC DOY is 171.
- HEMCO uses **two** solar calculations: `ExtState%SUNCOS` gates day/night,
  while `GET_GAMMA_PAR_PCEEA` calls MEGAN `SOLAR_ANGLE`, using local time
  from `HcoClock_GetLocal`. Without a TIMEZONES file, that clock uses
  15-degree longitude bins, not continuous longitude/15 hours.
  `HCO_SUNCOS` alone is therefore not the PAR-response sine.
- For CECE's effective `solar_cosine` input, supply the MEGAN PAR sine when
  the outer HEMCO day/night gate is positive, and a nonpositive value otherwise.
  This combines the upstream gates without changing the emission equations.
  For arbitrary reference cases, preserve/export both upstream fields and
  their local-time settings; do not assume this case's reconstruction applies.
- Restart output is updated after emissions. It is not automatically the
  history state used for the just-completed emission calculation.

HEMCO's mean flux diagnostics cast to REAL(sp), accumulate flux times the
3600-second interval in REAL(sp), then divide by the interval and store
REAL(sp). For this one-step case, compare raw CECE double output as well as
that explicit diagnostic-storage projection. Report whether projection is
exact; do not label raw double fluxes bitwise identical to float diagnostics.

## Native MEGAN and MEGAN3 effective-history controls

These scalar options apply to C++ native `megan` and C++ `megan3`:

| Option | Default | Meaning |
|---|---|---|
| `temperature_history_k` | 297 | Effective historical temperature, K |
| `par_history_wm2` | 400 | Effective total historical PAR, W/m2 |
| `days_between_lai` | 30 | Positive interval for leaf-age calculation, days |
| `day_of_year` | 180 | Integer in 1-366 |
| `leaf_age_uses_temperature_history` | false | Use history temperature for leaf age instead of current temperature |

Defaults preserve the previous behavior. Histories are supplied, not evolved.
HEMCO mode keeps its separate `hemco_*` options. These additions do not change
Fortran bridges or claim full MEGAN3 canopy/PFT/state-history functionality.
MEGAN3 remains selected with `name: megan3`, `language: cpp`; it is not a new
`megan_method` string. Its AEF is kmol class m-2 s-1, whereas `megan` AEF is
kg compound m-2 s-1. Match the configured molecular weight before comparing.
MEGAN3 consumes the existing Soil-NO export, with the tested mass-to-amount
conversion before speciation; it does not run a duplicate Soil-NO algorithm.

Compare native MEGAN and MEGAN3 at both the historical defaults and matched
effective histories. They are diagnostic comparisons, not a second independent
standalone MEGAN3 reference. Remaining differences include LDF, normalization,
and low-solar-elevation behavior and must not automatically be called bugs or
improvements. Full canopy integration, evolving histories/restarts, PFT/AEF
generation, multi-time real meteorology, and coupled testing remain separate
implementation and validation work.
