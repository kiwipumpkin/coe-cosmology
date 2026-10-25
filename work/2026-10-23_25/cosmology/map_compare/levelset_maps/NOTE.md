# Level-set θ=75° cone A/A_i maps (NO noise)

**Created:** 2026-09-18 03:13 PM CDT (UTC-0500)

## What is physics vs Amp×sin diagnostic

### Physics (ST+Gauss two-scale)
- Profile: twoscale with ρ(D)=Ω_m ρ_crit=0.3 ρ_crit
- Dials: **D/L**, **f_u**, **α**, **d/D**
- AF for these cards: **invsq** — A_raw ∝ 1/(H_obs/H)², mean-renormed on the sky,
  then θ=75° solid-angle cone average, then mean-renorm → A/A_i
- Also stored: **hratio** companion maps in the npz (A_raw = H_obs/H)

### Amp×sin diagnostic (NOT first-principles)
Applied on the sky as a function of Galactic longitude ℓ:

```
A' = (1 − Amp_A · sin(2π ℓ/360 + φ)) · A_physics
```

then cone-averaged and mean-normalized. This is the same sheet/local-density
diagnostic used in `cone_scan_75` / `fullsky_AB` — **not** derived from Gauss.

## Families

### A) Amp×sin diagnostic (baseline physics fixed)
Baseline: D/L=5.9, f_u=0.94, α=1, d/D=0.1 (recommended profile).
Only Amp and φ change — matches the level-set CURVE family on A vs ℓ.

| # | name | Amp | φ |
|---|------|-----|---|
| 01 | amp0_physics | 0 | — |
| 02–07 | ±0.03, ±0.06, ±0.10 | at φ=0.35 rad (sheet phase) |
| 08 | local-density | +0.07 | ≈212.82° |

### B) Contour level-set peaks (profile dials)
| # | peak | dials | Amp/φ |
|---|------|-------|-------|
| 09 | baseline recommended | D/L=5.9, f_u=0.94 | +0.07 / 212.82° |
| 10 | D/L vs f_u | D/L=5.69, f_u=0.98 | +0.02 / 217° (from NOTE) |
| 11 | α vs d/D | α=0.62, d/D=0.22 | physics only (no Amp in NOTE) |
| 12 | D/L vs α | D/L=5.69, α=0.50 | physics only (no Amp in NOTE) |

## Plotting convention
- Migkas Galactic: **ℓ increases LEFT** (`lon_disp = −ℓ_wrapped`)
- magma colorbar fixed **[0.85, 1.20]**
- Cards **1400×900**, full Mollweide visible, ℓ labels 90 / 0 / 270
- Celestial equator via astropy SkyCoord (if available)
- **NO noise** on any map

## Files
- Individual cards: `map_NN_<shortname>.png` (1400×900)
  - `map_01_amp0_physics.png`
  - `map_02_amp_p003_phi035.png`
  - `map_03_amp_m003_phi035.png`
  - `map_04_amp_p006_phi035.png`
  - `map_05_amp_m006_phi035_sheet.png`
  - `map_06_amp_p010_phi035.png`
  - `map_07_amp_m010_phi035.png`
  - `map_08_amp_p007_phi213_local.png`
  - `map_09_peak_baseline_rec.png`
  - `map_10_peak_DL_fu.png`
  - `map_11_peak_alpha_eps.png`
  - `map_12_peak_DL_alpha.png`
- `GRID_amp_family.png` — Amp family summary grid
- `GRID_profile_peaks.png` — profile peaks grid
- `levelset_cone75_maps.npz`, `manifest.json`

Preferred axis for sky rotation: (ℓ, b) = (280°, −20°).
