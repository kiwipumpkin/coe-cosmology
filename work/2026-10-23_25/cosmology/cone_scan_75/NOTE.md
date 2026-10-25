# θ=75° cone scan of full-sky A/B vs Migkas Fig.6

**Date:** 2026-09-18 01:22 PM CDT (UTC-0500)

## Cone definition (matched to Migkas)

Migkas+2020 §4.3: *“we consider a cone of a given **radius** θ (we use θ = 45°, 60°,
75° and 90°)”* — θ is the **angular separation** from the cone center (not diameter).
We use **θ = 75°**.

Scanning: cone centers on a lat–lon grid (Δℓ=Δb=10°). For each center ĉ, average
map values over pixels with angle(n̂, ĉ) < 75°, with **cos(lat) solid-angle weights**
from the parent full-sky maps. Equal-weight idealization — **not** Migkas cluster
selection or distance-from-center statistical weights.

## Honest labels

| Item | Meaning here |
|---|---|
| **A_hratio** | Physics AF = H_obs/H mean-renorm (primary) |
| **A_invsq** | Physics AF = 1/(H_obs/H)² mean-renorm (secondary) |
| **B** | Mass-proxy AN = m(r_src)/M_enc mean-renorm — **NOT** Migkas L_X–T slope B |
| **sine_*** | Amp×sin sheet / local-density **DIAGNOSTIC** — **not** derived from Gauss |
| **cone avg** | Idealized solid-angle mean — not full cluster pipeline |

## Axes scanned

- Sample preferred: (ℓ, b) = (280°, −20°)
- Joint 3-sample: (ℓ, b) = (303°, −27°)

Physics maps re-rotated from model frame for each axis. Sine diagnostic applied
as sky function of Galactic ℓ after rotation.

## R² vs digitised Migkas Fig.6 (θ=75° cone strip @ Migkas b)

Primary comparison: A_cone(ℓ), B_cone(ℓ) for centers along b = Migkas b, vs Fig.6.
Secondary: Δℓ=90° boxcar of latitude-weighted full-sky longitude profile (known strip match).

| axis | layer | kind | R²_A | R²_B | R²_min |
|---|---|---|---|---|---|
| sample_280m20 | physics | physics | 0.6457 | 0.4360 | 0.4360 |
| joint_303m27 | physics | physics | 0.2915 | 0.6781 | 0.2915 |
| joint_303m27 | sine_sheet | sine_diagnostic | 0.0363 | 0.5365 | 0.0363 |
| joint_303m27 | sine_local_density | sine_diagnostic | 0.0567 | 0.4727 | 0.0567 |
| sample_280m20 | sine_sheet | sine_diagnostic | 0.0090 | 0.3294 | 0.0090 |
| sample_280m20 | sine_local_density | sine_diagnostic | 0.0282 | 0.2691 | 0.0282 |

### Best PHYSICS-ONLY
`physics` @ `sample_280m20` — R²_A=0.6457, R²_B=0.4360

### Best SINE DIAGNOSTIC
`sine_sheet` @ `joint_303m27` — R²_A=0.0363, R²_B=0.5365

### Overall (includes diagnostic)
`physics` @ `sample_280m20` — R²_A=0.6457, R²_B=0.4360

## Cone A/B min/max (over cone-center sky)

| axis | layer | map | min | max |
|---|---|---|---|---|
| sample_280m20 | physics | A_hratio | 0.9844 | 1.0365 |
| sample_280m20 | physics | A_invsq | 0.9241 | 1.0271 |
| sample_280m20 | physics | B | 0.8788 | 1.1225 |
| sample_280m20 | sine_sheet | A_hratio | 0.9526 | 1.0660 |
| sample_280m20 | sine_sheet | A_invsq | 0.8912 | 1.0587 |
| sample_280m20 | sine_sheet | B | 0.8449 | 1.1666 |
| sample_280m20 | sine_local_density | A_hratio | 0.9411 | 1.0688 |
| sample_280m20 | sine_local_density | A_invsq | 0.8891 | 1.0713 |
| sample_280m20 | sine_local_density | B | 0.8431 | 1.1686 |
| joint_303m27 | physics | A_hratio | 0.9845 | 1.0363 |
| joint_303m27 | physics | A_invsq | 0.9248 | 1.0270 |
| joint_303m27 | physics | B | 0.8791 | 1.1222 |
| joint_303m27 | sine_sheet | A_hratio | 0.9442 | 1.0563 |
| joint_303m27 | sine_sheet | A_invsq | 0.9006 | 1.0694 |
| joint_303m27 | sine_sheet | B | 0.8533 | 1.1552 |
| joint_303m27 | sine_local_density | A_hratio | 0.9348 | 1.0558 |
| joint_303m27 | sine_local_density | A_invsq | 0.9023 | 1.0802 |
| joint_303m27 | sine_local_density | B | 0.8546 | 1.1532 |

## Blocker

Migkas **2D** θ=75° cone A/B significance maps (Figs. 8–9, 21–23) are **not digitised**.
This deliverable compares to **Fig.6 longitude strips** + produces our own cone Mollweides.
Digitising those 2D maps (or re-running on the public cluster table) is the next step for
a true apples-to-apples cone↔cone comparison.

## Recommendation for letter

Overall best Fig.6-strip match: physics @ sample_280m20 (R²_A=0.646, R²_B=0.436). Best PHYSICS-ONLY: physics @ sample_280m20 (R²_A=0.646, R²_B=0.436) — expect modest R² (quadrupole-led, two troughs). Best SINE DIAGNOSTIC: sine_sheet @ joint_303m27 (R²_A=0.036, R²_B=0.537) — Amp×sin is NOT first-principles; for letter use physics cone Mollweides as primary geometry result and sine only as labeled diagnostic vs Fig.6. Blocker: Migkas 2D θ=75° cone maps not digitised — strip R² is secondary check.

## Regenerable

```bash
/workspace/cosmology/.venv/bin/python /workspace/cosmology/cone_scan_75/run_cone_scan.py
```

Inputs: `/workspace/cosmology/fullsky_AB/fullsky_AB_maps.npz`,
`/workspace/cosmology/migkas/fig6_digitized.csv`.

## Amplitude washout (important)

At best physics @ (280°, −20°), θ=75° cone A peak-to-peak is **0.0519** vs
Fig.6 **0.3102** (ratio **0.17**). Cone averaging of the
quadrupole strongly suppresses A contrast; reported R² measures **shape**
correlation only, not amplitude match.

B cone ptp **0.2156** vs Fig.6 **0.1846** (ratio **1.17**).

## Sine diagnostic transfer

Amp×sin parameters tuned on 1D circle + Δℓ=90° boxcar **do not transfer** to
full-sky θ=75° cone averages (sine R²_A ≈ 0). Physics-only wins the Fig.6
strip comparison among all layers tested. Keep sine maps labeled DIAGNOSTIC only.

## Updated recommendation

LETTER: lead with PHYSICS θ=75° cone Mollweides (geometry result); best Fig.6-strip shape match is physics @ sample (280°,-20°) (R²_A=0.646, R²_B=0.436). Caveat: cone A amplitude washed to 0.17× Fig.6 ptp (quadrupole averaged over 75°). Sine Amp×sin is retained only as labeled DIAGNOSTIC — it loses Fig.6 match under cone averaging (best sine R²_A=0.036). Prefer sample axis (280°,-20°) over joint (303°,-27°) for A; joint wins B alone. Blocker: Migkas 2D cone maps (Figs 8–9) not digitised.
