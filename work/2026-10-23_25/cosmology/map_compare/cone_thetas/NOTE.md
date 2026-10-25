# Cone-averaged A/A_i at Migkas cone sizes (θ = 90°, 75°, 60°, 45°)

**Created:** 2026-09-18 06:14 PM CDT (UTC-0500)

## Deliverables
- `GRID_4x4_cone_thetas.png` — 4×4 grid (rows θ large→small; cols Physics invsq | Sheet invsq | Physics hratio | Sheet hratio)
- `GRID_ours_like_fig8_sheet_invsq.png` — 2×2 Migkas-style (sheet-rotated invsq only)
- `cone_thetas_maps.npz` — 16 cone maps + L/B centers
- `cards/` — optional individual Mollweide cards

## Method
- Full-sky maps from `fullsky_AB/fullsky_AB_maps.npz`
- Cone geometry: same 72×35 L/B centers as `cone75_maps.npz`
- For each θ: pixels with n·n_center ≥ cos(θ), weight-average, then mean-renorm → A/A_i
- Plot: Migkas Galactic (ℓ increases LEFT), magma [0.85, 1.20]

## θ=75° sanity vs cone75_maps.npz
- phys_invsq max|Δ| = 1.689e-05
- sheet_invsq max|Δ| = 1.555e-05
- ours phys [0.920232, 1.022748] vs ref [0.920232, 1.022748]
- ours sheet [0.885448, 1.053989] vs ref [0.885448, 1.053989]

## Ranges (min, max) for all 16 maps

| θ | Physics invsq | Sheet invsq | Physics hratio | Sheet hratio |
|---|---------------|-------------|----------------|--------------|
| 90° | [0.986, 1.014] | [0.954, 1.046] | [0.994, 1.006] | [0.963, 1.037] |
| 75° | [0.920, 1.023] | [0.885, 1.054] | [0.987, 1.039] | [0.955, 1.069] |
| 60° | [0.848, 1.054] | [0.811, 1.091] | [0.969, 1.079] | [0.936, 1.121] |
| 45° | [0.778, 1.090] | [0.746, 1.136] | [0.951, 1.126] | [0.914, 1.181] |

## Morphology vs θ
- **θ=90°** (largest): smoothest, lowest contrast — nearly whole-sky average washes features.
- **θ=75°**: Migkas Fig. 8 matching size; moderate contrast.
- **θ=60°**: sharper hot/cold structure.
- **θ=45°** (smallest): highest contrast / sharpest features (least averaging).

Same trend as Migkas Fig. 8: smaller cones → more contrast.

Done in 9.6s.
