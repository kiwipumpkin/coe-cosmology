# NOTES — ST+CoE figure-flow draft (2026-09-22)

**PDF:** `/workspace/cosmology/papers/figure_flow_draft_2026-09-22/main.pdf`  
**TeX:** `main.tex` · **Figures:** `figures/` (self-contained copies)  
**Do not overwrite:** `new_draft.tex`, stage_setting paper.

Telegram figure flow ~11:30 CDT 2026-09-22 → section order below.

## Figure-to-file map (included order)

| § | Topic | File in `figures/` | Source |
|---|--------|-------------------|--------|
| 1 | Homologous blast | `hubble_from_coe_relative_velocity.png` | `cosmology/figures/` |
| 2 | Decreasing ρ + m(<r) | `sedov_decreasing_density_toy_journal.png` (+ `.pdf`) | `cosmology/figures/` |
| 3a | Geometry (optional) | `geometry_CoE_observer.png` | `map_compare/amp_mean_param_noise/` |
| 3b | Ideal AF strip | `A_vs_angle_ideal_integration_sheet_AF.png` | same |
| 4 | Axisym + sine analogue | `A_vs_angle_axisymmetric_before_rotation_with_residuals.png` | `cosmology/figures/` |
| 5 | Amp×sin + 2ℓ | `A_vs_angle_after_2ell_with_residuals.png` | `cosmology/figures/` |
| 6 | CoE rings 2MAX+2MIN | `coe_rings_migkas_pin_axis.png` | `cosmology/figures/` |
| 7 | Migkas Fig.8 overlay | `migkas_fig8_theta75_with_coe_rings.png` | `cosmology/figures/` |
| 8a | φ level sets | `rotation_levelsets_phi.png` | `amp_mean_param_noise/` |
| 8b | Amp level sets | `rotation_levelsets_Amp.png` | same |
| 8c | R² table + heatmaps | `R2_table_best_A_B_with_heatmaps.png` | same |
| 8d | Ranked Amp×φ | `rotation_Amp_phi_R2_ranked.png` | same |
| 8e | Shared AB heatmaps | `rotation_Amp_phi_AB_shared_heatmaps.png` | same (present) |
| 9a | Best A + assoc B + 2ℓ | `best_A_vs_angle_with_associated_B_2ell.png` | same |
| 9b | Best B + assoc A + 2ℓ | `best_A_vs_angle_with_associated_B_2ell_fromB.png` | same |
| 10 | Amp family Mollweide | `GRID_amp_family.png` | `levelset_maps/` |
| 11a | Patchwork sources | `map_09_patchwork_sources.png` | `levelset_maps/` |
| 11b | Patchwork 4 draws | `GRID_patchwork_4draws_ALL_MARKERS.png` | same |
| 12 | CMB shell schematic | `coe_geometry_cmb_dipole_schematic.png` (+ `.pdf`) | `cosmology/figures/` |
| 13 | PLACEHOLDER (1/2) | empty `\fbox` | Chris to add |
| 14 | PLACEHOLDER (2/2) | empty `\fbox` | Chris to add |

Also copied for reference: `R2_best_A_B_Amp_phi_summary.csv`.

## Placeholder slots

1. **FIGURE TO ADD (1/2)** — blank caption underline; near end before acknowledgments.  
2. **FIGURE TO ADD (2/2)** — same.

## Verified numbers (sidecars; no invention)

| Item | Value | Source |
|------|-------|--------|
| Axisym R² / RMSE | 0.307→0.553 / 0.059→0.047 | `A_vs_angle_axisymmetric_before_rotation_with_residuals.json` |
| Amp_A, φ (A peak) | −0.06, 11.46° | `R2_best_A_B_Amp_phi_summary.csv` |
| Shared-opp Amp_B | +0.06; R²_A=0.746, R²_B=0.790 | same CSV |
| 2ℓ A-tuned ε, ψ | −0.04325, 175.81° | `best_fit_2ell_multiplier.npz` |
| 2ℓ A-tuned R²_A / RMSE | 0.746→0.875 / 0.038→0.029 | same npz |
| 2ℓ A-tuned R²_B after | **TBD** (CSV has pre-2ℓ only; screenshot flow quoted ~0.663) | — |
| 2ℓ B-tuned | ε≈+0.00643, ψ≈37.71°; R²_B≈0.790, R²_A≈0.766 | `best_fit_2ell_multiplier_from_B.npz` |
| Strip extrema | MAX 30/185; MIN 115/265 | `A_vs_angle_after_2ell_extrema_bands.json`, `coe_rings_*.json` |
| Band width | ±9° = 18° e2e = 5% of 360 | `coe_rings_migkas_pin_axis.json` |
| D/L toy | 3 | `sedov_decreasing_density_toy_journal.json` |
| κ_natural | Ω_m/2=0.15 | NOTES_FOR_CHRIS / captions ontology |
| Migkas A extrema | (17,−9), (274,−22) | `coe_rings_*.json`, `migkas_fig8_*.json` |
| Vision shell | r; Sn→Sn+1 = 2r | `coe_geometry_cmb_dipole_schematic_changes.json` |

## Amp×sin convention (verified, do not invent)

```
A'(ℓ) = (1 − Amp_A · sin(2π ℓ/360 + φ)) · A_physics
```
(ℓ in degrees; φ in radians in the trig arg — sheet/cone_scan/levelset_maps NOTE.)

Associated B: same sine, independent Amp_B (shared-opp often Amp_B = −Amp_A).

## 2ℓ formula (verified from npz)

```
A_after = (1 + ε · sin(2·ℓ_rad + ψ)) · A_before
```

## Fallbacks not used (preferred files present)

- `map_09_patchwork_random.png`
- `GRID_patchwork_4draws_OWN_EXTREMA.png` / `_EXTREMA.png`

## Open / skeptical notes (in PDF tone)

- Amp×sin and 2ℓ are fit templates until derived from time-integrated gravity.
- Bare Gauss+decreasing-ρ leading multipole is quadrupole, not pure dipole.
- CMB deferred as adversary.
- B has no published 2D map extrema (slope flat; Ball-normalized strip).
