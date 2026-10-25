# NOTES FOR CHRIS — Stage-setting paper (2026-09-20)

**PDF:** `/workspace/cosmology/papers/stage_setting_2026-09-20/main.pdf`  
**TeX:** `main.tex` · **Figures:** `figures/` (self-contained copies)

## Section ↔ figure map

| § | Topic | Figures included |
|---|--------|------------------|
| 1 | Hubble vector derivation | `hubble_from_coe_relative_velocity.png` |
| 2 | A vs angle model | `geometry_CoE_observer.png`, `toy_density_H_profile.png`, `A_vs_angle_ideal_integration_sheet_AF.png` |
| 3 | Migkas observable | `AB_vs_angle_migkas_fig6.png` (+ `fig6_digitized.csv` copied) |
| 4 | Residuals before rotation | `residual_before_rotation.png` *(generated from AF + digitized CSV)* |
| 5 | Amp×sin rotation | `A_vs_angle_ideal_plus_chris_rotation.png`, `R2_table_best_A_B_with_heatmaps.png`, `rotation_Amp_phi_R2_ranked.png`, `migkas_bars_rotation_overlay.png` |
| 6 | Residuals after rotation | `residual_after_rotation.png` *(from 2ℓ npz A_before)*, `best_fit_residual_vs_angle_with_sincos.png` |
| 7 | 2ℓ residual | `best_fit_with_2ell_multiplier.png`, `best_A_vs_angle_with_associated_B_2ell.png`, `best_A_vs_angle_with_associated_B_2ell_fromB.png`, `residual_after_2ell.png` |
| 8 | Scatter / R² | `migkas_sheet_scatter_amp006.png`, `migkas_bars_rotation_scatter.png` |
| 9 | CMB shell toy | `coe_geometry_cmb_dipole_schematic.png`, `four_scenarios_overview.png`, `scenarios_024_R4r_shift2r.png` |
| 10 | Mollweide / sky | `sky_mollweide.png`, `mollweide_amp_sin_coe_axis.png` *(new dipole-proxy)*, `sbs_SHARED_SCALE_both_mollweide.png`, `CARD_theta75_sheet_hratio.png`, `Cl_A_migkas_smooth_both.png` |
| 11 | Work remaining | (text bullets) |
| 12 | Relativity required | points to `papers/gr_ch2_ltb_2026-09-20/main.pdf` |

## Numbers used (on-disk / your list — no invention)

- Geometry grid: M/M0=3, d/D=0.5, L≈0.169, D/L≈5.9, κ∼0.7, R²≈0.543  
- Ideal continuous npz: d/D≈0.76, D/L≈7.62, κ≈0.458, R²≈0.78  
- Amp×sin: Amp=−0.06, φ≈11.46°, R²_A≈0.746; Amp_B=−Amp_A → R²_B≈0.79  
- 2ℓ A-tuned: ε≈−0.043, ψ≈176°, R²_A≈0.875 (RMSE 0.038→0.029)  
- CMB: T0=2.72548 K; A_dip=3362.08 μK; (δT/T)≈1.23e-3; δF/F≈4.93e-3; R/D≈203  
- CoE axis: (258.5°, +48.3°) / anti (78.5°, −48.3°)

## Open work (also in §11–12 of the PDF)

1. Joint A+B likelihood; 3D axis (ℓ,b).  
2. Retire or derive Amp×sin and 2ℓ from physics.  
3. Time-integrated gravity inside the fit.  
4. Selection / mask / absorption for Migkas-like samples.  
5. Real Mollweide of full ST+Gauss A (current new map is Amp×sin dipole proxy only).  
6. LTB (or cousin) chapter — draft already at `gr_ch2_ltb_2026-09-20`.  
7. CMB: Alnes/Kodama bounds as adversary; shell toy is not closure.

## Generated for this paper

- `residual_before_rotation.png`, `residual_after_rotation.png`, `residual_after_2ell.png`  
- `mollweide_amp_sin_coe_axis.png` (provisional; caption says so)

## Acknowledgments in PDF

Wolfgang Rindler; Luke Ackerman (genesis).

## Google Drive
Uploaded to My Drive as `stage_setting_CoE_Hubble_2026-09-20.pdf`
- id: 1BNb9n6Tz3NGsTI-fAb22lxAbCgOcp4D2
- link: https://drive.google.com/file/d/1BNb9n6Tz3NGsTI-fAb22lxAbCgOcp4D2/view?usp=drivesdk
