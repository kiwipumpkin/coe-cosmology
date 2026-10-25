# CMB ΔT/T on Chris's recombination track (variable R_see)

## Setup
- Reuses `../lib4.py` (Newtonian integrator, same Doppler convention) and termslice shell constants (KIN≈1.028, KOUT≈1.475, Ω_s=0.049). Gravity only, Λ=0. ⟨H⟩=71.6. Compact core (point); Mc re-solved at each track point so cluster quadrupole a2=0.084 at Rs=0.31 Gpc.
- Emission sphere is **observer-centred of radius R_see**, not the old fixed R=14.3 Gpc.
- Angle θ=0 toward CoE, 180 away. ΔT/T_Doppler = −n·(V_emit−V_obs)/c (linear). Static = (Φ_emit−Φ_obs)/c². Kill = −(Doppler − ⟨Doppler⟩).

## Formulas (documented track)
- Continuous Sedov-like / coasting homologous cooling: T ∝ 1/a ∝ 1/t for each shell (R ∝ t).
- Recombination at T_rec = 3000 K; today's monopole T0 = 2.725 K sets the ratio
  - f = T0/T_rec = t_rec/t_now = D_rec/D_now ≈ 9.083×10⁻⁴.
- Age t_now from the integrator: age = Ht × TH0 (TH0≈13.967 Gyr).
- Time since recombination: Δt = t_now − t_rec = t_now (1 − f).
- Vision radius: R_see = c × Δt ≈ 4.22–4.26 Gpc on this grid (nearly flat; set by age, not by D).
- Shell displacement in the same Δt: D_now − D_rec = V_shell × Δt with V_shell = D_now/t_now.
- Lock: R_see / (D_now − D_rec) = c / V_shell. Fraction of (R+D) we can see: R_see/(R_see+D_now).
- Shell rule on R_see: R_in = D − KIN·R_see, R_out = D + KOUT·R_see. Shell minimum D_min = KIN·R_see ≈ 4.34–4.38 Gpc. No track point was run below D_min.
- Scan parameter: D_now (equivalently D_rec = f·D_now). Endpoints: D_now = 4.5 → 94 Gpc, D_rec = 0.00409 → 0.0854 Gpc, R_see = 4.225 → 4.259 Gpc.

## Kill a2 vs D_rec (cluster fit kept)
- D_now=4.5 Gpc (D_rec=4.09 Mpc, near shell-min): kill a2 = +0.0220; peak |dE/E| = 10.0%; |static a2| = 9.5e−4. Note: vrmin<0 (turnaround oscillations near CoE axis).
- D_now=5 Gpc: kill a2 = +0.0562; peak |dE/E| = 28.9%.
- D_now=6 Gpc: kill a2 = +0.1087; peak |dE/E| = 69.9% (largest peak on the grid; still turning-around near axis).
- D_now=8–12 Gpc: kill a2 = +0.111 → +0.113; peak |dE/E| = 28% → 16.6%.
- D_now=14.7 Gpc: kill a2 = +0.1144; peak |dE/E| = 14.5%.
- D_now=22.28 Gpc: kill a2 = +0.1151; peak |dE/E| = 11.9%.
- D_now=49.17 Gpc (headline-like): kill a2 = +0.1149; peak |dE/E| = 9.69%; static a2 = +0.0182.
- D_now=70–94 Gpc: kill a2 ≈ +0.115; peak |dE/E| ≈ 9.2–8.9%.
- Plateau: for D_now ≳ 8 Gpc, kill |a2| saturates near 0.11–0.115 — about 10,000× the CMB limit 1e−5. Best point on the physical track is near shell-min (kill a2≈0.022), still ≈2,200× over.
- Compared with fixed R=14.3 (cycle4/dtt): headline kill a2 drops from ≈0.385 to ≈0.115 because R_see is ≈3.4× smaller, but the verdict does not change.

## Does any physical D_rec pass CMB (~1e−5) while keeping cluster fit?
- No. Every track point keeps cluster a2=0.084 by construction, and every point has |kill a2| ≥ 0.022 ≫ 1e−5.
- Static |a2| alone is 4.5e−4 to 2.1e−2 — also above 1e−5 except it never cancels the Doppler term on this coasting track (Meff flips sign near D≈7 Gpc; Doppler does not).
- R_see≈c·t_now is fixed by age + T0/3000; shrinking R_see further would require a younger universe or later recombination (larger f), which is outside this continuous T∝1/a track with today's T0=2.725 K.

## Files
- `compute_track.py` → `track.json`, `track.log`
- `plot_track.py` → `kill_a2_peak_vs_Drec.png`, `kill_a2_vs_Dnow_Rsee.png`
- `plot_kill_0_360.py` → `kill_term_0_360_track.png` (D_now = 4.5, 14.7, 22.28, 49.17)
