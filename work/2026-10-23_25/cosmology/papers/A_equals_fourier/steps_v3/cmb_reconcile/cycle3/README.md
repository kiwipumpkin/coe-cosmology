# Cycle 3: redshift-split Migkas test and far-core two-mode mix (exact GR)

## Files
- `data/migkas2020_tablec1.tsv`: per-cluster data for Migkas et al. 2020, from VizieR J/A+A/636/A15 (313 clusters: z, l, b, T, L_X). I could not find a Migkas 2021 (A&A 649 A151) per-cluster table on VizieR, so the 2021 numbers used here are the published ones.
- `migkas_A.py`: rebuilds the Fig. 6 A(l)/A_all profile (Δl = 90° bins, all latitudes, 5° steps). On the full sample it correlates with the digitised Fig. 6 at 0.956 (rms 0.022).
- `slices.py` → `slices.json`: first pass. The L–T slope was free in every bin, and the z<0.03 bootstrap blew up.
- `slices2.py` → `slices2.json`: robust pass. It adds:
  - the L–T slope fixed per slice;
  - 16–84 % bootstrap errors;
  - a 500-shuffle isotropy test;
  - the 2ℓ component aligned with the full-sample phase;
  - common-orientation model tests.
- `validity.py`: regularity check over the full domain, r = 0.001–60 Gpc: 1+2E > 0 and no shell crossing. The old `rm_core.validity` only scanned the shell interval.
- `cmbfix.py`: exact-LTB CMB with event-stopped rays, emitted at a fixed local age since the local bang. The per-ray null-condition error is reported.
- `mix3.py` → `mix3_{edge,lin,null,null2,null3,null4}.json`: two-mode scan. `bump_validity*.json`, `bump_edge.json` and `validity_ladder.json` hold the regularity scans.
- `hotspot.py` → `hotspot.json`: 96 Gauss rays plus 96 rays over ψ = 0.1–60°.
- `make_fig.py` → `cycle3_summary.png`.

## Part 1: redshift split (real per-cluster data, Migkas 2020)
- **Slices:**
  - z<0.03: n = 32, median z 0.018.
  - 0.03–0.07: n = 113, median z 0.049.
  - z>0.07: n = 168, median z 0.113.
- **2ℓ amplitude of A(l), B fixed per slice. p_iso is the chance that an isotropic sky gives a 2ℓ at least this large.**
  - z<0.03: 0.12 ± 0.08, p_iso 0.54 (isotropic 95th percentile 0.29). Not detected; this slice can't see anything below about 0.3. It has a large 1ℓ (0.40 ± 0.16, p 0.07).
  - 0.03–0.07: 0.15 ± 0.04, p_iso 0.012. Its phase is 154°, against 8° for the full sample, so only 0.055 of it lines up with the full-sample pattern.
  - z>0.07: 0.064 ± 0.031, p_iso 0.16. Its phase is 15°, aligned with the full sample (aligned part 0.062).
  - All: 0.084 ± 0.029, p_iso 0.03, phase 8°.
- **Answer to the key question:**
  - The 2ℓ pattern is not concentrated at z<0.03. It is undetectable there, and the full-sample phase is carried by the 0.03–0.07 and z>0.07 slices.
  - The evidence is weak, about 2σ per slice. This is suggestive, not decisive.
- **Near-core variant (exact LTB, mcc 0.08), no-sine 2ℓ predicted at the slice depths:**
  - D = 0.2 Gpc: 0.088 / 0.059 / 0.0017.
  - D = 0.3 Gpc: 0.088 / 0.083 / 0.009.
  - D = 0.5 Gpc: 0.088 / 0.087 / 0.058.
  - Against the z>0.07 data (0.064 ± 0.031), D = 0.2–0.3 sits about 1.8–2σ low. D = 0.5 is consistent, but D = 0.5 fails the CMB (cycle 2: a2 ≈ 1e-4).
- **R² per slice, with l0 and the sine refit in each slice (Chris's convention: model curve vs binned data), for z<0.03 / 0.03–0.07 / z>0.07:**
  - Far core (headline, sampled at the slice depth): 0.50 / 0.81 / 0.72.
  - Near core D = 0.2: 0.54 / −0.24 / 0.18.
  - Near core D = 0.3: 0.53 / 0.83 / 0.32.
  - Near core D = 0.5: 0.52 / 0.87 / −1.00.
  - Dipole only: 0.42 / 0.44 / 0.14.
- **R² per slice with one orientation for all depths (l0 and sine taken from the full-sample fit, as a single physical CoE direction requires):**
  - Far core: −0.16 / 0.11 / 0.01. On the full sample it is 0.86.
  - Near core (all three D): about −0.2 / −2.0 to 0.16 / −1.0 to −2.4.
  - No model reproduces the slices with a fixed orientation. The full-sample R² is not stable under a redshift split.
- **Where the near-core curves spike:** lines of sight that pass the near core (D ≤ z) give a sharp spike, visible in the figure. The D = 0.2 failure at z>0.07 is not caused by the spike; its 2ℓ amplitude there is only 0.002.
- **Published cross-checks:**
  - Migkas 2021 find a 9 % H0 dipole toward (l, b) ≈ (280°, −15°), or equivalently a bulk flow of about 900 km/s extending to at least 500 Mpc. The direction is stable over z ≈ 0.05–0.15.
  - CosmicFlows-4 (Watkins et al. 2023): 387 ± 28 km/s at 150 Mpc/h and 419 ± 36 km/s at 200 Mpc/h, toward l ≈ 297°, b ≈ −6°. The ΛCDM expectation is 120–139 km/s.
  - Watkins, Clark & Feldman 2026 attribute the CF4 flow mainly to mass beyond 200 Mpc/h.
  - All of these point to a signal that reaches beyond 0.3 Gpc. That does not favour a D ≤ 0.3 Gpc near core that produces everything.

## Part 2: far-core two-mode mix in exact GR
- **New finding: the single-bang far-core headline is not a valid spacetime.**
  - With the cavity filled (core only), 1+2E < 0 over r ≈ 0.2–18 Gpc around the CoE (minimum −24). The core's GM/c² is about 81 Gpc, larger than the 33 Gpc to Earth.
  - Single bang is regular only for core mass ≤ 0.005 × headline (Mc/Mcrit ≤ 4e-4), where the A(ℓ) R² is 0.28.
  - Other masses (× headline): 0.01 gives R² 0.29, irregular; 0.1 gives 0.44, irregular; 0.3 gives 0.67, irregular; 1 gives 0.88, irregular.
- **Correction to cycle 2.**
  - The headline a2 = 0.61 and the CoE hot spot were computed through the singular region. 25 of 96 rays (ψ = 1.4–46°) have null-condition error of order 1. The null check sampled only every 8th ray and missed them.
  - The cycle-2 linear law (a2 ≈ 34·Mc/Mcrit, masses ≤ 0.1 × omb) and all near-core models (D = 0.2–1.6 Gpc) are regular and stand.
- **The only regular high-mass models have a young core.**
  - This means a bang-time bump t_B(r) of about 0.999–0.99999 × the age, with width 11–14.5 Gpc. The core region banged within the last 0.1 % of cosmic time.
  - Best regular R²: 0.82 (core × 0.7, width 12 Gpc; core × 1.0, width 13.5 Gpc). No regular model reaches R² 0.85.
- **CMB in those regular models (96 rays, fixed-age emission):**
  - a2 ≈ 1.0, a3 ≈ −0.07, a4 ≈ 0.64 to 0.72. With 24 rays the values are 0.64–0.76, so they are not converged.
  - Nulling the dipole needs an observer speed of 0.99c.
  - The hot spot survives and is enormous: ΔT/T ≈ 500 at 0.1°, 50 at 1°, 9 at 5°, 4 at 10°, 1.5 at 20°, 0.4 at 40°. Half-maximum width is about 0.7°.
  - The failure is about 1e5 above the 1e-5 bound.
- **Two-mode nulling works only in the small-core regime (core × 0.005, R² 0.28), using an older-core bump (negative bang-time amplitude).**
  - With 48 rays: best max(|a2|, |a3|, |a4|) is 7e-5 (width 21 Gpc, amplitude −3.78e-5 of the age), from 1.25e-2 unmixed.
  - With 96 rays, re-tuned: 1.3e-4 (width 20 Gpc, amplitude −4.7e-5).
  - Ray-count convergence of the unmixed a4: 1.15 / 1.25 / 1.38 / 1.35 e-2 at N = 24 / 48 / 96 / 144, so the residual carries about 1e-4 of resolution noise.
  - Two knobs cannot null three multipoles. At the optimum, a3 and a4 are left with opposite signs.
- **Tuning precision:**
  - The slopes are da_ℓ/d(amplitude/age) ≈ 200–330.
  - Holding |a_ℓ| < 1e-5 would need the amplitude fixed to about 4e-8 of the age, which is 0.08 % of the amplitude itself.
  - Even with that tuning, the result misses 1e-5 by about 10×.
- **Positive (younger-core) small bumps only increase a2–a4.**
  - An older-core bump larger than the recombination age (about 9e-5 of the age) leaves rays toward the CoE with no last-scattering surface.

## Figure (`cycle3_summary.png`, viewed)
- **Top row:** A(l) per slice (the rebuilt Migkas profile ± bootstrap), with the far core and near core D = 0.2 / 0.5 fits.
- **Middle row, left to right:**
  - 2ℓ amplitude vs slice depth, with the isotropic 95 % line and the model predictions.
  - R² vs core mass, marking regular vs irregular models.
  - The mode-mix residual.
  - Regular models placed on the residual–R² plane.
- **Bottom row:** the CoE hot-spot profile, and whole-sky temperature for the regular young-core models.

## Caveats
- Data are the Migkas 2020 L_X–T sample only (313 clusters), not the 2021 ten-relation sample. The profile rebuild uses a simplified fit (see `migkas_A.py`).
- Slices are small: z<0.03 has only 32 clusters, and every per-slice 2ℓ is 2–2.5σ at best.
- R² compares a model point curve with binned data (Chris's convention). Binning the model changes values by about 0.1. Refitting l0 and the sine per slice inflates R².
- The LTB code uses an internal H = 70 unit convention; only ratios are used.
- The decaying mode is represented by the bang-time function, the exact LTB decaying mode.
- Emission at fixed local age is a proxy for recombination. Fixed local density is undefined on rays toward a young core, because the density there never reaches the recombination value before the local bang.
- Hot-spot multipoles are not converged even at 96 rays, but the conclusion (order-1 residual) is resolution-proof.
