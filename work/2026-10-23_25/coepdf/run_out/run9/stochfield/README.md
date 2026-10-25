# Stage A: stochastic density field → ensemble temperature C_ℓ, regularized inverse fit (CoE flow, R = 4 Gpc)

## Model (5 parameters)
- Gaussian random field, net-zero mean. P_δ(k) = A·k^n·exp(−(k/k_c)²) [Gpc³]. Φ_k = −4πGρ̄ δ_k/k² (ρ̄ = rho_m from cosmo.py).
- Radial envelope on the potential: Φ = W(χ)·φ_g, with W = Gaussian(centre χ_w, width σ_w) in distance from the observer. Parameters: lgA, n, k_c, χ_w, σ_w.
- The field rides the CoE flow. The kernels are copied from los_flow_R/flowR.py (grow variant, A(t) = D/a, path rescaled so χ/R = η/η_rec, R = 4 Gpc):
  SW at the source, the ISW-like term (2/c²)ΔA·Φ, and the moving lens using the **full CoE radial velocity** (mlf). That velocity is superluminal (Ṙ ≈ 4c today, far larger early on; HEURISTICS #14), and that is why a density contrast as small as σ_δ ~ 1e-4 is enough. The moving lens is ≥ 99.99% of the field's power (scalar SW+ISW / total ≤ 9e-5).
- Lensing is left out (HEURISTICS #3: 1e-3–1e-2 relative at low ℓ). Ṙ is evaluated at x = 1. The Lagrangian displacement and the per-direction variation in the CoE unit vector are ignored in the analytic model but included in the MC check.
- The deterministic CoE background (smooth red + RS_cc map) is either added to the ensemble mean ("bg") or left out ("nobg").

## Forward model (fwd.py): exact spherical-Bessel projection, no Monte Carlo
C_ℓ = (2/π)∫k²dk P_Φ[ |∫K_s W j_ℓ|² + moving-lens term ]. The moving-lens term is u·∇(Wφ) = W u·∇φ + φW′(n·u). It is evaluated exactly in its m-averaged form using cosθ Gaunt couplings (ℓ±1).
**Validation** (mcval.py, validation_mc.png, mc_validation_chosen.json): 8 Gaussian realizations on the fine grid (0.078 Gpc), passed through the actual los_flow_R `comp()` at R = 4.
MC/analytic ratio: mean 0.97 over ℓ 2–30 (range 0.75–1.31, consistent with 8-realization scatter), 0.82 over ℓ 31–60, 0.58 over ℓ 61–100. Above ℓ ~ 50 the drop is grid smoothing and interpolation, so the analytic curve is not MC-validated there.

## Inverse (fit.py)
Objective: χ²_CV + λ·Reg, with χ²_CV = Σ(D−D_P)²/(σ_P,side² + (2/(2ℓ+1))·D_field²), using the Planck 2018 TT full file with asymmetric errors.
Reg = σ_δ² (L2 amplitude, 0.078 Gpc smoothing) + power-weighted ⟨(d²lnP/dlnk²)²⟩ + R²⟨W′²⟩/⟨W²⟩ (radial roughness).
Nelder–Mead from multiple starts, with an L-curve over λ = 0 and 1e-2…1e3. I also report χ²(Planck errors only), which uses the same definition as ΛCDM's 20.1, and Pearson r². These are labels only.

| fit | λ chosen | lgA | n | k_c | χ_w, σ_w (Gpc) | σ_δ | χ²_CV | χ²_Planck | Pearson r² |
|---|---|---|---|---|---|---|---|---|---|
| ℓ2–30 bg | 100 | −9.26 | −0.19 | at bound (200) | flat (σ_w ≈ 20) | 1.3e-4 | 10.8 | 17.6 | 0.032 |
| ℓ2–30 nobg | 3.16 | −9.19 | −0.24 | at bound | flat | — | 10.2 | 17.3 | 0.044 |
| ℓ2–30 nobg, λ=0 | 0 | −9.70 | 0.85 | — | 4.02, 0.079 (thin shell at source) | 3e-4 | 9.2 | 18.3 | 0.053 |
| ℓ2–100 bg | 100 | −9.79 | 0.48 | at bound | flat | — | 76.0 | 114.5 | 0.721 |
| ℓ2–100 nobg | 1 | −9.78 | 0.47 | at bound | 2.48, 5.5 | — | 78.1 | 117.0 | 0.725 |
| ΛCDM ℓ2–30 | | | | | | | 10.3 | 20.1 | 0.012 |
| ΛCDM ℓ2–100 | | | | | | | 57.2 | 90.7 | 0.794 |

(Exact σ_δ values for each λ are in fit_*.json and chosen_*.json.)
- At ℓ ≤ 30 the data constrain only the amplitude and the slope n. k_c and the radial window are unconstrained, and the L-curve is flat: regularizing to a uniform window costs Δχ² < 0.1. The simplest consistent statistics are a near-white-noise density field (n ≈ −0.2, no cutoff) filling the whole path.
- ℓ ≤ 100 pulls the slope to n ≈ +0.48.
- In ℓ2–100 nobg, the λ = 10 and 100 entries are optimizer failures (λ = 100 collapsed). Ignore them.
- Degeneracy (d_degeneracy_l30_bg.png): with a thin shell (σ_w = 0.3), shell centres from 0.8 to 3.8 Gpc all reach χ²_CV ≈ 13–18 if the amplitude is rescaled. The required lgA falls from about −7.0 to −9.0 as the shell moves outward, because the moving-lens kernel rises toward the source. A tight n–lgA ridge also appears.
- ℓ ~ 220: the analytic code has no grid limit (an ℓ ≈ kχ structure needs k ≈ 55/Gpc, i.e. a ~0.11 Gpc wavelength feature at χ ≈ 4). MC validation would need grid dx ≲ 0.02 Gpc (k_Nyq ≳ 150/Gpc), a nested grid or ~1000³ cells. A peak would also need a feature in P(k), which this smooth parametric form cannot make. Not attempted.

## Axial modes (axial.py, f_axial_m0.png, axial.json)
The a_ℓm are rotated so that the CoE (Galactic l=180°, b=−15°) is the z axis; a self-test with an axisymmetric blob gives m=0 fraction = 1. The quantity is the excess (2ℓ+1)|a_ℓ0|²/Σ_m|a_ℓm|², where an isotropic sky has 1 and the 95% level is ≈ 3.84.
- Planck SMICA (unmasked, nside 256), ℓ 2–10: 1.32, 1.30, 1.17, 0.00, 0.32, 0.01, 1.34, 1.65, 0.15. Over ℓ 2–30 the only excess above 3 is at ℓ = 29. **ℓ = 5 has essentially zero m=0 power about the CoE axis.**
- The CoE background is purely axisymmetric (excess = 2ℓ+1 at every ℓ, though its power is ≈ 0 above ℓ = 3).
- One field realization: excess > 3 at ℓ = 9, 22, 24, 28. Field+background adds ℓ = 3.
- **No 1:2:3 harmonic ladder (±15%) was found in any sky**, so no axial box length is implied. Had one existed, L = 2πR·j/ℓ_j (with R = 4 or 4.23 Gpc) is the conversion coded in axial.json.

Files: fwd.py, fit.py, figs.py, deg.py, mcval.py, axial.py; fit_l{30,100}_{bg,nobg}.json, chosen_*.json; figures a_cl_*, b_pk_w_*, c_lcurve_*, d_degeneracy_l30_bg, e_skymap, f_axial_m0, validation_mc.
Note: /workspace/joint/lcdm.txt is flat beyond ℓ = 30, so ΛCDM for ℓ > 30 uses the Planck minimum-theory file.

## Acoustic-peak inverse stage (acoustic_model.py, acoustic_fit.py, acoustic_resume.py, acoustic_figs.py → acoustic_fit.json)
Density spectrum: P_δ = P0·k^n·[1 + A·cos(k·r_s + φ)]·exp(−(k/k_D)²). The window W is a Gaussian (centre, width). The physics terms, R = 4 Gpc, and the fixed background (ℓ 2–30) are the same as in stage A.
Parameter count: 8 for free W, 6 for the fixed thin shell at the source (centre 4.0, width 0.1 Gpc). No odd/even asymmetry term was fitted.

**Forward model.** Exact projection for ℓ ≤ 50, Limber above that (vectorized, every ℓ to 2500, about 0.1 s per call).
- Limber vs exact (acoustic_limber_check*.json): for interior or broad shells it agrees within 3% for ℓ ≥ 150. A broad shell (width 0.5 Gpc) agrees within 0.3% even at ℓ = 20.
- For shells cut off at the path end (the source surface), Limber overestimates by about 1.5–5× below ℓ 150 and by up to ±20% around ℓ 300.
- So the best fits use an exact-to-Limber ratio computed at 13 ℓ nodes between 51 and 400 and interpolated. I iterated this twice and refit each time. It ends at ℓ 400, where it leaves a small step in the curve (the ratio there is about 0.92).

**Fitting.** The fit measure is χ² with Planck errors: ℓ 2–29 unbinned plus the binned file for ℓ ≥ 30 (111 points). On this measure ΛCDM scores 98.7.
- Multi-start: r_s ∈ {0.042, 0.063, 0.084, 0.105} Gpc × φ ∈ {0, π/2, π, 3π/2}, with 2 window seeds for free W. Nelder–Mead from each start.
- L-curve: λ ∈ {0, 10, 100, 1000}, with Reg = A² + (0.1/width)² + σ_δ². The best entry is λ = 0 for both fits. Regularizing only raises χ².
- Checkpoints go to acoustic_progress.json.

| | free W (8 params) | fixed thin shell at source (6) | ΛCDM |
|---|---|---|---|
| r_s | 0.0819 Gpc = 81.9 Mpc | 0.0814 Gpc = 81.4 Mpc | |
| φ, A, n | 1.22, 0.19, −0.33 | 1.46, 0.18, −0.02 | |
| k_D | 345 /Gpc | 308 /Gpc | |
| W centre, width | 3.944, 0.050 Gpc (at the lower bound) | 4.0, 0.1 (fixed) | |
| σ_δ (0.078 Gpc smoothing) | 4.2e-4 | 2.5e-4 | |
| fit-measure χ² (Planck err / +CV) | 11453 / 7861 | 17103 / 12299 | 98.7 |
| ℓ 2–30: χ²(Planck err), χ²(+CV), Pearson r² | 142.3, 50.1, 0.010 | 1575, 47.2, 0.152 | 20.1, 10.3, 0.012 |
| ℓ 2–2500 (unbinned) | 11772, 8704, 0.871 | 16380, 12187, 0.588 | 2558, 2027, 0.983 |
| ℓ 30–2500 | 11634, 8656, 0.874 | 14805, 12140, 0.746 | 2538, 2017, 0.983 |
| first 5 peaks | 227, 530, 831, 1131, 1431 | 233, 529, 829, 1129, 1428 | 220, 536, 813, 1127, 1421 (Planck binned spline: 224, 535, 821, 1122, 1412) |

Reading the results:
- A single cosine with r_s ≈ 82 Mpc puts the peaks at an even spacing of about 300 in ℓ, within about 10 of Planck's peak positions.
- The first peak's height is badly missed: the model reaches about 3,600 μK² against Planck's 5,750. Above ℓ ≈ 30 the model also sits above Planck, and that excess dominates the χ².
- The smooth cosine can't give the odd/even peak-height pattern. The free-W fit pushes the window to its 0.05 Gpc minimum at the source surface.
- Plausibility: the density contrast is only about 1e-4 because the moving-lens term is boosted by the superluminal full-velocity frame (HEURISTICS #14).

Figures: acoustic_cl_{free_W,fixed_shell}.png (top: ℓ 2–30 linear; bottom: log-ℓ to 2500, with a ±1σ cosmic-variance band), acoustic_pk_w_*.png, acoustic_resid_*.png, acoustic_lcurve_*.png.

## Single-realization figures (real_figs.py)
Each figure shows one random sky (seed 42) as the main line. Each ℓ is drawn as D_bg + D_field·χ²_{2ℓ+1}/(2ℓ+1), so the scatter is cosmic variance. The deterministic CoE background has no scatter. The ensemble mean stays on the plot as a faint line.

χ²(Planck err) and Pearson r² for this one sky:
- acoustic free W: ℓ 2–30 → 140.1, 0.000; ℓ 2–2500 → 12426, 0.866
- acoustic fixed shell: ℓ 2–30 → 1508, 0.168; ℓ 2–2500 → 16967, 0.595
- stage A: l30_bg 37.4, 0.000; l30_nobg 37.6, 0.029; l100_bg 177.3, 0.607; l100_nobg 179.5, 0.616

A different seed gives different numbers.

Files: acoustic_cl_free_W_real.png, acoustic_cl_fixed_shell_real.png, a_cl_{l30,l100}_{bg,nobg}_real.png

## Acoustic ladder (acoustic_ladder.py → acoustic_ladder.png, checkpoint acoustic_ladder.json)
Setup: free_W best-fit P(k), a shell at R−0.056 Gpc with width 0.05 Gpc, and multipliers ×0.01–×100 (9 log steps). Each level is drawn as one sky (seed 42+j, χ² scatter) with a faint mean line. Field only, with no background.

**Caveat:** there is no exact-node correction, so Limber overestimates by about 1.5–5× for ℓ 50–150. That shows as a visible jump at ℓ = 50, worst at R = 8 and 14.

Planck first-peak height ≈ 5,750 μK² (binned).

| R (Gpc) | first peak height at ×1 | levels bracketing Planck's 1st peak | ×1 multiplier needed | model peaks (first 5) |
|---|---|---|---|---|
| 2 | 6,823 | ×0.316–×1 | 0.85 | 262, 412, 560, 708 |
| 4 | 4,403 | ×1–×3.16 | 1.32 | 218, 530, 831, 1131, 1431 |
| 8 | 1,553 (at ℓ ~ 438) | ×3.16–×10 | 3.73 | 438, 1065, 1672, 2276 |
| 14 | 824 (at ℓ ~ 769) | ×3.16–×10 | 7.03 | 769, 1869 |

Planck: 224, 535, 821, 1122, 1412. The peak positions scale as ℓ ∝ R for fixed r_s. Only R = 4 lines up.

## Corrected ladder at R = 4, 5, 6 Gpc (acoustic_ladder456.py → acoustic_ladder_456.png, acoustic_ladder_456.json, checkpointed per node)
The exact/Limber ratio is computed at 13 nodes (ℓ 51–400) for each R. The ℓ = 50 seam is now gone: D51/D50 = 1.005, 1.002, 1.004.

At ℓ = 400, where the correction stops, a step of about 9–19% remains. It creates spurious "peaks" at ℓ 412 (R = 5) and 419 (R = 6), which are excluded below.

| R (Gpc) | 1st peak ℓ, height at ×1 (μK²) | bracketing levels | multiplier for Planck's 1st peak (5,750) | model peaks |
|---|---|---|---|---|
| 4 | 228, 3,576 | ×1–×3.16 | 1.62 | 227, 530, 831, 1131, 1431 |
| 5 | 282, 2,661 | ×1–×3.16 | 2.18 | 284, 663, 1042, 1417 (four found) |
| 6 | 339, 2,061 | ×1–×3.16 | 2.81 | 341, 797, 1252, 1703 (four found) |

Planck: 224, 535, 821, 1122, 1412.
