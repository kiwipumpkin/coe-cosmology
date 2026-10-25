# CMB reconciliation, cycle 2

Scripts (Chris's integrators and LTB code are used unchanged):
- `A_ltb_cmb.py`: exact LTB CMB for the far core and the mass ladder. Writes A_ltb_cmb.json.
- `A_full.py`: full or scaled mass. Rays are traced only to z ≈ 1090 so they don't run into the core singularity. Writes A_full_omb_core_x0.5.json.
- `B_hybrid_fit.py`: the geometry + tilt ladder, and the core-only (no-cavity) refits. Writes B_hybrid_fit.json and B_hybrid_fit_nocav.json.
- `C_near_core.py`, `C_small.py`, `C_depth.py`: the near-core variant in exact LTB, covering the A fit and the CMB. Write C_near_core.json, C_small.json and C_depth.json.
- `make_fig.py`: writes cycle2_summary.png.
- Shared helpers are in `ltb_common.py`.

Exact GR setup: LTB dust with Λ = 0, plus a point-like core (Lc = 0.005 D), Ω_s = 0.049.
- Single bang time, so this is the pure growing mode, the same construction as Chris's relativistic fit.
- E(r) comes from the areal label today, and t0 from a local mean H = 70 at Earth. This is the code's unit convention; only dimensionless ratios are used.
- The dust is extended outward so that the last-scattering surface sits in matter.
- CMB = ΔT/T of 1/(1+z). Emission is at fixed proper time with the z monopole set to 1089. The dipole is removed with an axial boost. a_ℓ are Legendre coefficients about the CoE axis.

## A. Far core, nonlinear GR
- **Does the potential stay constant? Yes, in the growing mode.**
  - The ℓ=2 of the emission-epoch redshift is flat to 5–10% from z ≈ 1000 down to z ≈ 40.
  - It changes only after z ≈ 20, once curvature dominates.
  - Linear growth gives Φ(emission)/Φ(today) = 1/0.104 = 9.6, not 1100.
  - The 1/a deepening from cycle 1 belongs to a constant-mass core on a fixed homologous flow, which is a decaying mode.
- **The self-similar "Mc ∝ D ∝ a" picture.**
  - Holding Mc/Mcrit(D) fixed in a coasting flow does mean Mc ∝ t, so GMc/D stays constant.
  - Dust in GR conserves the mass inside every comoving shell, so real core growth needs an energy source (flagged extension).
  - The growing mode achieves the same thing without one. At early times the dust is arranged so the density is uniform, and the effective perturbation grows with the structure.
  - ISW does not vanish, though: in an open Ω = 0.049 universe Φ decays by about 10× after z ≈ 20.
- **Exact numbers, far core, Ω_s = 0.049 geometry (D = 49 Gpc).** a2 ≈ 34 × Mc/Mcrit in the small-mass regime.
  - 5.8e-5 at Mc/Mcrit = 1.7e-6; 6.3e-4 at 1.7e-5; 3.9e-3 at 1.7e-4.
  - That is about 6× the static piece (5.6 × Mc/Mcrit) and about 1/90 of cycle 1's constant-mass linear estimate (about 3160 × Mc/Mcrit).
  - Emitting at fixed local density gives identical numbers, so the self-regulating front does nothing. This is now exact.
  - Null test: Mc/Mcrit = 1.7e-12 gives 5.8e-11 (null-condition error 1e-13).
- **Honest far-core number: a2 = +0.61** (a3 = −0.33, a4 = +0.12) for the one self-consistent exact-GR headline, the relativistic fit.
  - That fit has Mc/Mcrit = 0.083, R² = 0.885, H_r/H_⊥ = 1.20.
  - 0.61 is about 6e4× the 1e-5 bound. The static linear piece alone is 0.47.
- **Why the Newtonian masses can't be carried over.**
  - In the growing mode the A pattern needs Mc/Mcrit ≈ 0.06–0.08.
  - The Newtonian fits' 0.0017–0.0027 is bought by the τ_i = 0.01 homologous start. That is a decaying mode, which is exactly where Φ does grow into the past as 1/a.
  - Either way the CMB fails. The growing mode is the less bad route.
- **Hot spot at the CoE (new, exact GR).**
  - With Λ = 0 and baryons only, the last-scattering surface lies about 19 Gpc away along the ray and at about 180 Gpc areal radius today. The core (areal D = 33–94) sits only 12–17 Gpc away along the ray, so it is inside the last-scattering sphere.
  - Photons from behind the core cross its field and arrive blueshifted. At half the Ω_s = 0.049 Newtonian mass, z = 409 instead of 1097 at 5.6° from the CoE, and 935 at 12.9°. The resulting a2 = 0.078, a3 = 0.10, a4 = 0.12 are dominated by this spot.
  - Outside a 25° cap the spread is 8.6e-3.
  - At the full Newtonian mass, rays within about 13° of the CoE stall in the core's strong field.
- **Cases that can't be traced.**
  - The as-fitted cavity cases: rays within about 40° of the CoE enter the vacuum strong-field region and the ray solver fails.
  - The Ω_s-free fit, core only: filling the cavity at Ω_s = 1.01 gives 1+2E < 0, i.e. a region closed on a scale below D. It cannot be embedded in LTB.

## B. Geometry + tilt
- The tilt enters as the multiplicative dipole factor, the same form as the sine term; a bulk flow changes the apparent rate as v·cos(angle)/(cz). Refit with Mc/Mcrit fixed; D, kin, kout, ℓ0 and the tilt are free; Ω_s = 0.049 rule.
- R² against the fixed Mc/Mcrit: 0.893 at 1.7e-3, 0.800 at 1e-3, 0.601 at 5e-4, 0.431 at 2e-4, 0.361 at 1e-4, 0.333 at 2e-5, 0.343 at 2e-6.
- Core only (no cavity; uniform dust + core) works just as well for A: Mc/Mcrit is free at 1.65e-3 with R² = 0.893. The no-cavity ladder gives 0.892, 0.601, 0.361, 0.270.
- **CMB from the remaining core.** Static piece 5.6 × Mc/Mcrit; exact growing-mode slope about 34 × Mc/Mcrit.
  - Headline: 9.6e-3 (static) or about 6e-2 (exact slope, an estimate).
  - At 1e-4: 5.6e-4 or 3.4e-3.
  - Passing needs Mc/Mcrit ≲ 3e-7 (exact) or 2e-6 (static), where R² = 0.27–0.34.
- **Verdict: no.** The tilt carries only the dipole. The 2ℓ term (amplitude 0.073) is what sets the core, so lowering the core need does not bring a2 toward 1e-5 while the fit survives.
- **Tilt size.** The dipole amplitude is 0.052–0.062.
  - At Chris's 14.3 Gpc sampling that is about 0.2c, i.e. not a physical bulk flow.
  - At the Migkas depth of 0.5 Gpc it is about 2000 km/s.

## C. Variant (not a replacement): near core, D ≈ 0.2–2.5 Gpc, sampled at Migkas depths
- **Newtonian.** Chris's integrator is scale-free. With the shell rule written in units of the sampling radius, the fit repeats exactly (R² = 0.893 at D ≈ 3.4× the sampling radius) for any scale.
- **Exact LTB, growing mode** (the mass grows with the structure; core + uniform dust, no cavity), A sampled at z = 0.12 (about 0.5 Gpc).
  - Best grid point: Mc/Mcrit = 0.08, D = 1.6 Gpc, with R² = 0.884 (0.870 at z = 0.075, 0.656 at z = 0.3) and H_r/H_⊥ = 1.195.
  - Traced CMB (SW + ISW + Rees–Sciama, all included in the exact ray trace): a2 = 3.2e-3, a3 = 2.0e-3, a4 = 1.1e-3, i.e. 320× over.
  - D = 1.0 Gpc: R² = 0.834, a2 = 7.7e-4.
- **Depth law.** At fixed Mc/Mcrit, a2 ∝ D³, approximately a2 ≈ 3.2e-3 × (D/1.6 Gpc)³. Shrinking D together with the sampling depth:
  - D = 0.5 Gpc (sampled at z = 0.04): R² 0.875, a2 9.5e-5.
  - D = 0.3 Gpc (z = 0.025): R² 0.871, a2 2.2e-5 (1.4e-5 at Mc/Mcrit = 0.05).
  - **D = 0.2 Gpc (z = 0.016): R² 0.871, a2 6.5e-6, a3 5.1e-6, a4 3.6e-6. This passes.**
- **The catch.** The passing variants fit only when sampled within about 100 Mpc.
  - D = 0.2: R² = 0.67, −1.7, 0.47, 0.31 at z = 0.03, 0.05, 0.075, 0.12.
  - D = 0.3: R² = 0.57, 0.71, 0.43 at z = 0.05, 0.075, 0.12.
  - Migkas report the anisotropy extending to at least 500 Mpc.

## D. New ideas, with quick numbers
- **Mode trade-off.**
  - Growing mode: Φ is flat (Chris's intuition is right), but the core must be about 50× bigger, giving a2 ≈ 0.6.
  - Decaying mode (the Newtonian fits): small core, but Φ ∝ 1/a, giving a linear a2 of about 5–8.
  - No mix of pure far-core modes gets near 1e-5 without exact tuning (the cycle 1 candidate 2 tuning, 1e-3 to 1e-6).
- **CoE hot spot as a clean falsifier.** Any far core closer than the last-scattering surface (areal D ≲ 180 Gpc for Ω = 0.049, H = 70) gives a hot spot of order unity within about 15° of the CoE.
  - The relativistic fit's R² stays flat for D from 33 Gpc to above 1e6 Gpc, so D > 200 Gpc is allowed and removes the spot.
  - It does not remove the smooth a2 ≈ 0.5, which does not depend on D.
- **Depth law as the only open door.** Keeping Chris's geometry (core + V = kr + off-centre observer) and passing the CMB needs D ≲ 0.2–0.3 Gpc. The 2ℓ Migkas pattern would then have to live in clusters at z ≲ 0.02–0.03.

## Ranked status
- **1. Near-core variant at D ≈ 0.2 Gpc.** Passes: a2 6.5e-6, R² 0.87 when sampled at z ≈ 0.016. Fails at Migkas depth: R² ≤ 0.47 at z ≥ 0.075. Same geometry shape, very different scale.
- **2. Near-core variant at D ≈ 0.3–0.5 Gpc.** Marginal: a2 1.4e-5 to 9.5e-5, R² ≈ 0.87 at z = 0.025–0.04, and 0.71 at z = 0.075 for D = 0.3.
- **3. Far core, any form: dead.**
  - Exact growing mode: a2 = 0.61 (relativistic fit).
  - Newtonian mass in GR: a hot-spot-dominated a2 ≈ 0.08 at half mass, with the fit not reproduced in the growing mode.
  - The tilt does not lower the core need.

## Proposed cycle 3
- **Test the depth law against the data.** Split the Migkas sample by redshift (z < 0.03, 0.03–0.07, > 0.07) and fit the near-core variant to each slice.
  - The passing variant predicts the 2ℓ pattern concentrated at z < 0.03 and fading as the sampling sphere approaches the CoE.
  - Use the actual cluster redshift distribution, not a single sampling sphere.
- **Near-core variant, properly optimised.**
  - Continuous LTB fit over D, Mc/Mcrit and shell edges (Chris's cavity rules scaled to the sampling radius).
  - Full-sky CMB with N = 96 rays and ℓ ≤ 10.
  - Check the predicted local density contrast and bulk flow against CosmicFlows-4 and 2M++.
- **If Chris keeps the far core.** An exact-GR bang-time (decaying-mode) plus growing-mode mix in LTB, to see whether any smooth mode mixture nulls ℓ = 2–4 at once. This is the exact version of cycle 1's candidate 2.
  - Also a resolved map of the CoE hot spot (dense rays near the core) for D < 180 Gpc.

## Caveats
- The LTB runs use the code's H = 70 unit convention, which is 2.5% below the model ⟨H⟩ of about 71.6; ratios are unaffected at this level.
- 24 Gauss rays do not resolve the hot spot, so a_ℓ for the hot-spot cases is only indicative.
- B's CMB values use the exact-GR slope from the Ω_s = 0.049 geometry, not a re-trace at each D.
- C's fits are grid points, not continuous optima, and use core + uniform dust without Chris's cavity shell rules.
- The as-fitted cavity cases and the Ω_s-free fit could not be traced in GR (strong field, and 1+2E < 0 respectively).
- The data are the digitised Migkas 2020 Fig. 6 profile. Migkas 2021 describe the effect as "nearly dipole" and extending to at least 500 Mpc; that is the claim the near-core variant must confront.
