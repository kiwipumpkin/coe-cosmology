# CMB reconciliation, cycle 1

Run: `python3 estimates.py` (about 7 s, writes estimates.json), `python3 c5_kinematic_fit.py` (about 6 min, runs Chris's Newtonian integrator unchanged and writes c5_kinematic_fit.json), `python3 make_fig.py` (writes cycle1_summary.png).

**Conventions.** θ is measured from the CoE direction. a_ℓ are Legendre coefficients of ΔT/T, and "pass" means |a_ℓ| ≲ 1e-5 for ℓ ≥ 2. The fits are the dscan headline points (free D=93.7, omb D=49.2, rel/LTB D=32.9 Gpc) with the model ⟨H⟩ = 71.79 / 71.61 / 71.63. R = 14.3 Gpc, emission at a = 1/1101, coasting null-ray shape. The core is the only source; the cavity is excluded from all CMB numbers. The data-driven shear is taken from the Migkas Fourier form: the 2ℓ amplitude 0.073 needs β = dln h/dln q = 0.16 at Earth, i.e. H_r/H_⊥ = 1.16.

## Priority question: how much net photon redshift survives, and how it compares with the matter velocity effect
- **Photons do not cross the core, so nothing cancels on the way in and out.** A CMB photon starts on the R sphere and ends at Earth. Even in a fully static potential the ΔT/T is Φ(emit, θ) − Φ(obs), and its ℓ=2 part is the core's tidal field across the sphere.
- **(1) Static piece, using today's field.** a2 = −1.5e-2 (free), −9.6e-3 (omb), −0.47 (rel). This matches a2 = (Mc/Mcrit)(HR/c)²/2 to 4 digits.
- **(2) Time-varying piece.** With constant Mc the potential goes as 1/a, which gives SW −16.3 / −10.5 / −513 and ISW +24.6 / +15.9 / +776. ISW/SW = −1.51 in every fit, which is a property of the ray geometry, not of the fit. Net: +8.3 / +5.4 / +262.
- **(3) Matter velocity effect.** The velocity pattern the A fit needs, n·Δu/c over the same R sphere, has a2 = −0.37 in every fit, because it equals (2β/3)(HR/c).
  - **Ratio of static photon shift to velocity effect:** 0.040 (free), 0.026 (omb), 1.27 (rel) at R. At the Migkas bulk-flow depth of 0.5 Gpc it is 1.4e-3 / 9e-4 / 0.044. The ratio grows linearly with distance, because the static piece scales as s² and the velocity as s.
  - **Ratio of the net (constant Mc) shift to velocity effect:** −22 / −15 / −710.
- **Why the static ratio is small, and why that doesn't save the model.** In the Newtonian fits the velocity shear is built by the core's strong early pull. The τ_i = 0.01 start amplifies it about 99×, so today's potential is weak compared with the velocity it left behind: the ratio is roughly τ_i × (HR/c). But that same history, with a constant-mass core, makes Φ about 1100× larger at emission. So the photons see 15–700× the velocity effect, not 3%.
- **Can the path be arranged to climb back out?** Only in one configuration, all three together:
  - Tolman equilibrium at emission (T·√(−g00) = const);
  - a potential that stays static from emission to today;
  - Mc growing ∝ a.
  - Tolman equilibrium needs sound to cross the core scale D before emission. For a coasting universe that means ln(a_e/a_i) ≳ 20–39 for one crossing and 200–390 for ten, i.e. a start at a_i ≲ e^-20.
  - An emission sphere around an off-centre observer cannot be an equipotential: Φ varies by about 1e-2 across it (omb).

## Candidates
- **1. Velocity-only anisotropy (Φ ≈ 0; inherited, pressure-driven or initial).**
  - (a) Mechanism: the flow shear is carried by the matter, not by a potential.
  - (b) Equation: for dust with Φ = 0, d(n·u)/dl = n_i n_j ∂_i u_j, so ΔT/T = n·(u_obs − u_emit)/c exactly at linear order. u_obs gives only a dipole, so ℓ ≥ 2 comes purely from the emitters' velocities.
  - No-go: for any flow profile V = k r h(r) centred on the CoE, the ℓ ≥ 2 part vanishes only if h is constant on [D−R, D+R]. That range contains Earth, so β = 0.
  - (c) Numbers: today's pattern at R gives a2 = −0.37, failing by 3.7e4. A frictionless (ballistic) flow has peculiar velocity ∝ 1/a, which makes it −405 at emission; that would need σ/H today ≲ 1.3e-9.
  - (d) Dead, and it doesn't matter whether the origin is pressure, initial conditions, or gravity that later switched off.
  - (e) Prediction: a kinematic CMB quadrupole aligned with the CoE axis, with a3/a2 ≈ 0.23–0.65.
- **2. Exact SW–ISW cancellation through the core's mass history.**
  - (b) Equation: in the tidal limit the net ℓ=2 is proportional to −Φ_e + 2⟨Φ⟩, where the average is along the ray weighted uniformly in χ².
  - Mc ∝ a^0.635 cancels ℓ=2 exactly (checked numerically, all fits). But p must be tuned to 3e-5 (Newtonian) or 1e-6 (LTB), and ℓ=3 is left at 5–6e-3, failing by about 600.
  - An instant halving of Φ just after emission, then a static potential, cancels every ℓ. The halving ratio must be 0.497–0.500. With the drop at 1.02 a_e the residual is still 1e-4. The drop has to happen everywhere on the sphere at once, which is not causal in GR (Birkhoff: outgoing mass keeps its exterior field until the radiation shell passes).
  - Tuning when the V = kr flow starts (a_i) cancels ℓ=2 at a_i = 3.6e-4 but ℓ=3 at 5.7e-4. At the ℓ=2 zero, a3 = −1.1 to −100.
  - (d) Unnatural: needs tuning of 1e-3 to 1e-6 and is acausal.
  - (e) Prediction: an ℓ=3/ℓ=4 residual with a CoE-axis symmetry.
- **3. Self-regulating recombination front (emission at local T = 3000 K).**
  - For an external (Laplace) core field with Φ ∝ 1/a, the Newtonian-gauge δ_γ = 0. The front then sits at the same coordinate time, so there is no timing compensation and the SW factor is 1.
  - The textbook 1/3 applies only to self-sourced adiabatic perturbations. Even granting it, a2 = 3.2e-3 (omb) to 0.16 (rel), failing by 300 or more.
  - Only the Tolman variant from the priority section cancels. It needs a_i ≲ e^-20, a static potential and Mc ∝ a.
  - (d) Dead as stated. The Tolman variant needs tuning of all three at once.
  - (e) Prediction: a temperature–polarisation (TE) signature of hydrostatic stratification at last scattering.
- **4. Late thermalisation, or a late-forming core.**
  - Last scattering would have to be at z < 0.24 (omb), 0.19 (free) or 0.03 (rel), in front of the Migkas clusters themselves. Dead.
  - A core switched on at a_f, with its mass rescaled to keep the fitted shear, gives an ISW of |a2| = 0.28 at a_f = 0.01 and still 3.8e-4 at a_f = 0.99 (omb). Roughly a2 ~ K·ln(1/a_f), with K = the fit's shear parameter. Dead.
- **5. Kinematic centre (Mc → 0).**
  - Without core and cavity, A is constant and only the sine term survives: R² = 0.27. That is the best single dipole on the 71 Migkas points.
  - Chris's integrator with mcc = 2e-6 and a refit: R² = 0.34 (Ω_s = 0.049 rule) and 0.82 (Ω_s = 1, D = 94, ℓ0 = 107°).
  - That 0.82 comes entirely from the shell's empty cavity (M_hole about 230× Mc). Its CMB a2 is about 3.4 even with the potential frozen, so it is excluded.
  - An mcc sweep at the omb headline geometry gives R² = 0.89, 0.76, 0.41, 0.28 and 0.33 at mcc = 1.7e-3, 1e-3, 3e-4, 1e-4 and 2e-6.
  - For reference, a free fit of harmonics up to 2ℓ gives 0.81, and up to 3ℓ gives 0.92.
  - (d) Passes the CMB but loses the fit: R² 0.91 → 0.27–0.34.
- **6. Tilt (matter frame ≠ CMB frame).**
  - A uniform tilt changes A only as a dipole that falls off as 1/z. The CMB cost is a dipole, which is removed anyway. So the tilt can supply the sine or dipole part: v ≈ 790–1600 km/s at z_eff = 0.1, depending on whether A ∝ H² or ∝ H. That matches Migkas's "~900 km/s to ≥500 Mpc".
  - A super-horizon tilt mode (Turner 1991) limits the shear to σ/H ≲ 1e-5 (c/HR)² = 9e-7. So a tilt cannot make the 2ℓ/3ℓ part.
  - (d) Survives for the dipole only, R² ≈ 0.27.
  - (e) Prediction: the matter number-count dipole (radio sources, quasars) differs from the CMB kinematic dipole along the Migkas direction (l, b) ≈ (280°, −15°).
- **7a. Rotation (transverse |ΔV| is invisible to redshift).**
  - Needs ω/k = 0.54. But Migkas measures redshift (line-of-sight), so it cannot see rotation either.
  - Transverse Doppler (v/c)² would be 3.4 at R. Vorticity scales as a^-2, which gives ω/H ≈ 600 at emission. Dead.
  - Caveat: Chris's |ΔV|/s equals the line-of-sight rate only to O(β²).
- **7b. Core field biasing the cluster distance indicator** (gravitational redshift of the X-ray temperature, lensing). The effect is ~ (Mc/Mcrit)(Hs/c)²/2, about 1e-5–1e-4 at 1 Gpc, against the ~5% needed. It is also always (R/s)² ≈ 200–800× smaller than its CMB counterpart. Dead.
- **7c. Special observer position, or core + cavity cancellation.** a2 does not depend on D, and a single monopole always has a P2 tidal pattern. Setting Mc = M_hole removes the field at Earth, and with it the A anisotropy: the same Meff drives both. Dead.
- **7d. New: local, gravitating, sheared flow (growing mode, Ω_b = 0.049, Λ = 0).**
  - The structure is local, so there is no SW at emission. The CMB sees only the late ISW, a2 ≈ 2|1−f|(HL/c)·Φ2 with Φ2 = (3Ω/2f)(σ/H)(HL/πc)², i.e. roughly δ·(HL/c)³.
  - With σ/H = 0.036–0.05 the CMB passes if L ≲ 0.56–0.63 Gpc. Examples: 7e-6 at L = 0.5 Gpc; 1.1e-4 (11×) at L = 1.25 Gpc.
  - Baryon-only gravity helps here: Φ per unit shear scales as Ω/f ≈ 0.3.
  - Cost: δ ≈ 0.7–1.0 over about 0.5 Gpc, and the "core" must be local, not at D = 33–94 Gpc.
  - This is an order-of-magnitude estimate, good to about a factor of 3.

## Ranked shortlist
- **1. Local growing-mode sheared flow (7d).** Passes for L ≲ 0.6 Gpc, and Migkas's own bulk-flow depth is ≥0.5 Gpc. Needs δ ≈ 1 on 0.5 Gpc: a KBC-void-like structure or its overdense counterpart. Test: Gpc-scale density contrast in galaxy counts aligned with the Migkas axis, and the Migkas anisotropy dying beyond z ≈ 0.15.
- **2. Tilt for the dipole or sine part (6),** plus kinematic geometry (5). Clean on the CMB, but R² only 0.27–0.34 on its own. Test: number-count dipole ≠ CMB dipole.
- **3. Instant potential halving, or Tolman + static potential (2/3).** Only with 1e-3–1e-6 tuning plus acausal or extreme early-universe assumptions (a_i ≲ e^-20). Kept only as a flagged extension.
- **Everything built on a far (D ≥ 33 Gpc), constant-mass core fails by 1e3–1e5.** The kinematic identity plus the spherical no-go means there is no velocity-only escape.

## Cycle 2 should compute
- The CoE geometry sampled at Migkas depth (0.3–0.6 Gpc) with a local core or overdensity (D ~ 1–2 Gpc).
  - The mass should grow in the Ω_b = 0.049, Λ = 0 growing mode, not stay constant.
  - Refit A(ℓ) against Migkas (does R² ≈ 0.8–0.9 survive?).
  - Ray-trace SW + ISW + Rees–Sciama to R = 14.3 Gpc, with emission far outside the structure.
  - Check δ against survey homogeneity and the KBC-void literature.
- Tilt + geometry joint fit: a dipole from the tilt plus 2ℓ/3ℓ from a local structure. Report the CMB number for each part separately.
- If Chris wants to keep the far core: a nonlinear (LTB/Szekeres) SW + ISW check that the −1.51 ISW/SW ratio and the 1/a history survive outside linear theory.

## Caveats
- All potential numbers are linear, while GMc/(c²D) is 0.1–2.5 today, so the far-core numbers are only order-of-magnitude. The frozen/static numbers are charitable lower bounds.
- β comes from the Fourier series of the digitised Migkas 2020 Fig. 6, where the 2ℓ term dominates. Migkas 2021 describes the anisotropy as "nearly dipole". If the true pattern is more dipolar, candidate 6 improves.
- The candidate 7d estimate uses a single-mode Poisson/continuity relation and a late-ISW factor, good to about 3×. It is not a ray-trace.
- Candidate 5 fits use Chris's code unchanged, including the shell cavity. The cavity-driven R² = 0.82 is reported but excluded from the CMB model per instruction.
