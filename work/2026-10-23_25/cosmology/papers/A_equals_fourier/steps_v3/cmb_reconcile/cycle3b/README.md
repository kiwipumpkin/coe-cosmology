# Cycle 3b: distance-dependent A(ℓ) test (replaces the cycle-3 fixed-orientation test)

## Retraction
- The cycle-3 "common orientation" test is retracted. It held one A(ℓ) curve fixed across all slices (far core R² −0.16 / 0.11 / 0.01).
- Chris is right that in his geometry A = |ΔV|/separation depends on where the sources sit. The correct test evaluates the model at each slice's own cluster distances, which is done below.

## Method
- **Code:** Chris's own Newtonian integrator, exactly as coe_steps_v3.ipynb Step 2 (ptshell + prof_core.run).
- **Parameters:** the full-sample fits, held fixed with nothing refitted: D, Mc, shell edges, l0 and the sine (amp, φ).
  - Ω_s free: D = 93.7 Gpc (thin_b_free_hug.json).
  - Ω_s = 0.049: D = 49.2 Gpc (universe_newt_omb.json).
- **Only the sampling distance changes** (DR = D/Rs).
  - Sanity check: at Rs = 14.3 Gpc the code gives R² 0.907 (free) and 0.893 (Ω_s = 0.049) against the digitised Fig. 6, matching the headline.
- **LTB / relativistic fit: not run.** Cycle 3 showed that fit is not a regular spacetime: 1+2E < 0 within about 18 Gpc of the CoE.
- **Cluster distances:** r = cz / 71.6 km/s/Mpc (the model ⟨H⟩), using Migkas 2020 per-cluster redshifts.
  - Slice medians: z<0.03 at 0.07 Gpc, 0.03–0.07 at 0.21 Gpc, z>0.07 at 0.47 Gpc, all at 0.31 Gpc.
- **Three model constructions, each normalised to the slice mean:**
  - (P) point curve averaged over the slice's cluster-distance distribution (Chris's convention: model curve vs binned data);
  - (PB) the same curve binned in 90° longitude windows like the data;
  - (C) per-cluster: each cluster gets A(ℓ_i; r_i), then averaged in the same 90° bins as the data.
- **Data:** the rebuilt Migkas 2020 slice profiles from cycle3/slices2.json (L–T slope fixed per slice).

## Results
- **The model barely moves at cluster distances.**
  - Over Rs = 0.02–2 Gpc, the predicted maximum stays at ℓ ≈ 18–19° (free) and 20–21° (Ω_s = 0.049). The minimum stays at 264–265° and 274–275°.
  - The 2ℓ phase stays at 12.8° and 13.0°.
  - The pattern starts to shift only when the sampling distance is several Gpc: at 10–14.3 Gpc the maximum moves to 24–30° and the minimum to 265°.
  - Reason: for Rs ≪ D, |ΔV|/s is set by the local velocity-gradient (tidal) tensor at Earth, which does not depend on Rs.
- **R² per slice, no refitting.** Values are point curve (P), with per-cluster binned (C) in brackets.
  - Ω_s = 0.049:
    - z<0.03: −0.28 (−0.19).
    - 0.03–0.07: 0.22 (0.24).
    - z>0.07: −0.46 (−0.62).
    - All: 0.73 (0.61).
  - Ω_s free:
    - z<0.03: −0.62 (−0.47).
    - 0.03–0.07: 0.39 (0.40).
    - z>0.07: −2.82 (−2.83).
    - All: −0.13 (−0.13).
  - The Ω_s-free headline's 0.909 depends on sampling at 14.3 Gpc. At real cluster distances its full-sample R² against the rebuilt profile is −0.13, and against the digitised Fig. 6 it is 0.10.
  - The Ω_s = 0.049 fit holds up better on the full sample: 0.73 against the rebuilt profile, 0.79 against the digitised figure.
- **2ℓ phase, predicted vs observed.** The observed errors are rough, from amp2 bootstrap errors.
  - z<0.03: predicted 13°, observed 60° ± 18°. But the 2ℓ is not significant in this slice (p_iso 0.54).
  - 0.03–0.07: predicted 13°, observed 154° ± 8°. That is 39° off (modulo 180°), several σ.
  - z>0.07: predicted 13°, observed 15° ± 14°. Consistent.
  - All: predicted 13°, observed 8° ± 10°. Consistent.
  - Per-cluster sampling alone shifts the model phase to 10–18°. It cannot produce 154°.
- **Max/min positions, predicted vs observed.** The data profiles have two nearly equal maxima about 180° apart (a 2ℓ shape), so the raw max/min flips between them; the 2ℓ phase is the robust measure.
  - Ω_s = 0.049 prediction: maximum about 20°, minimum about 275° for every slice.
  - Observed (raw / 3-harmonic fit):
    - z<0.03: max 220° / 210°, min 10° / 12°.
    - 0.03–0.07: max 0° / 340°, min 255° / 231°.
    - z>0.07: max 195° / 197°, min 115° / 285°.
    - All: max 175° / 179°, min 255° / 264°.
  - The model has a secondary maximum near 180°. So the full-sample and z>0.07 max positions are the model's second peak, not a shift.
- **Bottom line:**
  - At the cluster distances actually sampled (0.02–1.3 Gpc), Chris's geometry predicts essentially no movement of the extrema or the 2ℓ phase between slices (under 2°).
  - The observed slice-to-slice change (2ℓ phase 13° → 154° → 15°) is therefore not explained by distance dependence. It looks like noise or a real local structure in the 0.03–0.07 slice.
  - The deep slice and the full sample agree with the predicted phase. Their R² values are poor to moderate, because the predicted amplitude and shape (sampled at cluster distances) differ from the 14.3 Gpc headline.

## The sine term
- It was held fixed (multiplicative 1 − amp·sin(ℓ − l0 + φ), amp −0.124 free and −0.062 Ω_s = 0.049).
- If it represents a solid-body rotation, its |ΔV|/s contribution is independent of distance, so holding it fixed is correct.
- If it represents a tilt, i.e. an Earth peculiar velocity, its contribution scales as v/s ∝ 1/Rs. Moving from 14.3 Gpc to 0.2 Gpc then multiplies it by about 70, giving amplitudes of about 4–9. That is physically impossible for the data, so a tilt reading of the sine is not compatible with cluster-distance sampling.

## Files
- `model_rs.py` → `af_tables.json`: Chris's integrator at 50 sampling radii, 0.02–14.3 Gpc.
- `compare.py` → `compare.json`: per-slice curves, R², extrema, 2ℓ phases, and the continuous curves.
- `make_fig.py` → `cycle3b_summary.png` (viewed):
  - top row: slice data with the model at the slice's distances (point and per-cluster);
  - bottom row: predicted max/min vs sampling distance with the data points, 2ℓ phase vs distance, and R² per slice.

## Caveats
- Chris's model is 2D (Galactic plane). Clusters at all latitudes are assigned by longitude only, as in the headline.
- Distances use the low-z Hubble law with the model ⟨H⟩ = 71.6.
- Data are the Migkas 2020 L_X–T sample only. The slices are small (32 / 113 / 168), and phase errors are approximate.
- R² follows Chris's convention (model point curve vs binned data). Binned and per-cluster variants are also given, and they change R² by ≤ 0.15 except as noted.
