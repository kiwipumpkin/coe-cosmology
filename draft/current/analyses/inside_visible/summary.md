# Can masses and voids inside our visible sphere cancel the leftover?

**Short answer:** yes. In this model, structure inside 14.3 Gpc gets below the 1× flat-CMB target with fewer objects and about the same mass as the outside fits. The best results:
- 7 masses: 0.62×
- 4 masses + 4 deep voids: 0.18×
- 4 masses with a looser per-mass cap of 1e53 kg: 0.47×

It is **more** efficient than structure outside the sphere. Caveat: for objects between us and the CMB, the model leaves out effects (ISW/Rees–Sciama, lensing, line-of-sight potential) that would probably be as large as the effect it does include (see Caveats).

## Setup (same as ../solution_space except where noted)
- Same leftover R (237.3×), floor and kernel (`ss_core.unit_field` = `lattice_core.field3`). Fitter: least squares with 4 numbers per object and random starts.
- **Region:** the object's center at distance d from us, with d ≥ core + 0.1 Gpc (we are not inside any core or void) and d + core ≤ 14.3 − 0.1 Gpc (no core touches the ring sphere). Any direction is allowed. The CoE interior can't reach inside the sphere.
- **Objects:** compact masses (δ = +100 core, 1e49 to 2e52 kg) and deep voids (δ = −0.7, uncompensated, Rv 0.5–4 Gpc). Shallow voids were dropped because they did nothing outside.
- **Mass budget:** the total concentrated mass must stay under the visible sphere's mass, 1.04e54 kg (12,249 Gpc³ × 2.9e-27 kg/m³). It never came close to binding.
- **Starts:** 40 random starts for N ≤ 6 (every split) and for masses-only up to 12. The 12 mixed classes at N = 7–8 got 16 starts. The cap-1e53 and compensation variants got 12 starts each.

## 1. Kernel check (is it valid for sources inside the ring?)
- `field3` is plain Newtonian point-mass physics. It takes the gravitational pull of the mass on each ring point, subtracts its pull on us, projects onto the line of sight, and removes the mean. Nothing in it assumes the source is outside the ring.
- I checked it directly against an independent calculation of (pull at ring point − pull on us)·line of sight for 200 random interior positions between 0.5 and 13.5 Gpc. The correlation was ≥ 0.99999999 every time.
- So the **formula is valid at any position** except right on the ring or right at us, which the margins exclude.
- **What it leaves out** (it was built for sources beyond the ring, where these effects are small):
  - **Time-varying potentials (ISW / Rees–Sciama):** CMB photons cross these potential wells on the way to us.
  - **Gravitational lensing** of the ring.
  - **Potential difference (Sachs–Wolfe-type term):** between ring points and us.
  - **Calibration:** the conversion factor KG was calibrated on exterior clumps.
  - For objects of 1e52–1e53 kg sitting 4–8 Gpc away along the line of sight, these effects are probably **not** small. Treat the inside results as "what the tidal/velocity model alone says."

## 2. Compensation for net mass
- **Local** (main run): the net mass is removed smoothly from the visible sphere, following the exp-lock density.
  - Analytically, the uniform part of a ball centered on us gives **exactly zero** pattern (shell theorem: the ring is pulled equally everywhere and we feel nothing).
  - The ~1.5% density gradient across the sphere gives a pure dipole: −(9/10)·G·M/(R·L) along the CoE axis, about 50× floor per 1e53 kg.
  - A Monte Carlo check agrees (correlation 0.9998).
- **None:** gives the same best values for masses only (7.27, 5.98, 4.25, 0.61, 0.09 for N = 4–8). So local compensation barely matters.
- **Surround** (the mass is taken from 14.3–30 Gpc, as in the outside fits): this is harder. Masses only reach 15.4× at N = 4 and 6.35× at N = 8.
  - That is still better than outside masses only (35.6× and 13×).
- **Robust case:** 4 masses + 4 voids (0.18×) is almost exactly net zero (+5.4e52 / −5.4e52 kg). It scores 0.3× under all three compensation choices.

## 3. Results (`best_vs_N_inside_vs_outside_logy.png`)
Best leftover (× floor). For mixed classes, the best split is shown. Outside N > 8 comes from the other job's add-one-at-a-time runs (`../solution_space/extend_*.json`).

|N|inside masses|inside masses + voids|inside voids only|inside masses, cap 1e53|outside masses|outside masses + voids|
|-|-|-|-|-|-|-|
|3|156|159|162|156|121|92|
|4|**7.3**|55|15|**0.47**|36|75|
|5|6.0|6.4|11|0.43|29|22|
|6|2.2|2.9|5.5|0.09|23|16|
|7|**0.62**|1.08|2.9|0.09|17|7.4|
|8|0.52|**0.18** (4+4)|0.70|0.02|13|4.5|
|9–12|0.11 / 0.06 / 0.05 / 0.02|–|–|0.03 / 0.01 (N = 9, 10)|11 / 9.7 / 8.4 / 7.7|3.0 / 1.3 / 0.58 (N = 9–11)|

**How fast it drops:**
- **First below 1×:** 7 masses (0.62×), 8 voids alone (0.70×), 4 + 4 mixed (0.18×), or 4 masses with the 1e53 cap (0.47×). Outside needs about 11 objects (mixed), and masses alone still sit at 7.7× at N = 12.
- **The big step is again at 4 masses:** 156× → 7.3× inside, compared with 121× → 36× outside.
- **Mixed classes don't beat masses alone below N = 8.** The fit space is rugged, and with N ≤ 4 some splits stall (e.g. 3+1 → 61×). Voids work as negative masses. Deep voids alone get to 15× with 4 and 0.70× with 8.

**What's left in the best fits:** only m = 7–12, and below the floor. For example:
- 4 masses: m = 7 at 6.2× and m = 6 at 4.9×; m = 1–5 are all ≤ 1.9×.
- 4 + 4: every m ≤ 0.14×.
- 7 masses: every m ≤ 0.6×.

**Distinct solutions** (2 Gpc / 0.25 dex tolerance, mirror and up/down flips merged):
- 4 masses: **1** basin (unique).
- 5 masses: 16 distinct below 10×.
- 7 masses: 35 distinct below 10×, 4 below 1×.
- 8 masses: 12 below 1×.
- 12 masses: 26 below 1×.
- Past about 5 objects nearly every start finds its own solution, the same overfitting regime as outside: 11–24 effective constraints against 4 numbers per object.

**Objects at the bounds:**
- The mass cap binds often: 3 of 4 masses at 2e52 (N = 4), 2 of 7 (N = 7), 6 of 8 (N = 8).
- About half of the voids sit near the 4 Gpc size cap.
- Almost nothing sits on the ring margin (one mass in the N = 5 best fit). A few voids or masses touch the "not on us" margin in some mixed fits.
- With the cap raised to 1e53, 2 of the 4 masses still sit at the cap (0.47×). From 5 masses up, none do.

## 4. Mass needed vs the visible sphere's budget (1.04e54 kg)
|fit|concentrated mass|void deficit|fraction of sphere's mass|
|-|-|-|-|
|4 masses, 7.3×|7.2e52 kg|–|6.9%|
|7 masses, 0.62×|9.2e52 kg|–|8.8%|
|4 masses + 4 voids, 0.18×|5.4e52 kg|−5.4e52 kg (voids Rv 3.5–4 Gpc)|5.2% + 5.2%|
|4 masses, cap 1e53, 0.47×|3.0e53 kg|–|29%|
|outside, for comparison: 4 + 4, 4.5×|6.2e52 kg|−6.3e52 kg|–|

So inside needs about the same kind of mass as outside (several e52 kg, 5–10% of all matter within the sphere), but gets 10–25× lower leftovers with it.

## 5. Why inside works better (lever arm; `lever_arm.txt`)
- **Smooth multipole spectrum.** A mass inside the ring at distance d produces multipoles that fall off smoothly with m, roughly as (d/14.3)^m. The residual is a smooth cosine series peaked at m = 3–5.
  - Useful m = 3–6 signal per 1e52 kg (median over directions): 54× at 6 Gpc, 108× at 8, 148× at 10.
  - For comparison: about 70× at 20 Gpc and 46× at 25 Gpc outside.
  - Outside masses close to the ring give sharp spikes that spill into high m, and distant ones are weak.
- **Every direction is open.** Outside, the CoE interior blocks the ℓ ≈ 162–228° in-plane spots where the residual wants mass. Inside, nothing is excluded.
  - The unique 4-mass inside solution puts masses exactly on the CoE axis (ℓ = 195°, d = 7.6 Gpc, and ℓ = 15°, d = 5 Gpc) plus a mirror pair at ℓ = 128° / 263°. That fits the residual's symmetry about ℓ = 15/195° exactly.
- **The dipole cancels itself.** An interior mass gives a large m = 1 term through its pull on us: thousands × floor per 1e52 kg at 4–6 Gpc. Opposite masses (195° vs 15°) cancel it. So nearly all of the remaining signal goes into m = 3–6.

## Caveats
- **Missing physics for interior objects** (ISW/RS, lensing, SW potential difference; see section 1). These would probably change the numbers substantially.
- **Observational check:** in standard cosmology, 4–8 Gpc comoving is roughly z ≈ 1–3, inside the range of galaxy, quasar and Lyα surveys. These are δ = +100 cores of ~0.8 Gpc radius holding 1e52 kg (~5e21 M☉), and 3.5–4 Gpc voids. Such objects would have to be checked against those surveys and against CMB lensing maps. Real structure that big is not observed at that scale.
- **Compensation choice matters for masses-only fits.** Taking the mass from outside the sphere roughly doubles the leftover. Net-zero mass + void mixes don't depend on it.
- **Fewer starts at larger N.** N = 7–8 mixed got 16 starts and the cap/compensation variants got 12, so distinct counts there aren't comparable to the 40-start classes.

## What was cut
- Shallow voids.
- Mixed classes beyond N = 8.
- Cap-1e53 mixed fits beyond 4+1 and 4+2 (0.48× and 0.44×).
- The no-compensation and surround-compensation variants ran masses only up to N = 8.

## Files
- **Charts:** `best_vs_N_inside_vs_outside_logy.png`, `inside_universe_CCCCVVVV.png` (best mixed, 0.18×), `inside_universe_CCCCCCC.png` (fewest masses below 1×, 0.62×).
- **Data:** `results_{cap2e52,cap1e53,nocomp,surrcomp}.json`, `raw_*/`, `details.txt`, `lever_arm.txt`.
- **Code:** `inside_core.py`, `fit_in.py`, `analyze_in.py`, `plots.py`.
