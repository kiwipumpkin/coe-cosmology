# A smooth density field inside our visible sphere instead of point masses

**Question (Chris):** instead of single points, is there a distribution of mass inside our visible universe (slightly below average in some places, slightly above in others) that would smooth the CMB?

**Answer:** yes, in this model.

**How weak it can be:**
- A smooth field with **rms δ ≈ 1.0–1.4%** and **peak |δ| ≈ 5–7%** reaches the 1× flat-CMB target. That covers smoothing at 1–4 Gpc.
- The field's average is zero, so no mass is added overall.
- It moves only 2–5e51 kg, which is 0.2–0.5% of the sphere's mass. The point-mass fits needed about 1e53 kg.

**Where it sits:**
- Almost all of it lies just inside the CMB ring, close to the ring plane.
- It alternates under-dense and over-dense around the sky.

**Comparison with standard cosmology (ΛCDM):**
- That is **about 4–80× stronger** than ΛCDM's typical fluctuations today at 1–2 Gpc, and **about 60–370× stronger** at 4 Gpc.
- And it must be arranged in one specific pattern, not at random.

## Method
- **Same machinery:** same residual R (237×), same 1e-5 floor, same Newtonian kernel (`ss_core.unit_field`). Code reused from `../inside_core.py`.
- **Grid:** a continuous δ(x) on a 0.3 Gpc grid, kept between 0.5 Gpc from us and 13.8 Gpc (0.5 Gpc off the ring). Mass per cell = exp-lock density × δ × cell volume.
- **Data used:**
  - The ring pattern's Fourier modes m = 1–20. The checked difference from the full 360-point calculation is 0.01×.
  - The m > 20 part of R is 0.0014×.
- **Solver:**
  - The pattern is linear in δ, so this is a regularized linear inverse problem: minimize |R − pattern|² + λ·|w|².
  - The field is δ = S·w, where S is a Gaussian smoothing of σ = s Gpc. s = 0 means no extra smoothing (grid scale).
  - A **hard zero-net-mass row** keeps the average density unchanged, so no compensation is needed.
  - λ is swept to trace the trade-off curve, then bisected to hit exactly 1×.
- **Checks:**
  - **Grid convergence:** at s = 2, a 0.4 Gpc grid gives rms 0.0128, against 0.0124 on the 0.3 Gpc grid.
  - **Bounds:** δ never comes near −1 or +few at the 1× solutions, so the "cosmic-web" bounds [−1, +few] are never active. The solution is the unbounded one.
- **ΛCDM σ(R)** is computed, not approximated:
  - Spectrum: Eisenstein–Hu no-wiggle transfer function, Ωm 0.31, Ωb 0.049, h 0.70, ns 0.965, σ8 = 0.81 (checked: it reproduces 0.810).
  - It is linear and at z = 0. Both Gaussian and top-hat windows were computed (`lcdm.py`, `lcdm_sigma.txt`).

## Results at 1× (`tradeoff_leftover_vs_rms_delta.png`)
|smoothing|rms δ|peak |δ||mass moved|ΛCDM today, rms δ (Gaussian / top-hat at that scale)|needed ÷ ΛCDM|
|-|-|-|-|-|-|
|none (0.3 Gpc grid)|0.0083|0.20|1.3e51 kg (0.12%)|0.0058 / 0.022|0.4–1.4×|
|1 Gpc|0.0105|0.068|2.4e51 kg (0.23%)|5.9e-4 / 2.8e-3|4–18×|
|2 Gpc|0.0124|0.055|3.5e51 kg (0.33%)|1.5e-4 / 7.9e-4|16–82×|
|4 Gpc|0.0141|0.059|4.9e51 kg (0.47%)|3.9e-5 / 2.2e-4|64–370×|

- **The curve is a cliff.** Below about 0.5% rms the field does almost nothing. Between 0.5% and 1.4% it drops from 237× to 1×.
- **Going further is possible.** About 4% rms with 4 Gpc smoothing reaches 0.01×.

## What the field looks like (`map_minimal_field_2Gpc.png`, `sky_delta_vs_longitude.png`)
**Where it sits:**
- At 2 Gpc smoothing, 63% of the field's δ² is at r > 11 Gpc and 64% is within 2 Gpc of the ring plane.
- With no smoothing, it is 96% at r > 11 Gpc, packed right against the ring.
- Above and below the plane (|z| > 4 Gpc), δ is only 0.5% rms.

**Shape around the sky:** eight lobes around the ring, a mainly m = 3–5 pattern. That is simply the residual's own shape (m = 3, 4, 5 at 200, 185, 173× floor) painted on the inside edge of the ring.
- **Over-dense (red)** at ℓ ≈ 15°, 128°, 195° (the CoE direction) and 263°. These are exactly where the inside point-mass fits put their 4 masses.
- **Under-dense (blue)** at ℓ ≈ 75°, 162°, 228° and 315°. These are **exactly the four-mass longitudes of the outside fits.**

**Why the signs flip:** a mass beyond the ring pulls the nearby ring points outward. The same density placed just inside the ring pulls them inward. So an inside **deficit** does what an outside **excess** does in the same direction. The smooth field, the inside point masses and the outside four masses are three versions of the same pattern.

## How it compares with ΛCDM
- **At the Gpc scales where it lives, a "smooth" field must be 4–370× stronger than standard fluctuations today**, depending on the smoothing scale and window. And it must be coherent, phased to the residual's m = 3–5 pattern.
- **The unsmoothed version** (0.8% rms, 20% peaks, on ~0.3 Gpc scales) sits within ΛCDM's amplitude range at that scale. But it is a thin, finely arranged layer pressed against the ring, so it depends on how close we let it get:
  - 0.5 Gpc gap: rms 0.0083
  - 1 Gpc gap: 0.012
  - 2 Gpc gap: 0.025
  - 4 Gpc gap: 0.106, with peaks above 1 (δ < −1, unphysical)
  - With 2 Gpc smoothing over the same gaps: 0.012 → 0.017 → 0.032 → 0.128.
- **Distance makes it worse.** Read in standard cosmology, 11–13.8 Gpc comoving is roughly z ≈ 20–700 (approximate). Linear fluctuations there are smaller than today's by roughly the growth factor, about 1/(1+z). So the shortfall against ΛCDM is **much larger** than the z = 0 numbers in the table, by roughly 20–700×.
- **Random ΛCDM structure:** I ran 20 random ΛCDM linear fields (z = 0 amplitude, ≥ 0.3 Gpc scales) through this kernel.
  - They give about 11× in m = 1 (our own motion) and only about 1–2× floor in each m ≥ 2. The total is 15× median, or about 5× without m = 1.
  - So ordinary structure produces CMB wiggles around the floor level. It does not produce the specific 200× m = 3–5 pattern. Cancelling that pattern needs an arranged field.

## Caveats
- **Missing physics, as before:** the kernel is Newtonian tidal/velocity only. For a field next to the last-scattering surface along our line of sight, time-varying potentials (ISW / Rees–Sciama), lensing and the potential difference between ring and us (Sachs–Wolfe) are left out. They are probably not small.
- **Not unique:** the minimum-rms solution is just the smallest field that works. Many other fields also reach 1×.
- **Only the ring plane is constrained:** the fit sees only the in-plane ring, so nothing constrains the field far above or below the plane. A full-sky fit would likely need more structure.
- **Values are rounded:** the "needed ÷ ΛCDM" ratios are rounded and compare rms values over the whole sphere's volume, while the field is concentrated near the ring.

## Files
- **Charts:** `tradeoff_leftover_vs_rms_delta.png`, `map_minimal_field_2Gpc.png`, `sky_delta_vs_longitude.png`
- **Data:** `tradeoff.json`, `checks.json` (where the field sits, gap sensitivity, random-ΛCDM test), `shape.txt`, `lcdm_sigma.txt`, `delta_H0.3_s*.npy`
- **Code:** `sf_core.py`, `solve.py`, `run_tradeoff.py`, `checks.py`, `lcdm.py`, `plots.py`
