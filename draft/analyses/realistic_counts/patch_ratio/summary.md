# Domains made of fine-grained over- and under-dense patches

Chris: "The mixed domains should have a certain ratio of local over and under density, rather than such large pure over or under donations."

## The model
- **Domains.** Every domain (grid cell of 0.5, 1 or 2 Gpc) is a fine mix of small patches of scale r_p. Over-dense patches have δ₊ and under-dense patches have δ₋, and both kinds are present everywhere.
- **What varies.** Only the over-patch volume fraction f changes from place to place.
- **Net contrast.** A domain's net contrast is δ = f·δ₊ + (1−f)·δ₋. The cosmic mean f₀ = −δ₋/(δ₊−δ₋) gives zero.
  - Since δ = (f−f₀)(δ₊−δ₋), a cap |f−f₀| ≤ X is the same as a cap |δ| ≤ B = X·(δ₊−δ₋).
  - So I solve once for the smallest B that reaches 1×, then convert to f and to the over:under ratio for each patch choice.

**Patch values** (r_p is the radius scale of the fits; patch volume is 4/3·π·r_p³). The "extreme" values are the observed upper-envelope trends in `../../density_vs_scale/fits_v2.json`: over_noDef for δ₊ and under_noDisp for δ₋.

|patch choice|δ₊|δ₋|cosmic mean f₀|mean over:under volume ratio|
|-|-|-|-|-|
|extreme, 30 Mpc|+2.10|−0.55|0.21|0.26 : 1|
|extreme, 50 Mpc|+0.83|−0.39|0.32|0.47 : 1|
|extreme, 100 Mpc|+0.24|−0.25|0.51|1.04 : 1|
|typical ±0.3 (30/50/100 Mpc)|+0.30|−0.30|0.50|1 : 1|

**Solver.** Same kernel, residual and zero-net-mass row as before. The box problem (|δ| ≤ B) was solved through its 41-variable dual.
- The dual gives a **certified lower bound** on the leftover.
- A complementary-slackness polish gives an **actual field** that reaches it.
- Earlier runs (`lsq_linear`) stalled near the threshold. That made the old no-mirror threshold look far too high (see the correction below).

## Minimal ratio shift to reach 1× (`leftover_vs_ratio_shift.png`)
|case|domain contrast needed, B|f shift (50 Mpc extreme)|over:under ratio shift, 50 Mpc extreme|30 Mpc extreme|100 Mpc extreme|typical ±0.3|
|-|-|-|-|-|-|-|
|inside, mirror allowed, 0.5 Gpc domains|±0.030–0.0325 (bracket)|±0.027|**+12.7% / −11.8%**|±7.5%|+31% / −23%|+24% / −20%|
|inside, mirror allowed, 1 Gpc|±0.033 (0.88× at 0.0331)|±0.027|**+13.0% / −12.0%**|±7.6%|+31% / −24%|+25% / −20%|
|inside, mirror allowed, 2 Gpc|±0.0425 (0.72× at 0.0428)|±0.035|**+17% / −15%**|±10%|+43% / −30%|+33% / −25%|
|inside + outside, no mirror, 1 Gpc|±0.050–0.055 (bracket)|±0.045|**+22% / −19%**|±13%|+58% / −37%|+45% / −31%|

- **In words:** with 50 Mpc patches at the observed extreme contrasts, the over-patch fraction has to move from 32% to about 29% or 35%. That is the over:under ratio going from 0.47:1 to about 0.41:1 or 0.53:1, which is **0.88× to 1.13× the mean**.
- **Patch size matters:** smaller patches with stronger contrasts need less ratio change. Typical ±0.3 patches need about twice as much.
- **Correction to the "mixed" run:**
  - **No mirror:** the old run said 20–60% ratio deviation. That came from an under-converged solver. With the exact dual, no mirror needs only about **1.6×** the inside contrast (±0.05 vs ±0.033), which is about 9–10% in that cluster:void bookkeeping.
  - **Inside (mirror allowed):** the old figure (6.2–6.5%) is consistent with this run.

## Compared with natural scatter (`required_vs_natural_scatter.png`)
**Per domain.** Random placement of N_p = V_domain/V_patch patches scatters f by √(f₀(1−f₀)/N_p). Required f shift ÷ that scatter (inside, mirror allowed):

|domain|30 Mpc extreme|50 Mpc extreme|100 Mpc extreme|typical ±0.3, 50 Mpc|vs ΛCDM σ at domain scale|
|-|-|-|-|-|-|
|0.5 Gpc (N_p ≈ 240 for 50 Mpc)|1.0σ|0.9σ|0.7σ|1.7σ|1.6σ (σ_ΛCDM 0.020)|
|1 Gpc (N_p ≈ 1,900)|2.9σ|2.6σ|2.1σ|4.8σ|5.2σ (0.0064)|
|2 Gpc (N_p ≈ 15,000)|10.6σ|9.3σ|7.6σ|17.6σ|23σ (0.0019)|
|no mirror, 1 Gpc|4.8σ|4.2σ|3.5σ|8.0σ|8.6σ|

- **Small domains look ordinary.** In a single 0.5 Gpc domain the needed drift is within normal randomness (about 1σ).

**Whole pattern (the decisive number).** The drift has to have the same sign across thousands of neighbouring domains, arranged as one specific sky pattern. I compared the CMB pattern that random fields actually produce, through the same kernel, with the one needed (237× floor):

|random source (inside region)|CMB pattern rms it produces|needed ÷ produced|chance of the needed pattern (40 Fourier modes)|
|-|-|-|-|
|patch white noise, 100 Mpc extreme|42× floor|5.7×|≈44σ|
|patch white noise, 50 Mpc extreme|34×|7.0×|≈53σ|
|patch white noise, 30 Mpc extreme|30×|7.9×|≈59σ|
|patch white noise, typical ±0.3|18×|13×|≈88σ|
|ΛCDM linear field|12×|20×|≈128σ|

- **No mirror** gives similar or worse numbers: patch noise 14–33× floor, ≈63–120σ.
- **The σ values are approximate:** they come from a Mahalanobis distance converted with Wilson–Hilferty. They mean "not by chance", not a precise probability.

## Coherence: same eight lobes (`map_over_under_ratio.png`)
- **Same eight lobes as before:**
  - Over-rich at ℓ ≈ 15°, 128°, 195° and 263°, each at **over:under 1.13 : 1** relative to the mean.
  - Under-rich at ℓ ≈ 75°, 162°, 228° and 315°, each at **0.88 : 1**.
- **Every domain sits at the cap.** It's an on/off pattern that fills the whole inner sphere and widens inward from the ring.
- **Ring-shell power:** 74% is in m = 3–5 (m3 25%, m4 20%, m5 29%), 9% in m = 6, and almost none in m = 1–2.
- **Correlation length:** the f-pattern's autocorrelation halves at about 2.5 Gpc and first crosses zero at about 5.5 Gpc. So it is organised over Gpc scales, not random from domain to domain.

## Bottom line
- **The shift is small.** A ±12–13% shift in the local over:under patch ratio (50 Mpc patches, observed extremes) or about ±25% (typical ±0.3 patches) flattens the CMB to 1× with the mirror allowed. Without the mirror it is about ±20% / ±45%.
- **One domain at a time, that's normal.** For 0.5 Gpc domains it is within about 1σ of random patch scatter. At 1–2 Gpc it is 2–20σ per domain.
- **The whole pattern is not.** It must be a coherent 8-lobe, m = 3–5 arrangement about 5 Gpc across. Random patch scatter produces only 1/6–1/13 of the needed CMB signal, and ΛCDM produces 1/20 (about 50–130σ). So the ratio drift has to be organised; ordinary randomness can't supply it.

## Caveats
- **White-noise scatter is an upper limit on natural randomness at large scales.** Real structure is anti-correlated on Gpc scales, so ΛCDM gives less scatter than independent patches (at 0.5 Gpc: patch noise 0.037 against ΛCDM 0.020).
- **Patch values come from the observed envelopes,** which are extremes, not typical patches. The "typical ±0.3" rows show the difference.
- **Missing physics, as before:** time-varying potentials, lensing and the ring–us potential difference are left out.
- **The fields are not smoothed beyond the domain size.** The solution is on/off at the cap.
- **What was cut (time):**
  - The final bisection for 0.5 Gpc domains, which is left as a bracket.
  - The 0.5 Gpc no-mirror case; the 1 Gpc no-mirror case is used instead.
  - A polished actual field for the no-mirror threshold. Its value is dual-certified: certainly >1× at 0.050 and reachable by 0.055.
  - The ΛCDM chance number for the no-mirror case, where the covariance was near-singular; only the rms is quoted.

## Files
- **Charts:** `leftover_vs_ratio_shift.png`, `required_vs_natural_scatter.png`, `map_over_under_ratio.png`
- **Data:** `results.json` (thresholds, curves, per-domain table), `coherent_significance.json`, `pattern_in10.json`, `thr_*.json/.npz`, `log_dual_*.txt`
- **Code:** `ps_core.py` (model, dual solver, polish), `run_dual.py`, `run_thr.py`, `finish_thr.py`, `analyse.py`, `plots.py`, `map.py`
- **Obsolete:** `scan.py`, `run.py`, `lp.py`, `log_coarse.txt`, `log_run_*`, `sol_*.npz` and `lp_*` are earlier `lsq_linear`/LP attempts.
