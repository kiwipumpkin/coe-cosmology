# P₂ template + footprint dipole leakage

**Sedov–Taylor anisotropic Hubble / Cosmology agent**  
**Date:** 2026-09-17 (CT)  
**Parent folder (Drive):** `1wZYStJhpdowbcjTkkS9TOImAtt9_GvzD`  
**This Drive doc fileId:** `15osiNhP5b0xtHv3bksgKHL5dk6Gd4r7moMe1F80as6s`  
**Prior:** dipole fraction `13jZfhS54BCAL5mHA6XpsuVbBre1gd2i5W_BAPUPN--k`; enclosed mass `1uXglnM-4gfvX2E4Q41vxun6Cg0D9GUPxQMEyAvO1Ob0`  
**Code/data:** `p2_footprint_leakage.py`, `p2_footprint_leakage_results.json`  
**Local:** `/workspace/cosmology/p2_footprint_leakage.md`

**Binding:** κ_natural = 0.15; time-integrated O(1); parameters unconstrained; be skeptical. Migkas 2020: anisotropy toward (l,b)~(280°,−20°) / joint **(303°,−27°)**; ~16±3% DL; partial-sky cluster sample.

---

## 0. Executive recommendation for Chris

| Question | Answer |
|----------|--------|
| Does **footprint leakage alone** from a pure P₂ sky reach Migkas **A ~ 0.10–0.16** under **realistic** masks? | **NO** |
| Max fake \|A\| (realistic: \|b\|>20°, random 30–50%, crude cluster-like) | **~0.044** (α=1) / **~0.032** (drag) — **below** Migkas floor |
| When does fake A hit Migkas scale? | **Only** CoE-facing hemisphere/cone (selection artifact): A_moment ≈ **0.113** (α=1) / **0.075** (drag); A_Migkas(M+Aμ) ≈ **0.45** / **0.30** |
| Prefer what next? | **Quadrupole (P₂) template** comparison at κ=0.15 — not dipole MCMC on spherical M_enc |

**One paragraph:** A pure tidal quadrupole at natural κ (C₂ = 2ακ = 0.30 bolt-on / 0.199 drag) with axis at Migkas joint **(303°, −27°)** has **zero** full-sky Legendre dipole. Under galactic \|b\|>20° (parity-even cut), random 30–50% draws, and a crude CoE-preferring cluster-like mask (~16% sky), fitted fake dipole amplitudes stay **≲ 0.05** — **insufficient** for Migkas ~0.1–0.16. Only a **hard CoE-facing hemisphere** (or \|b\|>20° ∩ CoE hemi) produces Migkas-sized **fitted** A, with direction locked to the true axis and spuriously high R² — the selection artifact already quantified in `dipole_fraction.md`. **Leakage under realistic extragalactic footprints does not explain Migkas.** Do not retire the morphology NO-GO on that basis.

---

## 1. Pure P₂ forward template

### 1.1 Definition

Observer-centric fractional anisotropy (tidal limit of point-mass M_enc):

\[
\frac{H_\mathrm{obs}}{H}-1 = C_2\, P_2(\mu),\qquad
P_2(\mu)=\frac{3\mu^2-1}{2},\qquad
\mu=\hat n\cdot\hat a,
\]

with

\[
\boxed{C_2 = 2\alpha\kappa},\qquad \kappa=\kappa_\mathrm{natural}=0.15.
\]

Preferred axis \,\hat a = Migkas joint galactic **(l, b) = (303°, −27°)** (default).

| α | Value | C₂ | half-ptp of pure P₂ (= ¾ C₂) |
|---|------:|---:|-----------------------------:|
| bolt-on | 1.000 | **0.300** | **0.225** |
| drag (η→0) | ≈0.664 | **0.199** | **0.149** |

P₂ ∈ [−½, 1] ⇒ C₂ P₂ ∈ [−C₂/2, C₂]; half peak-to-peak = ¾ C₂.

### 1.2 Depths d/D ∈ {0.05, 0.1, 0.2} (α=1)

Pure P₂ is the ε→0 limit; forward point-mass adds octupole etc. C₂ stays ≈0.30; half-ptp grows mildly.

| d/D | half-ptp pure P₂ | half-ptp forward | C₂ forward | C₁ forward |
|----:|-----------------:|-----------------:|-----------:|-----------:|
| 0.05 | 0.225 | 0.237 | 0.300 | ~10⁻⁸ |
| 0.10 | 0.225 | 0.251 | 0.300 | ~10⁻⁸ |
| 0.20 | 0.225 | 0.287 | 0.300 | ~10⁻⁷ |

**Drag:** scale all amplitudes ×α_drag ≈ 0.664 (same morphology). At κ_natural, Migkas-depth variation is **O(0.15–0.25) half-ptp** — amplitude-natural, quadrupole-shaped.

---

## 2. Footprint leakage protocol

**Sky:** lat–lon grid (180×90), solid-angle weights ∝ cos(lat).  
**Field:** pure P₂ only (no octupole) — isolates **mask-induced** fake dipole.  
**Fits (unconstrained):**

1. **Vector dipole:** frac ≈ M + **A** · n̂ → \|A\|, recovered (l,b), R²  
2. **Fixed-axis Migkas:** frac ≈ M + A μ, μ = n̂·â_Migkas → A_Migkas  
3. **Moment:** A_moment = ⟨f μ⟩/⟨μ²⟩ (no free monopole)

**Masks:**

| Mask | Intent |
|------|--------|
| full sky | orthogonality check (expect A≈0) |
| random 50% / 30% | sparse/partial coverage (N=40 realizations for stats) |
| galactic \|b\|>20° | standard extragalactic cut |
| CoE-facing hemisphere μ>0 | worst-case selection (prior A~0.11) |
| anti-CoE hemisphere | control |
| \|b\|>20° ∩ CoE hemi | “realistic plane cut + directional selection” |
| crude cluster-like | \|b\|>20° + mild CoE preference + random thinning → ~16% sky |

---

## 3. Results — fake A by mask

### 3.1 Bolt-on α=1 (C₂=0.30)

| Mask | sky frac | \|A\|_vector | A_Migkas | A_moment | R²_dip (fixed) | dir. bias vs (303°,−27°) |
|------|----------:|------------:|---------:|---------:|---------------:|-------------------------:|
| full sky | 1.00 | **~0** | ~0 | ~0 | 0 | (noise) |
| galactic \|b\|>20° | 0.65 | **~0** | ~0 | ~0 | 0 | (noise) |
| random 50% (1 draw) | 0.50 | 0.002 | −0.002 | −0.002 | ~0 | large |
| random 30% (1 draw) | 0.30 | 0.004 | −0.003 | −0.003 | ~0 | large |
| **cluster-like crude** | 0.17 | **0.044** | **+0.046** | **0.040** | 0.04 | ~15° |
| CoE-facing hemi | 0.50 | **0.450** | **+0.450** | **0.113** | 0.94 | **~0°** |
| anti-CoE hemi | 0.50 | 0.450 | −0.450 | −0.113 | 0.94 | ~180° |
| \|b\|>20° ∩ CoE hemi | 0.33 | 0.448 | +0.447 | 0.103 | 0.93 | ~0° |

**Random multi-realization (α=1):**

| keep | median \|A\| | p95 \|A\| | max \|A\| |
|------|------------:|----------:|----------:|
| 50% | 0.0027 | 0.0057 | 0.0061 |
| 30% | 0.0050 | 0.0090 | 0.011 |

### 3.2 Drag α≈0.664 (C₂≈0.199) — scale ≈×0.664

| Mask | \|A\|_vector | A_Migkas | A_moment |
|------|------------:|---------:|---------:|
| full / \|b\|>20° | ~0 | ~0 | ~0 |
| random 30% (typ.) | ~0.003 | ~0.003 | ~0.003 |
| **cluster-like** | **0.032** | **0.032** | **0.029** |
| CoE-facing hemi | **0.299** | **0.299** | **0.075** |

Random 30% p95 \|A\| ≈ 0.0055 (drag).

### 3.3 Why \|b\|>20° gives A≈0

P₂(μ) is **parity-even** (n̂ → −n̂). The galactic plane cut \|b\|>20° is also parity-even. The product remains even ⇒ **no odd (dipole) moment**. Random thinning breaks parity weakly ⇒ tiny A. CoE-facing hemi is **parity-odd as a domain** ⇒ large leak.

Analytic CoE hemi (μ-uniform on μ>0):  
A_moment = C₂ × ⅜ = 0.30 × 0.375 = **0.1125** — matches numerics.

---

## 4. Does leakage reach Migkas ~0.1–0.16?

| Scenario | Max fake \|A\| (α=1) | vs Migkas 0.10–0.16? |
|----------|--------------------:|----------------------|
| Full sky / \|b\|>20° | ~0 | **No** |
| Random 30–50% (p95) | ≲0.01 | **No** |
| Crude cluster-like (~16% sky, mild CoE prefer) | **0.044** | **No** (≲0.4× Migkas floor) |
| CoE-facing hemisphere | A_moment **0.113** / A_Migkas **0.45** | **Yes — but artifact** |
| \|b\|>20° ∩ CoE hemi | A_moment **0.103** / A_Migkas **0.45** | **Yes — still requires CoE-facing cut** |

**Drag:** all numbers ×~0.66; CoE hemi A_moment ≈ **0.075** (below 0.10); cluster-like ≈ **0.03**.

**Clear answer:** Under masks that resemble **standard extragalactic / partial-sky without deliberately selecting the CoE hemisphere**, pure-P₂ leakage **does not** produce Migkas-amplitude dipoles. Claiming Migkas A from quadrupole leakage requires admitting a **strong directional selection** aligned with the preferred axis — which is then a **systematic**, not a CP detection.

---

## 5. Direction bias

- **CoE hemi / CoE∩\|b\|:** recovered dipole direction locks to true axis (bias ~0°) — fake A **points at Migkas by construction**.  
- **Random / tiny A:** direction is noise (bias meaningless when \|A\|~0).  
- **Cluster-like (mild CoE prefer):** bias ~15° toward true axis, A still small.

So if an analysis finds (l,b) near (303°,−27°) **and** A~0.1 from a quadrupole-only sky, the footprint must be **strongly** CoE-facing — check the actual Migkas footprint against μ>0 before interpreting.

---

## 6. Skeptical caveats

1. **Crude cluster mask ≠ real Migkas catalog.** Real X-ray cluster footprints can be more anisotropic; a dedicated HEALPix mask of the actual sample could raise fake A somewhat — unlikely to jump from 0.04 → 0.16 without CoE-facing weight.  
2. **A_Migkas vs A_moment:** M+Aμ on a hemisphere returns A~0.45 while moment projection is 0.11. Always specify which estimator; unconstrained A is easy to inflate.  
3. **Pure P₂ omits octupole.** Prior φ-uniform octupole leak gave A~0.02 at d/D=0.1 — still sub-Migkas; does not change the realistic-mask conclusion.  
4. **Parameters unconstrained:** fitting free (l,b,A) will happily absorb leakage; report R² vs Q template.  
5. **Amplitude of P₂ itself** at κ_natural is Migkas-*scale in half-ptp* (0.15–0.25) — the issue is **morphology** (quad vs dipole), not raw variation size.

---

## 7. Program implication

| Path | Status |
|------|--------|
| Migkas (l,b,A) dipole as smoking CP from spherical M_enc | **Still NO-GO** |
| “Footprint leakage explains Migkas dipole” under realistic cuts | **NO** — does not reach 0.1–0.16 |
| CoE-hemisphere selection artifact | **YES** — A_moment~0.11 (α=1); must be checked, not celebrated |
| Honest next step | Fit **P₂ template** (κ~0.15, α_drag, axis free or fixed Migkas) to cluster H_obs / DL residuals; compare χ² to dipole |

---

## 8. Artifacts

| Path | Role |
|------|------|
| `/workspace/cosmology/p2_footprint_leakage.md` | This note |
| `/workspace/cosmology/p2_footprint_leakage.py` | Scan code |
| `/workspace/cosmology/p2_footprint_leakage_results.json` | Full numeric dump |
| Drive fileId `15osiNhP5b0xtHv3bksgKHL5dk6Gd4r7moMe1F80as6s` | Google Doc copy |
| Drive parent `1wZYStJhpdowbcjTkkS9TOImAtt9_GvzD` | Folder |

---

## 9. Bottom line for Chris

- **C₂ = 2ακ = 0.30 (0.199 drag)** at κ_natural; **half-ptp pure P₂ = 0.225 (0.149 drag)**; forward at d/D=0.05–0.2: half-ptp **0.24–0.29** (α=1).  
- **Max fake A under realistic masks:** **~0.044** (cluster-like, α=1) — **does not reach** Migkas 0.10–0.16.  
- **Max fake A overall:** CoE hemisphere **A_amp = A_Migkas = 0.45**, **A_moment = 0.113** (α=1) — selection only.  
- **Recommendation:** Leakage alone under realistic footprints **fails** as a Migkas-dipole explanation. Keep morphology CONDITIONAL-NO / NO-GO for spherical M_enc; pursue **quadrupole forward template** comparison, not dipole smoking-gun claims.

*End of P₂ + footprint leakage note.*
