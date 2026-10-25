# Matter dipole notes (CatWISE / Secrest + our ptshell geometry)

Date: 2026-10-06 (America/Chicago). Geometry from `migkas-geometry-fit.json` **wls_fit**. No SR curve; no new Migkas refit.

---

## A. Literature (numbers with sources)

### A.1 Secrest et al. (2021) — CatWISE2020 quasar number-count dipole

- **Paper:** Secrest, von Hausegger, Rameez, Mohayaee, Sarkar, Colin, *A Test of the Cosmological Principle with Quasars*, ApJL 908, L51 (2021).  
  - arXiv: https://arxiv.org/abs/2009.14826  
  - Journal: https://iopscience.iop.org/article/10.3847/2041-8213/abdd40  
- **Sample:** 1,355,352 quasars (~1.36×10⁶); ⟨z⟩ ≈ 1.2; 99% with z > 0.1.
- **Dipole amplitude:** \(\mathcal{D} = 0.01554\) (quoted to five digits in the paper).
- **Direction (high-count pole):** Galactic \((ℓ,b) = (238.2^\circ,\,28.8^\circ)\).
- **Offset from CMB dipole:** \(27.8^\circ\).
- **CMB dipole used:** \(v = 369.82 \pm 0.11\) km s⁻¹ toward \((ℓ,b)=(264.021^\circ,\,48.253^\circ)\) (Planck 2020).
- **Kinematic expectation (Ellis–Baldwin):** \(\mathcal{D}_{\rm kin} \sim 0.007\) (≈ half the observed amplitude; “over twice as large”).
- **Significance vs purely kinematic null:** \(p = 5\times 10^{-7}\) → **4.9σ** (one-sided normal equivalent).
- **Clustering dipole (ΛCDM estimate in their App. A):** \(\mathcal{D}_{\rm cls} \approx 0.00024\) (b=1) — ~65× smaller than observed.

### A.2 Confirming Chris’s ~3% pole-to-pole reading

For \(N = \bar N\,(1 + \mathcal{D}\cos\theta)\):

| Quantity | Formula | With \(\mathcal{D}=0.01554\) |
|---|---|---|
| High pole | \((1+\mathcal{D})\bar N\) | \(1.01554\,\bar N\) |
| Low pole | \((1-\mathcal{D})\bar N\) | \(0.98446\,\bar N\) |
| **Pole-to-pole \(\Delta N/\bar N\)** | **\(2\mathcal{D}\)** | **0.03108 ≈ 3.11%** |
| High/low − 1 | \(2\mathcal{D}/(1-\mathcal{D})\) | 0.03157 ≈ 3.16% |
| Hemisphere mean contrast | \(\approx \mathcal{D}\) | ≈ 1.55% |

**Verdict:** Chris’s draft numbers (D≈0.0155, 2D≈3.1%, hemisphere ≈1.55%) are an accurate reading of Secrest et al. 2021. Use \(\mathcal{D}=0.01554\) when quoting the published figure.

**Caveat (source-count ≠ matter density):** The Ellis–Baldwin \(\mathcal{D}\) is a *number-count* dipole from aberration + Doppler boosting in a flux-limited sample. A kinematic FLRW contribution is *expected* at \(\sim 0.007\); the anomaly is the **excess** over that prediction. Converting to a matter-density dipole \(\delta\rho/\rho\) requires a tracer model (bias, luminosity function, selection), not a one-to-one assignment of gradient to CatWISE \(\mathcal{D}\).

### A.3 Follow-ups (radio + quasars; reviews)

1. **Secrest et al. (2022)** — joint NVSS radio + CatWISE quasars (independent samples, shared sources removed).  
   - arXiv: https://arxiv.org/abs/2206.05624  
   - Journal: ApJL 937, L31; https://iopscience.iop.org/article/10.3847/2041-8213/ac88c0  
   - NVSS: \(\mathcal{D}=(1.23\pm0.25)\times10^{-2}\), ~3× kinematic, 2.6σ alone.  
   - WISE (this cut): \(\mathcal{D}=(1.48\pm0.16)\times10^{-2}\), ~2× kinematic, 4.4σ alone.  
   - **Joint:** 5.1σ rejection of kinematic null.  
   - Shared dipole: \(\mathcal{D}=(1.40\pm0.13)\times10^{-2}\) toward \((233^\circ\pm6^\circ,\,+34^\circ\pm5^\circ)\), 27° from CMB.  
   - After subtracting kinematic pieces: residual common dipole \(\mathcal{D}=(0.86\pm0.14)\times10^{-2}\) toward \((217^\circ\pm10^\circ,\,+20^\circ\pm7^\circ)\).

2. **Dam, Lewis, Brewer (2023)** — Bayesian reanalysis of CatWISE.  
   - arXiv: https://arxiv.org/abs/2212.07733  
   - MNRAS 525, 231 (2023).  
   - \(D/D_{\rm CMB} = 2.46\pm0.18\) (free amplitude); \(2.68\pm0.23\) with free direction + ecliptic bias; ~5.7σ-class excess.

3. **2025 review (real):** Secrest, von Hausegger, et al., *Colloquium: The Cosmic Dipole Anomaly*, Reviews of Modern Physics 97, 041001 (2025).  
   - arXiv HTML: https://arxiv.org/html/2505.23526  
   - Summarizes >5σ tension between matter number-count dipole and CMB-inferred kinematic expectation; discusses systematics and alternatives.

4. **2025 reassessment:** arXiv:2511.00822 — CatWISE excess revised to ~3.3–3.6σ after clustering, shot noise, mask leakage (anomaly reduced but not removed).

### A.4 DESI DR1 dipole (2026) — velocity vs number-count

- **Wu & Xia (2026):** *Measuring Peculiar Velocity and Tomographic Redshift Dipole with DESI DR1 Catalogs*.  
  - arXiv: https://arxiv.org/abs/2608.30914  
  - Phys. Rev. D 114, 043544 (2026) [per INSPIRE].  
- **Redshift dipole (cleaner kinematic probe):** high-z QSO sample → \(v = 357.95^{+55.05}_{-48.47}\) km s⁻¹, **consistent with CMB** \(369.82\pm0.11\) km s⁻¹.  
- **Number-count dipole on same QSOs:** direction \((ℓ,b)=(276.64^\circ,\,44.20^\circ)\) near CMB, but amplitude \(\tilde A = 0.122\) ≫ kinematic \(\sim10^{-3}\); authors attribute excess to footprint / completeness / window-function leakage in DR1.  
- **Implication for us:** independent *velocity* probes can agree with CMB while *counts* remain anomalous — reinforces Chris’s caveat that source-count ≠ ρ, and that systematics matter for count dipoles.

### A.5 Kashyap, Singh, Jain (and related)

1. **Kashyap, Singh, Jain (2025):** *Matter Dipole and Hubble Tension due to Large Wavelength Perturbations*, Phys. Rev. D 112, 083524 (2025).  
   - arXiv: https://arxiv.org/abs/2504.14190  
   - Adiabatic modes \(k\sim(10^{-4}{-}4\times10^{-3})\,{\rm Mpc}^{-1}\) (Gpc scales): can contribute to number-count dipole out to z~1 but **cannot fully explain** CatWISE; superhorizon curvature \(k\lesssim0.3 H_0\) can raise local \(H_0\).

2. **Follow-on / extension:** MNRAS framework with correlated long-wavelength modes (builds on Kashyap+2025) — can enhance intrinsic matter dipole and address Hubble tension via cosmic variance without changing background ΛCDM expansion.  
   - https://academic.oup.com/mnras/article/548/2/stag582/8653934

3. **Earlier related:** Das, Sankharva, Jain (2021), arXiv:2101.11016 — superhorizon perturbation for NVSS excess dipole.

### A.6 Sign: high-count pole vs CMB vs our CoE axis

**Convention in our notes / fits** (`migkas-geometry-fit.py`, `smooth-symmetric-A.py`):

- \(\ell_0 \approx 15.16^\circ\) = **anti-CoE** in Galactic longitude (θ=0 in the integrator is anti-CoE).  
- **CoE** at \(\ell_0 + 180^\circ \approx 195.16^\circ\) (Chris’s “~193°” end).  
- Low / anti-CoE end ~13°–15° (matches Chris’s ~13°).  
- Migkas \(A/A_{\rm all}\) is a longitude-only diagnostic; we place CoE/anti-CoE at \(b=0\) for plane comparisons unless a latitude is fitted.

| Pair | Angular separation |
|---|---|
| Secrest high-count (238.2°, 28.8°) ↔ CMB (264.0°, 48.3°) | **27.8°** (paper) |
| Our CoE (195.2°, 0°) ↔ Secrest high-count | **~50°** |
| Our CoE ↔ CMB | **~76°** |
| Our anti-CoE (15.2°, 0°) ↔ Secrest | **~130°** |

**Signed match:** high-count pole is on the **CoE side** of the axis (closer to CoE than to anti-CoE). That is the favorable signed alignment Chris wants: more sources toward the CoE hemisphere. It is *not* identical to the CMB apex; Secrest itself is 28° from CMB, and our CoE (plane) is ~50° from Secrest / ~76° from CMB.

---

## B. Our matter distribution (ptshell, wls_fit)

### B.1 Fixed geometry

From `/workspace/cosmology-course/notes/migkas-geometry-fit.json` → `wls_fit`:

| Parameter | Value |
|---|---|
| \(D\) | 58.197 Gpc |
| \(m_{\rm cc}\) | 0.001712 |
| \(R_{\rm in}\) | 0.500 Gpc → \(r_{\rm in}/D \approx 0.00859\) |
| \(R_{\rm out}\) | 72.834 Gpc → \(r_{\rm out}/D \approx 1.2515\) |
| \(\Omega_s\) (Os) | 0.049 |
| \(L_c\) | 0.005 (D units) ≈ 0.291 Gpc |
| \(\ell_0\) | 15.163° (anti-CoE) |

Density law (`ptshell.py` / `density_plot.py` convention): exponential core (mass \(m_{\rm cc}\) in critical units) + **constant** \(\rho/\rho_{\rm crit} = \Omega_s\) between \(R_{\rm in}\) and \(R_{\rm out}\). No separate uniform background (Ob=0); the shell *is* the uniform component.

### B.2 Where Earth sits

**Earth is inside the thick shell:** \(R_{\rm in} \ll D < R_{\rm out}\) (0.50 ≪ 58.2 < 72.8 Gpc). The interior is nearly empty aside from the compact core (\(L_c\sim0.3\) Gpc). Across the vision sphere, both the near side (toward CoE) and the far side (anti-CoE) still sample the **same constant-\(\Omega_s\) slab** for the R values considered.

### B.3 Native contrast at R = 3 and 14.3 Gpc

Computed from the same profile the integrator uses (and cross-checked with `P.run("ptshell",…)` → built-in `contrast`). Files:

- `/workspace/cosmology-course/notes/matter-dipole-rho.csv`
- `/workspace/cosmology-course/notes/matter-dipole-contrasts.csv`
- `/workspace/cosmology-course/notes/matter-dipole-summary.json`
- `/workspace/cosmology-course/notes/matter-dipole-rho.png`

| R (Gpc) | \(r_{\rm in}=D-R\) | \(r_{\rm out}=D+R\) | In shell? | \(\rho_{\rm in}/\rho_{\rm crit}\) | \(\rho_{\rm out}/\rho_{\rm crit}\) | \((\rho_{\rm in}-\rho_{\rm out})/\bar\rho\) |
|---|---|---|---|---|---|---|
| 3.0 | 55.20 | 61.20 | both yes | 0.049 | 0.049 | **≈ 0** |
| 14.30 | 43.90 | 72.50 | both yes (out almost at Rout) | 0.049 | 0.049 | **≈ 0** |

Hemisphere averages and cos-weighted dipole of ρ on the vision sphere: also **≈ 0** (numerical noise ~10⁻¹⁶). Integrator `contrast` field: **0**.

**Interpretation:** With a flat shell and Earth deep inside it, the native model has **no few-percent matter dipole** at either vision radius. The Migkas \(A/A_{\rm all}\) anisotropy in this geometry comes from **velocity / expansion geometry**, not from a ρ gradient on the vision sphere. That is the opposite of assigning a CatWISE-like density contrast by hand.

*(If Rout were pulled inward so \(D+R > R_{\rm out}\), the outward pole would leave the shell and contrast would jump discontinuously — that is a lever for a future constraint, not the present fit.)*

### B.4 Caveat again

Even if we later match a few-percent \(\delta\rho/\rho\), CatWISE \(\mathcal{D}\) is a **tracer number-count** dipole. Path:  
\(\rho(r) \rightarrow [\rho(D-R)-\rho(D+R)]/\bar\rho \rightarrow D_{\rm matter}(z) \rightarrow\) bias/selection → \(D_{\rm count}\).  
Do **not** set \(\delta\rho/\rho = 2\mathcal{D}\) blindly.

---

## C. Fitting plan: CatWISE as an extra constraint (leave freedom)

Goal if the independent matter dipole is taken as real and CoE-aligned in sign: add a soft constraint on the **signed** local density contrast without wrecking the Migkas \(A/A_{\rm all}\) geometry fit.

### C.1 What to hold (Migkas geometry)

Hold (or tightly prior) the parameters that the WLS Migkas fit already pins:

- \(D\), \(\ell_0\) (axis), and the velocity-derived \(A/A_{\rm all}\) shape quality (\(R^2\), χ²).  
- Prefer **not** to free Amp/φ sine modulation here (geometry-only path).

### C.2 What can absorb a ~few-percent signed ρ contrast

Because the present shell is flat, contrast must come from **radial structure**, not from Os alone:

1. **\(R_{\rm out}\) (or \(k_{\rm out}\))** — move the outer edge near \(D+R\) so the anti-CoE vision pole samples lower density (outside or on a tapered edge). Note from fit JSON: “χ² is flat in Rout once Rout ≥ D+R”; so Rout has **freedom** on the Migkas side once large enough — ideal for a density-dipole prior.
2. **Soft outer taper / non-constant shell** — replace top-hat Os with a mild radial slope or exponential decline outside some \(r_*\), controlled by 1–2 extra params.
3. **\(R_{\rm in}\) / void depth** — less useful for R=3–14 (Earth already outside the tiny cavity); more relevant for very large R or early-time maps.
4. **Os** — overall normalization; **does not create a dipole by itself** if still spatially flat. Keep for mean density / critical matching.
5. **\(m_{\rm cc}\), \(L_c\)** — affect central mass and local Ht; tiny at Earth for current Lc. Secondary for count dipole; keep tied to Migkas/SR needs.

### C.3 Proposed procedure (next step, not done here)

1. Define a target **matter** dipole amplitude after tracer conversion, e.g.  
   \(D_\rho^{\rm target} \approx (\mathcal{D}_{\rm CatWISE} - \mathcal{D}_{\rm kin})/b_{\rm eff}\)  
   with \(\mathcal{D}_{\rm CatWISE}=0.01554\), \(\mathcal{D}_{\rm kin}\sim0.007\), and \(b_{\rm eff}\) a free effective bias (or a prior range ~1–2), **signed** so \(\rho\) higher toward CoE.
2. Observable from the model (same as Chris asked):  
   \(C(R) = [\rho(D-R)-\rho(D+R)]/\bar\rho_{[D-R,D+R]}\)  
   and/or cos-weighted \(D_\rho(R)\) on the vision sphere, for R ∈ {3, 14.3} (and optionally a redshift-weighted mix matching ⟨z⟩~1.2).
3. Joint objective:  
   \(\chi^2_{\rm tot} = \chi^2_{\rm Migkas}(A/A_{\rm all}) + \lambda\,(C(R) - C_{\rm target})^2/\sigma_C^2\)  
   with λ ≤ 1 and σ_C generous (CatWISE systematics + bias uncertainty). Start with λ small.
4. Free set for a first pass: **{Rout (or kout), one taper param}**; hold D, Rin≈floor, mcc, Os, ℓ₀.  
   Second pass: allow mild Os–Rout trade if mean density is constrained elsewhere.
5. **Freedom that remains:** overall Os scale; core parameters; any SR/temperature pipeline; sine multipoles; bias \(b_{\rm eff}\); which R (or n(z)) enters \(C(R)\); DESI redshift-dipole prior as a *velocity* check separate from counts.

### C.4 Success criteria for a later fit

- Signed \(C(R)>0\) (higher ρ toward CoE) at the few-percent level after tracer conversion, without dropping Migkas \(R^2\) much below the current ~0.55 (geometry-only) / without breaking the SR residual story.  
- Rout not forced below D+R in a way that reintroduces a hard edge artefact on the CMB last-scattering cartoon unless intentional.  
- Explicit statement of source-count vs ρ conversion assumptions.

---

## File index

| Path | Content |
|---|---|
| `/workspace/cosmology-course/notes/matter-dipole-notes.md` | This note |
| `/workspace/cosmology-course/notes/matter-dipole-rho.csv` | ρ(r)/ρ_crit along radial line |
| `/workspace/cosmology-course/notes/matter-dipole-contrasts.csv` | Contrasts at R=3 and 14.3 |
| `/workspace/cosmology-course/notes/matter-dipole-summary.json` | Machine-readable summary + angles |
| `/workspace/cosmology-course/notes/matter-dipole-rho.png` | One-page ρ(r) with D±R marked (legend ≥14 pt) |
