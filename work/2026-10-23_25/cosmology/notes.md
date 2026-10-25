# Project Status & Development Plan
**Sedov–Taylor anisotropic Hubble — Cosmology agent / Chris**  
Last updated: 2026-09-17 (America/Chicago)

**Framing (Chris clarification):** The Universe has a **real geometric center** (center of expansion / CoE), argued to lie **beyond the observable horizon**. We are an **off-center observer** in a **centered Sedov–Taylor blast-wave cosmos**. This is **not** “ΛCDM + a small dipole patch.” The scientific goal is to **locate that center** (direction + distance scale) and constrain the enclosed mass / similarity structure — or to **kill the idea** if CMB isotropy (and related probes) make a real center untenable.

---

## 1. Honest status (no oversell)

### What is solid
- **Core ontology is clear:** preferred CoE + ST-like radial expansion (V ∝ r from the center) + observer offset D; anisotropy is geometric, not a weak FLRW perturbation.
- **Toy kinematics exist:** Google Sheet `16bDVnqOx_gaEJcDTJcuYwkkaY3uOdIZxCBdUyAHjHys` — Newtonian plane model with k, G, M, D, R; recession velocity vs angle with/without central gravity.
- **Visualization demos:** Python/Plotly recession curves, XY circles, 3D surfaces, Mollweide maps rotated toward Migkas-like axes, with noise/binning/smoothing.
- **Partial lit + outline:** Migkas et al. (2020), Sedov (1946), Taylor (1950); chapter outline for a long thesis-style treatment.

### What is sketchy
- **“Match” to Migkas is visual only** — maps rotated by hand; no χ²/MCMC; parameters unconstrained.
- **Eqs 3.11–3.12 are heuristic** — gravity → peculiar velocity → H_obs(θ) mixes length scales without a consistent derivation.
- **Beyond-horizon center is asserted, not constructed** — outline notes Hubble sphere vs particle horizon, but there is no GR matching showing how a CoE outside the particle horizon can still set local expansion anisotropy while preserving (or explaining) CMB near-isotropy.
- **ST → cosmology mapping is incomplete** — classical blast wave ≠ metric cosmology; E, ρ₀, ξ need a defined cosmological meaning.

### What is unfinished
- Manuscript ≈ **abstract + outline + partial Ch. 2–3 + demos**; §§4–8 mostly empty “≈25 pages” bullets.
- **No center localization** yet (the stated goal): no posterior on (l, b, D) or M(r).
- **CMB adversary not confronted quantitatively** — this is the make-or-break, not a side note.
- Spreadsheet = **toy Newtonian only** (often G=k=M=1).

**Bottom line:** Clear research *program* (find the CoE in an ST cosmos). Current artifacts are a Newtonian intuition pump + outline — **not** yet a falsifiable, peer-reviewable localization of that center.

---

## 2. Core claim summary (from draft extract)

**Title:** *Anisotropic Hubble Expansion from a Central-Mass Sedov–Taylor Model of the Universe*

**Claim:** Cosmic expansion is a Sedov–Taylor blast wave about a real center beyond our horizon. Isotropic ST recovers v = H d. An off-center observer plus central mass yields angle-dependent recession velocity, motivated as the origin of Migkas et al. (2020) cluster LX–T / DL anisotropy.

**Assumptions:** dust-like ST medium; preferred CoE; Earth offset by D; central mass (point or Eq. 3.13 profile); Newtonian first pass; Migkas signal treated as cosmological target for center-finding.

**Eqs 3.1–3.14 (abbrev.):** R(t)=ξ(E/ρ₀)^{1/5}t^{2/5}; v=(2/5)R/t; v⃗=H(t)r⃗; LOS projection → H_obs; g=GM/r² → dipolar H_obs(θ); extended m(r).

**Migkas role:** empirical anisotropy to *aim* the center-finding (direction ~ (280°,−20°) / joint ~(303°,−27°); ~4–5σ; ~16±3% DL in joint analysis) — not yet a likelihood for (D, M, ξ).

**Horizon argument:** CoE placed beyond observable horizon so we do not “see” the center as a local object; anisotropy is the clue. **CMB isotropy is the main adversary of that claim** (off-center observers in centered cosmologies generically imprint CMB multipoles).

---

## 3. Spreadsheet (toy only)

Id `16bDVnqOx_gaEJcDTJcuYwkkaY3uOdIZxCBdUyAHjHys`: Newtonian V=kD from CoE, observer at D, test point at R, gravity from M → recession v(θ).  
**Use for intuition / Phase A cross-check. Do not quote as cosmological parameter constraints.**

---

## 4. Phased development plan

### (A) Clean derivation + consistent notation (eqs 3.1–3.14)
Self-consistent Newtonian kinematics of an **off-center observer in a centered ST velocity field**, then gravity. Fix symbols (r vs d), vector recession, and units. Deliverable: cleaned note + spreadsheet agreement tests.

### (B) Falsifiable predictions (center-finding observables)
Predict full-sky, redshift-dependent H_obs(n̂,z) / DL(n̂,z) for a true CoE. Confront:
- Migkas-like cluster dipoles (target for direction),
- SN Ia anisotropic Hubble diagram,
- BAO directional isotropy,
- bulk-flow vs true anisotropic expansion degeneracy,
- **especially CMB** (see go/no-go).

### (C) Parameter estimation — locating the center
Fit **direction (l,b), D (or D/R), M(r) or (M, r_s), ξ / H normalization** to anisotropy data. MCMC forward-model matching Migkas binning; report posterior on CoE coordinates. First concrete step: cleaned forward model → fit amplitude + direction to Migkas A(n̂); ask whether any (D, M) gives ~O(10%) DL variation without absurd parameters.

### (D) Relativistic / light-travel corrections
Light cones, evolving positions, geodesic DL; whether a beyond-horizon CoE can be formulated in GR (LTB-like / blast-wave metric) without destroying the program.

### (E) Go / no-go (center localization program)
**Scientific goal:** locate the CoE. **Main adversary:** CMB isotropy.

**NO-GO if:**
1. **CMB (primary):** Any CoE location that fits late-time ~O(10%) anisotropy forces CMB dipole/quadrupole (or Bianchi shear) beyond Planck / Saadeh et al. (2016) limits, with no viable screening/mechanism at z∼1100.
2. **No causal GR construction** for “center beyond particle horizon” that still sources the local ST field we inhabit.
3. **Amplitude conflict:** center parameters needed for Migkas destroy BAO/SNIa isotropy at the measured level.
4. **Target evaporates:** Migkas-level signal shown to be mostly scatter/systematics (e.g. simulation-based reassessments), removing the empirical hook for center-finding.
5. **Degeneracy:** bulk flow + FLRW explains the data with better evidence than a real CoE.

**GO only if:** a CoE posterior exists that (i) fits late-time anisotropy better than bulk-flow nulls, (ii) has an explicit CMB-compatibility mechanism or prediction, (iii) forecasts distinctive eROSITA/Roman/CMB-S4 tests.

---

## 5. Immediate next actions
1. Phase A: clean eqs 3.1–3.14 + sheet consistency.
2. One-page **CMB adversary memo** (off-center centered cosmos → multipoles; Alnes & Amarzguioui; Saadeh+2016).
3. Minimal MCMC center/amplitude fit to Migkas dipole.
4. Go/no-go memo after (2)–(3).

## 6. Assets
Draft: `/workspace/cosmology/draft.txt`, `/workspace/sedov-taylor-universe-extract.txt`  
Folder: parentId `1wZYStJhpdowbcjTkkS9TOImAtt9_GvzD`  
Companion: *Literature Map*


---

# Literature Map
**Sedov–Taylor anisotropic Hubble — off-center observer in a centered blast-wave cosmos**  
Last updated: 2026-09-17 (America/Chicago)

**Framing:** Chris’s working hypothesis is a **real cosmic center** (beyond the horizon). Literature below is organized for a **center-localization** program. **CMB isotropy is the main adversary** of that claim. No cheerleading: supportive probes are listed honestly alongside critiques and null results.

---

## 1. Primary empirical hook — Migkas et al. and follow-ups

### Migkas et al. 2020, A&A 636, A15
- **Title:** Probing cosmic isotropy with a new X-ray galaxy cluster sample through the LX–T scaling relation.
- **Result:** Directional variation in LX–T normalization; ≳4σ toward (l,b)∼(280°,−20°) for 313 MCXC clusters; joint with other samples ~5σ toward ~(303°,−27°); interpreted as possible ~16±3% luminosity-distance variation (or systematics / X-ray effects).
- **Role for this project:** Best all-sky late-time anisotropy map to *aim* CoE direction. Not proof of a geometric center.

### Migkas et al. 2021, A&A 649, A151
- Ten cluster scaling relations; apparent ~9% H0 dipole toward ~(280°,−15°), >5σ in their MC tests; **or** ~900 km s⁻¹ bulk flow to ≳500 Mpc — **degenerate** until higher-z clusters (eROSITA).
- **Implication:** Even if real, signal may be peculiar velocity, not anisotropic expansion from a CoE.

### Follow-ups / related claims
- Pandya et al. 2024 (cluster σv scaling): H0 dipole claims consistent in direction with Migkas, amplitude sensitive to scatter.
- Boubel et al. 2025 (Tully–Fisher): smaller ~3% H0 variation, direction roughly nearby — weaker than cluster LX–T.

### Critiques / simulation reassessment
- **FLAMINGO-based analyses (e.g. A&A 2025 work on scaling relations as isotropy tracers):** isotropic ΛCDM lightcones + realistic scatter can produce large apparent anisotropies; reported significance of Migkas-like signals **drops** (order ~3σ class, much of “signal” from noise/scatter). Treat Migkas amplitude as **upper target**, not gospel.
- Always separate: (a) detection of directional variation in fitted normalizations vs (b) cosmological anisotropic expansion vs (c) bulk flow vs (d) X-ray absorption / calibration / selection.

---

## 2. Other cosmic dipole / anisotropy claims

| Probe | Claim (peer-reviewed / widely cited) | Tension with CoE program |
|-------|--------------------------------------|---------------------------|
| Cluster bulk flows / “dark flow” | Kashlinsky et al. (ApJ 2008+): large coherent flows via kSZ; controversial, contested systematics | Could mimic CoE dipole without a blast-wave center |
| Galaxy peculiar velocities | Watkins, Feldman & Hudson (MNRAS 2009): ~407 km s⁻¹ in ~50 h⁻¹ Mpc window — large vs vanilla ΛCDM but << Migkas 900 km s⁻¹ scale | Sets scale of known flows |
| SN Ia dipoles | Mixed; some residual dipoles after CMB-frame correction; often consistent with local structure | Must predict μ(z,n̂) for off-center ST observer |
| CMB large-angle anomalies | Planck 2018 VII (A&A 641, A7): ΛCDM-consistent stats + known anomalies (hemispherical asymmetry, low-ℓ alignments) | Anomalies ≠ license for O(10%) late H dipole from a CoE |
| Quasar / radio dipole | Several claims of excess dipole vs kinematic expectation; contested | If real, another preferred axis — must align with CoE or explain misalignment |
| Bianchi / anisotropic expansion limits | **Saadeh et al. 2016, PRL 117, 131302:** Planck T+P; anisotropic expansion strongly disfavored (odds ~121000:1); (σ_V/H)_0 < 4.7×10⁻¹¹ (95%) | **Central adversary** for global anisotropic expansion |

---

## 3. Centered cosmologies & off-center observers (theory map)

These are the closest peer-reviewed precedents to “real center + we are offset”:

- **Alnes & Amarzguioui 2006, Phys. Rev. D** (astro-ph/0607334): Off-center observer in spherically symmetric inhomogeneous (LTB-like) universe → **CMB dipole constrains offset to ≲~15 Mpc** for their void models; induced quad/oct too small to explain low-ℓ alignment.
- **Alnes & Amarzguioui 2006/2007** SN Hubble diagram for off-center observers: SN data do not tighten the CMB dipole bound much.
- **Kodama et al. 2010** (arXiv:1004.3089): Analytic dipole/quadrupole for off-center observers in general spherical symmetry; update WMAP constraints on observer location.
- **García-Bellido & Haugbølle** and related void literature; Biswas, Notari, et al. void tests (CMB+BAO+SN+H0): large voids that replace dark energy face **kSZ, BAO, CMB** limits; fine-tuning to sit near center.
- **Classical ST:** Sedov 1946; Taylor 1950, Proc. R. Soc. A 201, 159 — blast-wave similarity; **not** cosmological metrics. Mapping ST → expanding universe remains an open modeling step for this project.

**Honest reading for Chris’s thesis:** LTB/void literature shows that **centered models with off-center observers are already tightly squeezed by the CMB dipole**. A CoE *beyond the horizon* is a stronger claim than a Gpc void with us near center — it needs a **new** mechanism for why last-scattering looks nearly isotropic while z≲0.3 clusters show ~O(10%) effects. Without that mechanism, CMB wins.

---

## 4. Why ΛCDM / CMB still dominate

1. **Precision isotropy:** Planck temperature and polarization are consistent with statistical isotropy at high accuracy; global anisotropic expansion is disfavored at extreme odds (Saadeh+2016).
2. **Distance ladder + BAO + CMB concordance structure:** Even with H0 tension, the FLRW background + perturbations framework fits a huge dataset; alternatives must beat that evidence, not only Migkas maps.
3. **Causality / Copernican cost:** A real center beyond the horizon that still organizes our local Hubble flow requires either fine-tuned location or exotic pre-horizon physics; void models already paid a Copernican price and largely lost to kSZ/BAO.
4. **Degenerate explanations:** Bulk flow + selection/scatter can produce apparent H0 dipoles without abandoning FLRW.
5. **No ST-cosmology standard model:** There is no established peer-reviewed pipeline from Sedov–Taylor + beyond-horizon mass → CMB spectra + BAO + light curves. Until built and tested, ΛCDM remains the default.

---

## 5. What a serious paper needs to survive peer review

1. **Explicit metric / dynamical system** (not only Newtonian V=kD): state how ST similarity sits in GR; matching across horizons.
2. **CMB-first consistency section:** compute predicted CMB multipoles for the same (l,b,D,M) that fit Migkas; compare to Planck / Saadeh limits. If incompatible, the paper fails.
3. **Quantitative fit** with uncertainties: posterior on CoE direction and D; evidence ratio vs isotropic H0 and vs bulk-flow model.
4. **Redshift dependence:** CoE anisotropic expansion vs bulk flow diverge with z — predict and test (eROSITA clusters, SN Ia).
5. **Systematics control:** X-ray NH, temperature calibration, Malmquist, hemispheric exposure — engage FLAMINGO-style nulls.
6. **Falsification criteria stated up front** (see Status doc Phase E).
7. **No oversell of dark-energy elimination** until the model fits SN+BAO+CMB jointly.

---

## 6. Key citations (journal / year)

1. Migkas et al., A&A 636, A15 (2020) — LX–T isotropy probe.
2. Migkas et al., A&A 649, A151 (2021) — multi-relation H0 dipole / bulk flow.
3. Saadeh, Feeney, Pontzen, Peiris, McEwen, Phys. Rev. Lett. 117, 131302 (2016) — Bianchi / anisotropic expansion limits.
4. Planck Collaboration, A&A 641, A7 (2020) — Planck 2018 isotropy & statistics.
5. Alnes & Amarzguioui, Phys. Rev. D (2006) — off-center CMB in spherical inhomogeneity (astro-ph/0607334).
6. Kodama et al., related analytic off-center CMB formulae (2010).
7. Watkins, Feldman & Hudson, MNRAS (2009) — large-scale bulk flow.
8. Kashlinsky et al., ApJ (2008+) — cluster kSZ bulk flow claims (contested).
9. Sedov (1946); Taylor, Proc. R. Soc. A 201, 159 (1950) — blast-wave similarity.
10. FLAMINGO / A&A (2025) scaling-relation isotropy mocks — reassessment of apparent cluster anisotropies.
11. Pandya et al. (2024); Boubel et al. (2025) — related late-time dipole claims (directions nearer Migkas, mixed amplitudes).

---

## 7. Working conclusion for the center-localization program

- **Supportive enough to motivate work:** Migkas 2020/2021 give a concrete sky direction and amplitude scale for a CoE search.
- **Not supportive enough to claim discovery:** systematics, bulk-flow degeneracy, and simulation nulls weaken the hook.
- **Main adversary:** CMB isotropy / off-center multipole bounds. The project’s credibility hinges on facing that adversary with calculations, not rhetoric about a center beyond the horizon.
