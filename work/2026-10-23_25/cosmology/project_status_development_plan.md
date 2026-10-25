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
