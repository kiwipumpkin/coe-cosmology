# Module 5 — Late-universe anisotropy claims (Migkas / Secrest)

**Audience.** Chris: chemistry → applied physics.

**Stance.** Claims-with-caveats. No cheerleading of ΛCDM (Lambda cold dark matter) or of the Center-of-Expansion (CoE) sketch. Amp\(\times\sin\) / \(2\ell\) / \(R^2\) to digitised Migkas are **diagnostics**, not proofs.

### Acronyms

| Acronym | Meaning |
|---------|---------|
| **CP** | Cosmological principle |
| **CoE** | Center of Expansion |
| **FLRW** | Friedmann–Lemaître–Robertson–Walker |
| **CMB** | Cosmic microwave background |
| **SN** / **SNIa** | Supernova / Type Ia supernova |
| **GR** | General relativity |

**Spine.** \(D\), \(r\), \(R_{\mathrm{CMB}}\); decreasing \(\rho\); Gauss/shell; bulk-flow degeneracy from Module 4.

---

## 1. Goal

Understand what cluster \(L_X\)–\(T\), Type Ia supernova (**SNIa**) directional fits, and radio/quasar dipoles *measure*, including systematics. Separate “anisotropy of the expansion” from “anisotropy of the tracer.”

---

## 2. Chem / applied-physics bridge

Calibration curves in analytical chemistry: if your standard curve’s intercept depends on which rack of samples you use, you have either real sample heterogeneity or a lab systematic. Migkas-style directional \(A\) in \(L_X\)–\(T\) is that problem at cosmic scale.

---

## 3. Migkas et al. — what was measured (verified)

- **Migkas et al. (2020),** *A&A* **636**, A15: X-ray galaxy clusters through the \(L_X\)–\(T\) scaling relation; directional anisotropy reported at \(\gtrsim 4\sigma\); maximum \(D_L\) variation \(\sim 16\pm 3\%\) in their assessment.
- **Migkas et al. (2021),** *A&A* **649**, A151: ten scaling relations; apparent \(\sim 9\%\) spatial variation in local \(H_0\); alternatively \(\sim 900\,\mathrm{km\,s^{-1}}\) bulk flow out to \(\sim 500\,\mathrm{Mpc}\).

**Why \(L_X\)–\(T\) cares about distance.** Temperature is roughly distance-insensitive compared to X-ray luminosity. If normalization \(A\) depends on sky direction after cleaning, something directional sits in distance, calibration, absorption, selection — or expansion.

**Caveats (recite).** Galactic absorption maps; X-ray temperature calibration; Malmquist / flux-limit anisotropy; bulk flow vs true \(H\) degeneracy; preferred axis need not match the cosmic microwave background (**CMB**) dipole apex.

Digitised Fig. 6 style points: `/workspace/cosmology/migkas/fig6_digitized.csv`.

---

## 4. Secrest et al. — what was measured (verified)

- **Secrest et al. (2021),** *ApJL* **908**, L51: CatWISE quasars; dipole amplitude \(\sim 0.0155\), over twice the kinematic CMB expectation in their analysis.
- **Secrest et al. (2022),** *ApJL* **937**, L31: joint radio (NVSS) + quasar challenge; combined significance \(\sim 5.1\sigma\) against purely kinematic expectation in their assessment.

**What they measure:** anisotropy in *number counts / angular distribution*, not the same as directional \(H\). Systematics, clustering bias, and flux calibration remain live. Claims-with-caveats — not automatic proof of a CoE, not automatic dismissal either.

---

## 5. Project thread — diagnostics ≠ proofs

Research-pipeline fits of a toy angular template (Amp\(\times\sin\), \(2\ell\) residuals, \(R^2\)) to digitised Migkas \(A(\ell)\) are **diagnostics**:

- They do not derive Einstein’s equation or fix gauge in **general relativity (GR)**.
- They do not satisfy the CMB dipole budget by themselves.
- They have no survey selection function attached in this course pack.

A CoE hypothesis that “fits Migkas with \(R^2\approx 0.5\)” has done homework, not discovery. Module 9 forbids using \(R^2\) as pass/fail.

SNIa directional residuals are a parallel probe: powerful distances, sparse sky, extinction-sensitive; check whether the analysis assumes **Friedmann–Lemaître–Robertson–Walker (FLRW)** distances while testing FLRW symmetry (from the **cosmological principle (CP)**).

---

## 6. Readings

- Migkas et al., *A&A* **636**, A15 (2020); *A&A* **649**, A151 (2021).
- Secrest et al., *ApJL* **908**, L51 (2021); *ApJL* **937**, L31 (2022).
- Local CSV: `/workspace/cosmology/migkas/fig6_digitized.csv`.

---

## 7. Exercises

1. Why does \(L_X\)–\(T\) probe distance (and thus inferred \(H\)) if temperature is roughly distance-independent?

2. **Project thread: Diagnostics ≠ proofs.** If research notes report an \(R^2\) (or Amp\(\times\sin\) / \(2\ell\) residual) between a toy angular template and Migkas \(A(\ell)\), treat that number as a diagnostic. Explain in two sentences why it cannot establish a GR CoE model. Point to `/workspace/cosmology/migkas/fig6_digitized.csv`.

3. List three non-cosmological effects that could fake a directional cluster anisotropy. For each, say whether it would also fake a Secrest-style number-count dipole.

---

## 8. Takeaway

Migkas and Secrest are serious empirical pressure *if* systematics fail — and they are not automatic CoE proof. Module 6 budgets the CMB geometric dipole against the same \((D,(\Delta H/H)_{\mathrm{eff}})\).
