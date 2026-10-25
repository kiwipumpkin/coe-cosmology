# Module 2 — Symmetry and the FLRW metric

**Audience.** Chris: chemistry → applied physics. Plain speech; intuition before Killing-vector algebra.

**Stance.** No cheerleading of ΛCDM (Lambda cold dark matter) or of the Center-of-Expansion (CoE) sketch. Same evidentiary standard both ways.

### Acronyms

| Acronym | Meaning |
|---------|---------|
| **CP** | Cosmological principle |
| **CoE** | Center of Expansion |
| **FLRW** | Friedmann–Lemaître–Robertson–Walker |
| **CMB** | Cosmic microwave background |
| **LTB** | Lemaître–Tolman–Bondi |
| **GR** | General relativity |

**Spine reminder.** \(D\) = CoE→us; \(r\) = vision shell from CoE; \(R_{\mathrm{CMB}}\) = last-scattering / cosmic microwave background (**CMB**) shell; decreasing \(\rho\) (not uniform ball); Gauss/shell.

---

## 1. Goal

See the **Friedmann–Lemaître–Robertson–Walker (FLRW)** metric as what you get when you demand maximal spatial symmetry (homogeneous + isotropic) on each constant-time slice. State what the metric form / Killing structure **forbids** (a preferred center). State what a **Center-of-Expansion (CoE)** ansatz **breaks**.

---

## 2. Chem / applied-physics bridge

Continuum limit of kinetic theory → Navier–Stokes: once you assume local equilibrium and spatial uniformity of the “background,” the equations collapse. FLRW is that collapse for gravity + fluid. A spherical blast with a real center is the symmetry-breaking cousin — still continuum physics, different symmetry.

---

## 3. What “maximal spatial symmetry” buys you

The **cosmological principle (CP)** asserts homogeneity + isotropy of spatial sections. In **general relativity (GR)**, that symmetry is implemented by requiring enough **Killing vectors** (infinitesimal isometries) that every point and every direction on a spatial slice looks the same.

Under those assumptions (plus a perfect-fluid matter model), Einstein’s equation reduces the geometry to FLRW form:

\[
ds^2 = -c^2 dt^2 + a(t)^2 \left[ \frac{dr^2}{1-kr^2} + r^2 d\Omega^2 \right].
\]

Here \(a(t)\) is the scale factor; \(k\) is spatial curvature. Coordinates can be chosen so that “fundamental observers” are at fixed spatial coordinates.

**What this forbids in the geometry of the slices:**

- A preferred **center** (no point is special).
- A preferred **axis** (no direction is special).
- Large-scale observables that depend on sky direction *because of the background metric itself*.

**What this does not forbid:**

- Local structure, peculiar velocities, and bulk flows as *perturbations* or as motions relative to the FLRW frame.
- Choosing coordinates that *look* centered on us for convenience — that is gauge, not a physical CoE.

---

## 4. Project thread — what a CoE ansatz breaks

A CoE sketch with decreasing density \(\rho(r)\) measured from a real center keeps **spherical symmetry about the CoE** but drops **isotropy about every point**. The observer at distance \(D\neq 0\) sees a preferred axis (toward / away from the CoE).

| Feature | FLRW | CoE + off-center observer |
|---------|------|---------------------------|
| Preferred center | Forbidden | Built in |
| Isotropy about us | Exact (background) | Approximate if \(D/R_{\mathrm{CMB}}\) and contrast small |
| Homogeneity | Exact (background) | Broken by construction |
| Natural cousin in GR | FLRW family | **Lemaître–Tolman–Bondi (LTB)** and related spherical models (Module 7) |

**Alnes & Amarzguioui (2006)** study CMB anisotropies for an off-center observer in a spherically symmetric inhomogeneous universe — that is the LTB-class world, not FLRW. Reading their setup is enough for this module: different ansatz, different allowed observables.

Leaving FLRW is allowed in GR. Surviving the CMB dipole budget, kSZ, supernovae, and baryon acoustic oscillations is the hard empirical question — not settled by writing a non-FLRW metric.

---

## 5. Readings

- Steven Weinberg, *Gravitation and Cosmology* (1972), derivation of the Robertson–Walker metric from homogeneity and isotropy.
- Barbara Ryden, *Introduction to Cosmology* (2nd ed., 2017), Robertson–Walker chapter.
- Optional: Killing vectors / maximal symmetry (Weinberg is enough).
- **Project thread:** H. Alnes & M. Amarzguioui, *Phys. Rev. D* **74**, 103520 (2006); arXiv:astro-ph/0607334 — abstract/intro: spherical symmetry *with* a center ≠ FLRW.

---

## 6. Exercises

1. If space is isotropic about *every* point, why must it also be homogeneous? (One-paragraph physical argument.)

2. **Project thread: What FLRW forbids.** List two things FLRW cannot describe without leaving the ansatz (preferred center; bulk flow with a preferred axis). For each, say which symbol in the CoE sketch (\(D\), preferred axis toward/away from CoE) marks the break.

3. Analogy: what continuum assumption in fluid dynamics is closest to “perfect fluid on FLRW slices”? When does a blast-wave description refuse that assumption?

---

## 7. Takeaway

FLRW is a **symmetry ansatz**, not a theorem. The CoE spine breaks the “no preferred center” clause on purpose. Module 3 asks what that does to distances and redshift for an off-center observer.
