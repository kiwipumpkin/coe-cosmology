# Module 7 — Inhomogeneous models (LTB) and Alnes–Kodama / kSZ

**Audience.** Chris: chemistry → applied physics.

**Stance.** State bounds with correct scope. No cheerleading of ΛCDM (Lambda cold dark matter) or of the Center-of-Expansion (CoE) sketch.

### Acronyms

| Acronym | Meaning |
|---------|---------|
| **CP** | Cosmological principle |
| **CoE** | Center of Expansion |
| **FLRW** | Friedmann–Lemaître–Robertson–Walker |
| **CMB** | Cosmic microwave background |
| **LTB** | Lemaître–Tolman–Bondi |
| **kSZ** | Kinematic Sunyaev–Zeldovich effect |
| **SN** / **SNIa** | Supernova / Type Ia supernova |
| **GR** | General relativity |
| **BAO** | Baryon acoustic oscillations |

**Spine.** \(D\), \(r\), \(R_{\mathrm{CMB}}\gg D\); decreasing \(\rho\); Gauss/shell; dipole budget from Module 6.

---

## 1. Goal

Meet **Lemaître–Tolman–Bondi (LTB)** as the spherically symmetric dust generalization of **Friedmann–Lemaître–Robertson–Walker (FLRW)**. Understand Alnes/Kodama-class **cosmic microwave background (CMB)** bounds as honest limits on *spherical off-center + large contrast*. Say what those bounds do **not** settle. Include the **kinematic Sunyaev–Zeldovich (kSZ)** effect.

---

## 2. Chem / applied-physics bridge

Spherical blast / Sedov–Taylor and Gauss’s law problems: high symmetry buys analytic control and also buys strong observational vetoes once you leave the center.

---

## 3. LTB in one page

FLRW assumes homogeneity + isotropy (the **cosmological principle (CP)**). LTB keeps spherical symmetry about a center but allows the density, curvature, and bang time to depend on radial coordinate. An observer at \(D\neq 0\) from that center is the Alnes-class setup.

**Smoller & Temple (2000),** *Commun. Math. Phys.* **210**, 275: FRW with a spherical shock cutoff — ontology example of a center-compatible relativistic continuum, not a finished observational CoE model.

**Honesty from local notes.** Some research-pipeline “LTB” numerics still evaluate FLRW-shell \(H(r)\) from enclosed mass on an observer-centered circle — Newtonian scaffolding wearing **general relativity (GR)** clothes. A finished model needs past-directed null geodesics, \(d_L(z,\hat n)\), and a real temperature map. See `/workspace/cosmology/papers/gr_ch2_ltb_2026-09-20/NOTES_FOR_CHRIS.md`.

---

## 4. Alnes–Kodama — what the bounds do and do not settle

**Alnes & Amarzguioui (2006),** *Phys. Rev. D* **74**, 103520: for **Gpc-scale** voids with **order-unity** density contrast, an off-center observer must lie within \(\sim 15\,\mathrm{Mpc}\) of the void center or the induced CMB dipole exceeds what COBE/WMAP-era data allow.

**Kodama, Saito & Ishibashi (2010),** arXiv:1004.3089 (and Saito, Ishibashi & Kodama, arXiv:1002.3855): analytic formulae for off-center CMB anisotropy; bounds soften to **tens of Mpc** in some setups — still \(\ll\) Gpc.

| Correct scope | Overstatement to avoid |
|---------------|------------------------|
| CMB constrains **spherical off-center + large contrast** very hard | “Every CoE idea is dead” |
| Migkas-sized contrast (\(\sim 6\)–\(9\%\)) changes the arithmetic (Module 6: \(D/R_{\mathrm{CMB}}\lesssim 0.02\)) | “Alnes is irrelevant if contrast is 6%” — the *lesson* (dipole \(\propto (D/R)\times\)contrast) still kills large \(D/R\) |
| Bounds are for the model class studied | Different geometry / cancellation mechanisms need their own full-sky \(\delta T/T\) calculation — not assumed |

---

## 5. kSZ — García-Bellido & Haugbølle (2008)

**García-Bellido & Haugbølle (2008),** *JCAP* **09**, 016: even an **on-center** observer sees off-center clusters; those clusters see a large CMB dipole that appears to us as a **kSZ** signal. Deep Gpc voids with significant underdensity are strongly pressured by kSZ — a parallel kill channel to the off-center CMB dipole.

Project implication: a decreasing-\(\rho\) CoE that implies a *deep* Gpc underdensity must face kSZ limits even if \(D\) is small. Shallow (Migkas-sized) contrast is a different regime — still not a free pass without a profile↔kSZ map.

Other pressures on Gpc voids (even if Type Ia supernovae (**SNIa**) look acceptable on-center): CMB, kSZ, baryon acoustic oscillations (**BAO**) / large-scale structure consistency — list qualitatively in exercises.

---

## 6. Readings

- Alnes & Amarzguioui (2006); Kodama et al. (2010); Saito, Ishibashi & Kodama, arXiv:1002.3855 (2010).
- Smoller & Temple (2000); arXiv:astro-ph/9812063.
- García-Bellido & Haugbølle (2008).
- `/workspace/cosmology/papers/gr_ch2_ltb_2026-09-20/NOTES_FOR_CHRIS.md`.

---

## 7. Exercises

1. List three observables that pressure Gpc-scale voids even if SNIa look acceptable on-center.

2. **Project thread: Scope sentence.** Rewrite “CMB constrains spherical off-center + large contrast hard” so a lab partner cannot confuse it with “every CoE idea is dead.” Include Alnes \(\sim 15\,\mathrm{Mpc}\) vs Migkas-sized contrast arithmetic.

3. **Project thread: kSZ + dipole gate.** What would you need (qualitatively) to keep \(R_{\mathrm{CMB}}\gg D\) *and* Migkas-sized late anisotropy *without* overshooting the CMB dipole *or* implying a deep void ruled out by kSZ?

---

## 8. Takeaway

LTB allows a center; Alnes–Kodama and kSZ police large-contrast spherical voids hard. Module 8 asks whether “backreaction” language is the right home for a CoE spine or whether explicit inhomogeneity is required.
