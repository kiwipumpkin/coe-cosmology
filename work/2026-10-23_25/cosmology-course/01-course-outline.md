# Course outline — Cosmological principle via CoE / off-center spine

Audience: Chris (undergrad chemistry → masters applied physics). Plain speech, physical intuition, chem/applied-physics bridges. Peer-reviewed and textbook sources preferred — use cites already verified in this pack; do not invent citations.

### Acronyms (define before use; modules repeat)

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
| **ΛCDM** | Lambda cold dark matter |

**Pedagogy.** The **Center-of-Expansion (CoE)** / off-center project (spherical decreasing-density expansion; observer at \(D\); vision shell \(r\); \(R_{\mathrm{CMB}}\gg D\); Gauss/shell; dipole budget; Migkas/Secrest as claims-with-caveats; Alnes–Kodama; **kinematic Sunyaev–Zeldovich (kSZ)**) is the **spine through every module**, not a Module 9 bolt-on. No cheerleading of **ΛCDM** (Lambda cold dark matter) or of the CoE sketch — same evidentiary standard both ways. Amp\(\times\sin\) / \(2\ell\) / \(R^2\) to digitised Migkas are **diagnostics**, not proofs. Expand **cosmological principle (CP)**, **Friedmann–Lemaître–Robertson–Walker (FLRW)**, **cosmic microwave background (CMB)**, **Lemaître–Tolman–Bondi (LTB)**, and **supernova (SN)** on first use in each module draft.

---

## Module 1 — The cosmological principle

**Goal.** State what the CP asserts and does not assert; separate Copernican “no center” rhetoric from homogeneity/isotropy as geometric assumptions; introduce the decreasing-\(\rho\) off-center CoE sketch as the clarifying example that runs the whole course.

**Chem / applied-physics bridge.** Mean-field chemistry and continuum hydro: you replace a molecular mess with a smooth field once a scale exists where fluctuations average. CP is the cosmological version of that continuum limit — useful when it works, dangerous when the averaging scale is unclear. Off-center in a Sedov–Taylor-like blast is the lab-frame offset from the explosion center.

**Readings.**
- Barbara Ryden, *Introduction to Cosmology* (2nd ed., 2017), early chapters on homogeneity/isotropy and the Robertson–Walker metric.
- Steven Weinberg, *Gravitation and Cosmology* (1972), Ch. 14–15; *Cosmology* (2008), opening chapters.
- P. J. E. Peebles, *Principles of Physical Cosmology* (1993), cosmological principle and large-scale structure.
- Migkas et al., *A&A* **636**, A15 (2020); Migkas et al., *A&A* **649**, A151 (2021); arXiv:2103.13904.
- Secrest et al., *ApJL* **908**, L51 (2021).
- A. J. Fennelly, *MNRAS* **181**, 121 (1977) — historical off-center / Hubble-anisotropy example (**project thread**).
- E. Bertschinger, *ApJ* **268**, 17 (1983) — Sedov → cosmological self-similar blast (background for decreasing-\(\rho\) intuition).

**Exercises.**
1. Write two sentences a chemist would accept: one for *isotropy about us*, one for *homogeneity of space*. Where do they diverge?
2. **Project thread: Uniform-ball autopsy.** Using Gauss / shell theorem only, explain why a uniform-density sphere is a bad pedagogical CoE for directional sky structure; contrast with centrally concentrated decreasing \(\rho\). Introduce symbols \(D\), \(r\), \(R_{\mathrm{CMB}}\).
3. **Project thread: Rhetoric vs physics.** “There is no center anywhere” vs “spatial sections are homogeneous and isotropic.” Which is data, which is add-on? How does an off-center decreasing-\(\rho\) sketch force the split?

**Draft:** `modules/01-cosmological-principle.md`

---

## Module 2 — Symmetry and the FLRW metric

**Goal.** See FLRW as what you get when you demand maximal spatial symmetry on each constant-time slice; state what the metric form / Killing structure **forbids** (a preferred center); state what a CoE ansatz **breaks**.

**Chem / applied-physics bridge.** Continuum limit of kinetic theory → Navier–Stokes: once you assume local equilibrium and spatial uniformity of the “background,” the equations collapse. FLRW is that collapse for gravity + fluid. A spherical blast with a real center is the symmetry-breaking cousin.

**Readings.**
- Weinberg, *Gravitation and Cosmology* (1972), derivation of the Robertson–Walker metric from homogeneity and isotropy.
- Ryden, *Introduction to Cosmology*, Robertson–Walker metric chapter.
- Optional depth: Killing vectors and maximal symmetry (Weinberg is enough).
- **Project thread:** revisit Alnes & Amarzguioui, *Phys. Rev. D* **74**, 103520 (2006) abstract/intro — spherical symmetry *with* a center is a different ansatz class than FLRW.

**Exercises.**
1. If space is isotropic about *every* point, why must it also be homogeneous? (One-paragraph physical argument.)
2. **Project thread: What FLRW forbids.** List two things FLRW cannot describe without leaving the ansatz (preferred center; bulk flow with a preferred axis). For each, say which symbol in the CoE sketch (\(D\), preferred axis toward/away from CoE) marks the break.
3. Analogy: what continuum assumption in fluid dynamics is closest to “perfect fluid on FLRW slices”? When does a blast-wave description refuse that assumption?

**Draft:** `modules/02-flrw-symmetry.md`

---

## Module 3 — Distance measures and redshift

**Goal.** Master Hogg’s distance ladder language (comoving, luminosity, angular-diameter) and see how an off-center observer on a CoE-centered vision shell confuses directional distance and inferred \(H\).

**Chem / applied-physics bridge.** Lab frame vs center-of-mass frame: the same collision looks different once you boost. Cosmological redshift mixes expansion, peculiar velocity, and path-dependent gravity; frame choice matters — especially when you are not at the symmetry center.

**Readings.**
- David W. Hogg, “Distance measures in cosmology,” arXiv:astro-ph/9905116 (1999; rev. 2000).
- Ryden, distance and redshift chapters.
- Weinberg, *Cosmology* (2008), luminosity distance and related measures.
- **Project thread:** course figures `/workspace/cosmology/figures/coe_geometry_cmb_dipole_schematic.png` (path length \(R\pm D\) intuition carries to distance as well as CMB).

**Exercises.**
1. In words: why can luminosity distance and angular-diameter distance disagree once expansion history is nontrivial?
2. **Project thread: Off-center shell → directional \(H\).** You sit at distance \(D\) from the CoE. Sources on a shell of radius \(r\) from the CoE are at different distances from *you*. Which observable can look like a sky dipole in \(H\)? Sketch with \(r\) and \(D\) labeled.
3. Name one quantity that is frame-dependent (lab-frame velocity) and one closer to an invariant (proper interval along a null geodesic). Map “lab offset from blast center” to \(D\).

**Draft:** `modules/03-distances-redshift.md`

---

## Module 4 — Hubble law locally vs cosmologically; bulk flows

**Goal.** Separate the local linear velocity–distance relation from cosmological expansion history; understand bulk flows as coherent peculiar motion that can mimic directional \(H\) variation — including the \(\sim 900\,\mathrm{km\,s^{-1}}\) Migkas alternative.

**Chem / applied-physics bridge.** Kinetic theory vs hydro: local random motions vs organized flow. A bulk flow is hydro-like; shot noise in galaxy peculiar velocities is kinetic-theory-like. CoE-driven directional expansion is a third story that must be separated from both.

**Readings.**
- Peebles, *Principles of Physical Cosmology* — peculiar velocities and large-scale flows.
- Migkas et al. (2021), *A&A* **649**, A151 — directional \(H_0\) vs bulk-flow interpretation (**project thread**).
- Optional: Heinesen multipole expansion of the Hubble law in general spacetimes (cite from local literature pack if assigning homework).

**Exercises.**
1. When is \(v = H_0 d\) a *definition* of local \(H_0\) and when is it a *prediction*?
2. **Project thread: Bulk flow vs CoE degeneracy.** How could a bulk flow of \(\sim 900\,\mathrm{km\,s^{-1}}\) out to hundreds of Mpc (Migkas 2021 alternative) look like a \(\sim 9\%\) \(H_0\) dipole? What observation would break the degeneracy with a true off-center expansion anisotropy?
3. Chemistry analogy: distinguish equilibrium thermal motion from a driven advection current. Map each to cosmology; say where a CoE sketch sits relative to both.

**Draft:** `modules/04-hubble-bulk-flows.md`

---

## Module 5 — Late-universe anisotropy claims (Migkas / Secrest)

**Goal.** Understand what cluster \(L_X\)–\(T\), SNIa directional fits, and radio/quasar dipoles *measure*, including systematics; treat Amp\(\times\sin\) / \(R^2\) to digitised Migkas as diagnostics only.

**Chem / applied-physics bridge.** Calibration curves in analytical chemistry: if your standard curve’s intercept depends on which rack of samples you use, you have either real sample heterogeneity or a lab systematic. Migkas-style directional \(A\) in \(L_X\)–\(T\) is that problem at cosmic scale.

**Readings.**
- Migkas et al. (2020), *A&A* **636**, A15 — \(L_X\)–\(T\) vs sky direction; digitised Fig. 6 style data at `/workspace/cosmology/migkas/fig6_digitized.csv`.
- Migkas et al. (2021), *A&A* **649**, A151 — \(\sim 9\%\) directional \(H_0\); bulk-flow alternative.
- Secrest et al. (2021), *ApJL* **908**, L51; Secrest et al. (2022), *ApJL* **937**, L31 — quasar + radio dipole.
- Caveats: X-ray absorption / calibration; apex may not match CMB dipole; bulk flow vs true \(H\) degeneracy.

**Exercises.**
1. Why does \(L_X\)–\(T\) probe distance (and thus inferred \(H\)) if temperature is roughly distance-independent?
2. **Project thread: Diagnostics ≠ proofs.** If research notes report an \(R^2\) (or Amp\(\times\sin\) / \(2\ell\) residual) between a toy angular template and Migkas \(A(\ell)\), treat that number as a diagnostic. Explain in two sentences why it cannot establish a GR CoE model. Point to `/workspace/cosmology/migkas/fig6_digitized.csv`.
3. List three non-cosmological effects that could fake a directional cluster anisotropy. For each, say whether it would also fake a Secrest-style number-count dipole.

**Draft:** `modules/05-late-universe-anisotropy.md`

---

## Module 6 — CMB isotropy and the dipole budget

**Goal.** State what CMB isotropy constrains; teach the geometric dipole budget \(D/R_{\mathrm{CMB}}\lesssim 0.02\) for Migkas-sized contrast; separate heating-dipole anti-CoE sketch from lensing; use shared figures.

**Chem / applied-physics bridge.** Blast waves / Sedov–Taylor: spherical explosion with **centrally concentrated, decreasing** density (not a uniform ball). Uniform \(\rho\) gives \(g\propto -r\) inside; decreasing \(\rho\) plus Gauss/shell \(M(<r)\) gives directional structure for an off-center observer. Outward homologous motion → geometric + kinematic heating dipole (hot on the anti-CoE side in the simple sketch).

**Teach symbols (spine).** \(D\): CoE→us. \(r\): vision shell. \(R_{\mathrm{CMB}}\): last-scattering sphere; \(R_{\mathrm{CMB}}\gg D\).

**Geometric dipole sketch.**
\[
(\delta T/T)_{\mathrm{geom}} \approx \alpha \,(D/R_{\mathrm{CMB}})\,(\Delta H/H)_{\mathrm{eff}},
\]
with \(\alpha\sim 0.5\)–\(2\). For \(\Delta H/H\sim 0.06\) and \(\alpha\sim 1\), keeping \((\delta T/T)_{\mathrm{geom}}\lesssim 1.2\times 10^{-3}\) needs \(D/R_{\mathrm{CMB}}\lesssim 0.02\). See `/workspace/cosmology/papers/gr_ch2_ltb_2026-09-20/CMB_DIPOLE_ESTIMATE.md`.

**Figures.** `/workspace/cosmology/figures/coe_geometry_cmb_dipole_schematic.png`; `/workspace/cosmology/figures/photon_hits_per_arc_vs_angle.png`.

**Readings.**
- Alnes & Amarzguioui, *Phys. Rev. D* **74**, 103520 (2006); arXiv:astro-ph/0607334 — \(\sim 15\,\mathrm{Mpc}\) for Gpc-scale order-unity contrast voids.
- Kodama, Saito & Ishibashi, arXiv:1004.3089 (2010) — tens of Mpc in some setups; still \(\ll\) Gpc.
- García-Bellido & Haugbølle, *JCAP* **09**, 016 (2008) — kSZ (forward to Module 7; flag now).
- Bertschinger, *ApJ* **268**, 17 (1983).

**Exercises.**
1. Why does a **uniform-density** ball fail as a pedagogical CoE for sky anisotropy? (Gauss / shell theorem.)
2. **Project thread: Dipole budget arithmetic.** For \(\Delta H/H=0.06\), \(\alpha=1\), max geom dipole \(1.2\times 10^{-3}\), compute max \(D/R_{\mathrm{CMB}}\). If someone proposes \(D/R_{\mathrm{CMB}}=0.3\), by what factor do they overshoot? How does allowed \(D\) compare to classic Alnes \(\sim 15\,\mathrm{Mpc}\) when contrast is only Migkas-sized?
3. Separate: “CMB is isotropic about us” vs “there is no center anywhere.” Which is closer to data?

**Draft:** `modules/06-cmb-dipole-budget.md`

---

## Module 7 — Inhomogeneous models (LTB) and Alnes–Kodama / kSZ

**Goal.** Meet LTB as the spherically symmetric dust generalization of FLRW; understand Alnes/Kodama-class CMB bounds as honest limits on *spherical off-center + large contrast*; state what those bounds do **not** settle; include kSZ.

**Chem / applied-physics bridge.** Spherical blast / Sedov–Taylor and Gauss’s law problems: high symmetry buys analytic control and also buys strong observational vetoes once you leave the center.

**Readings.**
- Alnes & Amarzguioui (2006); Kodama et al. (2010); Saito, Ishibashi & Kodama, arXiv:1002.3855 (2010).
- Smoller & Temple, *Commun. Math. Phys.* **210**, 275 (2000); arXiv:astro-ph/9812063 — FRW with spherical shock cutoff (ontology example).
- García-Bellido & Haugbølle (2008) — kSZ.
- Local honesty notes: `/workspace/cosmology/papers/gr_ch2_ltb_2026-09-20/NOTES_FOR_CHRIS.md` (Newtonian scaffolding caveats).

**Exercises.**
1. List three observables that pressure Gpc-scale voids even if SNIa look acceptable on-center.
2. **Project thread: Scope sentence.** Rewrite “CMB constrains spherical off-center + large contrast hard” so a lab partner cannot confuse it with “every CoE idea is dead.” Include Alnes \(\sim 15\,\mathrm{Mpc}\) vs Migkas-sized contrast arithmetic.
3. **Project thread: kSZ + dipole gate.** What would you need (qualitatively) to keep \(R_{\mathrm{CMB}}\gg D\) *and* Migkas-sized late anisotropy *without* overshooting the CMB dipole *or* implying a deep void ruled out by kSZ?

**Draft:** `modules/07-ltb-alnes-ksz.md`

---

## Module 8 — Averaging, backreaction, mean-field failure

**Goal.** Understand Buchert averaging / backreaction language and Wiltshire–Bolejko-type arguments for when a smooth FLRW mean field can misrepresent the real volume average; ask whether a CoE sketch needs backreaction language or **explicit** inhomogeneity.

**Chem / applied-physics bridge.** Mean-field theory in statistical mechanics: replacing \(\langle n^2\rangle\) by \(\langle n\rangle^2\) fails near critical points and in strongly inhomogeneous media. Backreaction asks whether that replacement is safe for gravity. A CoE with a real center is usually *explicit* inhomogeneity — different from “averaged FLRW plus small correction.”

**Readings.**
- Wiltshire, *Phys. Rev. Lett.* **99**, 251101 (2007); Wiltshire, *Phys. Rev. D* **80**, 123512 (2009).
- Bolejko & Korzyński, *Int. J. Mod. Phys. D* **26**, 1730011 (2017); arXiv:1612.08222.
- Optional: Clarkson, Bassett & Lu inhomogeneity / Copernican tests (assign from primary literature when needed).
- **Project thread:** revisit Module 1–7 geometry — is CoE a backreaction story or an LTB-class explicit model?

**Exercises.**
1. Give a chemistry example where a mean-field rate law fails, and map it to “light travels through voids and walls differently.”
2. What observation would convince you backreaction is *small* on Hubble scales? What would convince you it is *not*?
3. **Project thread: Backreaction vs explicit CoE.** Why is “we averaged the density” not the same as “we averaged the metric”? Does the CoE spine (preferred center, \(D\), decreasing \(\rho\)) live naturally in Buchert language, or does it demand an explicit inhomogeneous metric (LTB-class)? Defend in one paragraph.

**Draft:** `modules/08-averaging-backreaction.md`

---

## Module 9 — Spine review + final paper (one falsifiable claim)

**Goal.** **Review of the entire spine** — revisit every project thread from Modules 1–8; culminate in one final paper defending a **single falsifiable claim** about the CoE/off-center model (observable + dataset + pass/fail; Alnes–Kodama + kSZ stated plainly). This is **synthesis of the threaded course**, not “design a test from scratch as a standalone capstone.”

**Chem / applied-physics bridge.** Experimental design / pre-registration: fix kill criteria before tuning diagnostics; the notebook of prior labs (Modules 1–8) is what you synthesize, not a new unrelated experiment.

**Readings.**
- Revisit project threads and cites from Modules 1–8: Migkas (2020, 2021), Secrest (2021, 2022), Alnes (2006), Kodama (2010), García-Bellido & Haugbølle (2008), Hogg (1999), Wiltshire/Bolejko as needed.
- Local notes: `/workspace/cosmology/papers/gr_ch2_ltb_2026-09-20/NOTES_FOR_CHRIS.md`; `CMB_DIPOLE_ESTIMATE.md`.

**Structure of the module draft.**
1. Checklist revisiting project threads 1–8.
2. Final-paper prompt with required structure.
3. Worked example final paper / memo.
4. Rubric.

**Draft:** `modules/09-spine-review-final-paper.md` (replaces old `09-falsifiable-coe-test.md`)
