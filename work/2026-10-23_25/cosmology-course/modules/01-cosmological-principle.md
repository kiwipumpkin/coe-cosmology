# Module 1 — The cosmological principle

**Audience.** Chris: chemistry undergraduate background, masters-level applied physics. Plain speech. Physical intuition before symbols. Analogies from equilibrium / mean-field chemistry, continuum limits, symmetry breaking, blast waves, kinetic theory vs hydro, Gauss’s law / shell theorem, lab vs center-of-mass (**CM**) frame.

**Stance.** Rigorous. No cheerleading of \(\Lambda\)CDM (Lambda cold dark matter — the standard cosmological model with a cosmological constant \(\Lambda\) plus cold dark matter) and no cheerleading of any Center-of-Expansion (CoE) research program. Same evidentiary standard both ways. The CoE / off-center geometry is the **course spine**: it appears in sections and exercises from the start, not only as a late bolt-on.

### Acronyms (define once; then use freely)

| Acronym | Meaning |
|---------|---------|
| **CP** | Cosmological principle |
| **CoE** | Center of Expansion (the preferred center in the project-spine geometry) |
| **FLRW** | Friedmann–Lemaître–Robertson–Walker (the homogeneous + isotropic metric family) |
| **CMB** | Cosmic microwave background |
| **LTB** | Lemaître–Tolman–Bondi (spherically symmetric dust models; allow a center) |
| **kSZ** | Kinematic Sunyaev–Zeldovich effect |
| **SN** / **SNIa** | Supernova / Type Ia supernova |
| **GR** | General relativity |
| **BAO** | Baryon acoustic oscillations |
| **CM** | Center of mass (lab vs CM frame) |
| **HI** | Neutral atomic hydrogen |

**Spine symbols (learn now; reuse always).** \(D\): CoE → us. \(r\): vision shell (sources from the CoE). \(R_{\mathrm{CMB}}\): radius of the last-scattering sphere (CMB surface) in the CoE-centered picture; working sketch \(R_{\mathrm{CMB}}\gg D\). Density: **decreasing** outward (not a uniform ball). Gauss / shell: \(M(<r)\) sets the interior field.

---

## 1. What this module is for

By the end you should be able to:

- State the cosmological principle (CP) without mixing it with Copernican slogan-language.
- Separate **homogeneity** from **isotropy**, and **isotropy about us** from **homogeneity of space** — with the off-center CoE sketch as the clarifying example.
- Explain why “on large scales” is both necessary and slippery.
- See the Friedmann–Lemaître–Robertson–Walker (**FLRW**) metric as a **symmetry ansatz**, and name what a preferred center would break.
- List what late-universe isotropy tests actually constrain (clusters, Type Ia supernovae (**SNIa**), radio/quasar dipoles, cosmic microwave background (**CMB**)) — and what they do not.
- Carry \(D\), \(r\), \(R_{\mathrm{CMB}}\), decreasing \(\rho\), and Gauss/shell into later modules without re-deriving them from zero.

---

## 2. Copernican principle vs cosmological principle

### 2.1 Copernican language

The **Copernican principle**, in everyday talk, says: we are not in a special place. Historically it was a reaction against putting Earth at the center of everything. As a *philosophical* preference it is healthy. As a *physics assumption* it is underspecified: “not special” does not tell you the metric, the matter distribution, or the scale at which averaging is allowed.

### 2.2 Cosmological principle as physics

The **cosmological principle**, as used to build standard cosmology, is sharper. It asserts that on sufficiently large scales the universe is **spatially homogeneous** and **isotropic**. In practice that means: there exists a foliation into spatial slices on which the geometry and the mean matter distribution look the same at every point and in every direction.

That is a **symmetry assumption**. It is the input that selects the Friedmann–Lemaître–Robertson–Walker (FLRW) family of metrics (up to the usual curvature parameter and scale factor). It is not deduced from the mere fact that the CMB looks smooth on our sky.

### 2.3 Chris’s useful split (teach this explicitly)

Keep two sentences separate:

1. **Physics assumption:** on large scales, space is homogeneous and isotropic (CP → FLRW symmetry).
2. **Rhetoric / philosophy:** “there is no center” / “we cannot be special.”

Sentence (1) can be tested, pressurized, or replaced by a different ansatz. Sentence (2) is not a substitute for (1).

**Project thread (woven).** A spherical decreasing-density CoE with an observer at distance \(D\neq 0\) from the center is a concrete geometry where (1) fails while (2) becomes a slogan under pressure. The sketch does **not** prove a center exists; it shows that “no center” talk is not what FLRW’s math actually encodes. You can reject Copernican rhetoric and still ask careful questions about homogeneity; you can accept Copernican rhetoric and still misuse FLRW. This course trains (1), using the CoE geometry as the running counter-example.

**Chemistry analogy.** “Molecules are not special” is a worldview. “In the thermodynamic limit the free energy density is intensive and translationally invariant” is a physics assumption you can fail near a phase boundary or in a finite, strongly inhomogeneous sample. CP is closer to the second sentence. Off-center in a Sedov–Taylor blast is closer to “my thermometer is not at the explosion center.”

---

## 3. Homogeneity vs isotropy

| Word | Meaning (physics) | Sky / lab intuition |
|------|-------------------|---------------------|
| **Isotropy about a point** | No preferred direction as seen from that point | Your sky looks the same every way you turn |
| **Homogeneity** | No preferred *location*; every point is equivalent | Translational invariance of the mean fields |

**Important theorem-shaped fact.** If a space is isotropic about *every* point, it is homogeneous. The converse needs care: a homogeneous space need not be isotropic (crystal with preferred axis; Bianchi models). FLRW assumes both.

**What observations usually give you.** A single observer (us) can directly test **isotropy about us**. Homogeneity is harder: you need either (a) a reconstruction that turns redshift surveys into 3D maps under assumptions, or (b) a theoretical argument that isotropy about many observers implies homogeneity, or (c) consistency tests across independent probes. CMB isotropy about us is strong; it is not automatically a proof of homogeneity everywhere.

**Project thread (woven).** Sit at \(D>0\) from a CoE. Your sky can still look *nearly* isotropic if \(D/R_{\mathrm{CMB}}\) and \((\Delta H/H)_{\mathrm{eff}}\) are small (Module 6 budget) — so isotropy about us does **not** by itself prove \(D=0\). Homogeneity fails in the sketch by construction (preferred center). That is the pedagogical point: the two claims come apart.

**Chemistry / continuum analogy.** Local isotropy of a fluid at one thermometer location does not prove the whole tank is uniform. A radial temperature gradient can still exist. Homogeneity is the stronger claim.

---

## 4. What “on large scales” means (and why it is slippery)

Everyone agrees the universe is clumpy on small scales: people, planets, galaxies, clusters, filaments, voids. The CP says: pick a large enough averaging scale \(L\), and the coarse-grained density and geometry look FLRW.

Slippery parts:

1. **Which tracer?** Galaxies, clusters, neutral hydrogen (**HI**), quasars, and the CMB last-scattering surface do not automatically share one \(L\).
2. **Which statistic?** Two-point clustering, counts-in-cells, and directional Hubble fits can “look homogeneous” at different scales.
3. **Which epoch?** Structure grows; the scale that is “linear enough” changes with time.
4. **Operational circularity risk.** Distance and volume measures often assume FLRW to convert redshift to distance. Testing homogeneity with those distances requires care (Module 3).

**Project thread (woven).** Vision shell radius \(r\) (cluster scales, hundreds of Mpc) is not \(R_{\mathrm{CMB}}\) (last scattering). A CoE sketch can produce Migkas-sized contrast on shell \(r\) while the CMB shell at \(R_{\mathrm{CMB}}\gg D\) still demands a tight dipole budget. “Large scales” is probe-dependent inside the same geometry.

**Mean-field analogy.** In solution chemistry, activity coefficients approach 1 in the dilute limit — but “dilute enough” depends on the ion, the solvent, and the observable. “Large scales” in cosmology is the same kind of claim: useful, not automatic.

There is no single peer-reviewed number that is *the* homogeneity scale for all purposes. Treat any quoted scale (hundreds of Mpc, etc.) as probe-dependent and definition-dependent.

---

## 5. FLRW as a symmetry ansatz (and what a center breaks)

Once you assume spatial homogeneity and isotropy (and a reasonable matter model, usually a perfect fluid), Einstein’s equation plus symmetry reduce the metric to FLRW form:

\[
ds^2 = -c^2 dt^2 + a(t)^2 \left[ \frac{dr^2}{1-kr^2} + r^2 d\Omega^2 \right]
\]

(with standard coordinate choices). The scale factor \(a(t)\) carries the expansion history; \(k\) is spatial curvature.

**What the ansatz encodes.**

- No preferred center and no preferred axis *in the geometry of the spatial slices*.
- Worldlines of “fundamental observers” (comoving with the fluid) are interchangeable.
- All large-scale observables that respect the symmetry become functions of time (or redshift) only, not of direction.

**What the ansatz does not encode by itself.**

- Dark energy, inflation, or a particular \(H_0\). Those are extra physics or parameter choices *inside* or *on top of* FLRW.
- Permission to ignore bulk flows, local structure, or lensing when your survey lives in the real clumpy universe.
- A proof that averaging the real universe yields FLRW (Module 8).

**Project thread (woven).** A CoE ansatz with decreasing \(\rho(r)\) measured from a real center **breaks** the “no preferred center” clause. Spherical symmetry about the CoE can remain; FLRW’s *maximal* symmetry (isotropy about *every* point) does not. Module 2 develops Killing / metric-form language; Module 7 meets LTB as the dust cousin that *allows* a center. Flag now: leaving FLRW is allowed in GR; surviving CMB, the kinematic Sunyaev–Zeldovich effect (**kSZ**), supernovae (**SN**), and baryon acoustic oscillations (**BAO**) is the hard part.

**Continuum-limit analogy.** Navier–Stokes assumes a continuum field. It is not a denial that molecules exist; it is a claim that the continuum description is adequate for the question you are asking. FLRW is the continuum description for cosmology. Adequacy is empirical.

**Symmetry-breaking analogy.** Ferromagnets: above \(T_c\) the free energy is isotropic; below \(T_c\) a magnetization direction is chosen. Structure formation breaks FLRW symmetry *locally*. The CP says the *background* remains symmetric when averaged. Whether residual large-scale directional effects survive is exactly what late-universe isotropy tests ask — and what a CoE spine predicts in a controlled way.

---

## 6. What observations actually constrain

Be pedantic on purpose.

### 6.1 Isotropy about us

Maps on the celestial sphere test whether *our* sky is directional. Examples: CMB temperature (after dipole subtraction), number counts, cluster scaling residuals vs Galactic longitude, SN magnitude residuals vs direction, radio source dipoles.

### 6.2 Homogeneity

True homogeneity is a statement about 3D space. You need depth (redshift), a distance model, and control of evolution and selection.

### 6.3 Consistency of FLRW parameters

Measuring the same \(H_0\), \(\Omega_m\), etc., from independent probes tests the FLRW *parameter story*, not always the symmetry itself. Tensions can mean new physics *inside* FLRW, systematics, or a breakdown of the ansatz. Do not leap.

### 6.4 Honest one-liner

**Most “isotropy tests” constrain isotropy about us.** Promoting them to “the universe is homogeneous” always involves extra reasoning. Promoting them to “FLRW is uniquely true” involves still more. Promoting them to “a CoE is established” is a further, usually invalid, leap.

---

## 7. Late-universe isotropy tests (honest list) — project data, not scoreboard

These are tools. They feed the spine; they do not crown a winner.

### 7.1 Galaxy clusters — Migkas et al. (claims with caveats)

**What they did (verified).**

- Migkas et al. (2020), *A&A* **636**, A15: probe cosmic isotropy with X-ray clusters through the \(L_X\)–\(T\) scaling relation. Directional anisotropy reported at \(\gtrsim 4\sigma\) toward Galactic \((l,b)\sim(280^\circ,-20^\circ)\); combining samples (\(\sim 842\) clusters) intensifies toward \(\sim(303^\circ,-27^\circ)\); maximum \(D_L\) variation quoted \(\sim 16\pm 3\%\).
- Migkas et al. (2021), *A&A* **649**, A151 (arXiv:2103.13904): ten scaling relations; apparent \(\sim 9\%\) spatial variation in local \(H_0\) with a nearly dipole form, \(>5\sigma\) in their assessment; alternatively interpretable as a \(\sim 900\,\mathrm{km\,s^{-1}}\) bulk flow out to \(\sim 500\,\mathrm{Mpc}\).

**What \(L_X\)–\(T\) measures.** Temperature is distance-insensitive compared to luminosity. If the normalization \(A\) of \(L_X\)–\(T\) depends on sky direction after careful cleaning, something directional is present in distance, calibration, absorption, or sample — or in the expansion.

**Caveats you must recite.**

- X-ray absorption and calibration systematics can fake directionality.
- The preferred direction need not match the CMB dipole apex.
- Bulk flow vs true directional \(H\) anisotropy can be degenerate (Module 4).

**Project thread (woven).** Digitised Migkas Fig. 6 style points: `/workspace/cosmology/migkas/fig6_digitized.csv`. Research-pipeline \(R^2\), Amp\(\times\sin\), and \(2\ell\) residuals between a toy angular template and those points are **diagnostics**, not GR proofs. A fit quality number does not derive Einstein’s equation, does not fix gauge, and does not settle the CMB constraint. Module 5 deepens this; Module 9 forbids treating \(R^2\) as pass/fail.

### 7.2 Type Ia supernovae

Directional SN fits ask whether magnitude residuals prefer an axis. Powerful but sparse and calibration-sensitive. Always check sample, extinction, and whether the analysis assumes FLRW distances while testing FLRW symmetry.

### 7.3 Radio / quasar dipoles — Secrest et al. (claims with caveats)

**Verified.**

- Secrest et al. (2021), *ApJL* **908**, L51: CatWISE quasars (\(\sim 1.36\times 10^6\)); dipole amplitude \(\sim 0.0155\) toward \((l,b)\approx(238^\circ,29^\circ)\), over twice the kinematic CMB expectation; null rejected at high significance in their analysis (\(p\sim 5\times 10^{-7}\)).
- Secrest et al. (2022), *ApJL* **937**, L31: joint radio (NVSS) + quasar challenge; combined significance quoted \(\sim 5.1\sigma\) against the purely kinematic expectation.

**What they measure.** Anisotropy in the *angular distribution* of distant sources, compared to the dipole expected if we are simply moving relative to the CMB frame. Excess amplitude pressures the “kinematic only” story; systematics, clustering bias, and flux calibration remain live issues. Not the same observable as directional \(H\).

### 7.4 CMB (preview of Module 6 budget)

The CMB is isotropic about us at the \(\sim 10^{-5}\) level in temperature (after removing the \(\sim 10^{-3}\) dipole). The dipole itself is usually interpreted as our peculiar velocity relative to the CMB rest frame (\(\sim 370\,\mathrm{km\,s^{-1}}\)). Alternative **geometric** contributions are tightly limited.

**Project thread (woven) — dipole budget (order of magnitude).** Path asymmetry to a CoE-centered last-scattering sphere \(\sim D/R_{\mathrm{CMB}}\). Sketch:

\[
\left(\frac{\delta T}{T}\right)_{\mathrm{geom}} \approx \alpha \,\frac{D}{R_{\mathrm{CMB}}}\,\left(\frac{\Delta H}{H}\right)_{\mathrm{eff}},\qquad \alpha\sim 0.5\text{–}2.
\]

For \((\Delta H/H)_{\mathrm{eff}}\sim 0.06\) (Migkas-sized) and \(\alpha\sim 1\), keeping \((\delta T/T)_{\mathrm{geom}}\lesssim 1.2\times 10^{-3}\) needs \(D/R_{\mathrm{CMB}}\lesssim 0.02\). Details and figures: Module 6; `/workspace/cosmology/papers/gr_ch2_ltb_2026-09-20/CMB_DIPOLE_ESTIMATE.md`; `/workspace/cosmology/figures/coe_geometry_cmb_dipole_schematic.png`.

### 7.5 Historical cousin

Fennelly (1977), *MNRAS* **181**, 121: modeled reported Hubble anisotropy with large-scale inhomogeneity and an off-center Galaxy. Conceptual ancestor — not a modern constraint.

---

## 8. Decreasing density, Gauss / shell, and why not a uniform ball

**Project spine geometry (full statement).**

Imagine a real **Center of Expansion** and an observer who is **not** at that center. Expansion is homologous in the blast-wave sense (Sedov–Taylor intuition; Bertschinger 1983 is the classic cosmological self-similar shock reference — galaxy-formation context, not a CoE paper). Density is **centrally concentrated and decreases outward**.

**Not a uniform ball.** Uniform density gives enclosed mass \(M(<r)\propto r^3\) and gravitational acceleration \(g\propto -r\) inside (harmonic). That does **not** give a useful directional sky pattern of the Migkas type for an off-center observer in the way a decreasing-\(\rho\) profile can. Teach **Gauss / shell theorem**: exterior shells cancel; interior mass \(M(<r)\) sets the field. Decreasing \(\rho\) makes \(M(<r)\) and local expansion proxies direction-dependent on our sky.

Outward motion in this geometry produces, in the simple sketch, a **heating dipole hot on the anti-CoE side** from geometry + kinematics; **lensing is optional extra**, not required for the sketch (Module 6).

**Alnes–Kodama / kSZ (preview with correct scope).** For **Gpc-scale Lemaître–Tolman–Bondi (LTB) voids with order-unity density contrast**, Alnes & Amarzguioui (2006) found that an off-center observer must lie within \(\sim 15\,\mathrm{Mpc}\) of the void center or the induced CMB dipole exceeds what COBE/WMAP-era data allow. Kodama–Saito–Ishibashi (2010): same class of bound (tens of Mpc, still \(\ll\) Gpc). CMB **does** constrain **spherical off-center + large contrast** very hard. CMB does **not** automatically settle *every* CoE story if \(R_{\mathrm{CMB}}\gg D\) and contrast is only Migkas-sized — but \(D/R_{\mathrm{CMB}}\sim\mathcal{O}(0.1\)–\(0.5)\) still overshoots. The kinematic Sunyaev–Zeldovich effect (kSZ; García-Bellido & Haugbølle 2008) pressures deep Gpc voids even on-center. Modules 6–7 deepen this; Module 9 requires it in the final paper.

**Newtonian scaffolding warning.** Some local “LTB” numerics still evaluate FLRW-shell \(H(r)\) from enclosed mass on an observer-centered source circle — angular templates wearing general-relativity (**GR**) vocabulary. A finished model needs past-directed null geodesics, \(d_L(z,\hat n)\), and a real temperature map. Until then, \(R^2\) vs Migkas is a diagnostic plot, not a discovery. See `/workspace/cosmology/papers/gr_ch2_ltb_2026-09-20/NOTES_FOR_CHRIS.md`.

| Claim | Status in this sketch |
|-------|------------------------|
| “Isotropy about us is observed at high precision in the CMB” | Hard constraint on geometric dipoles |
| “Therefore there is no center anywhere” | Rhetoric, not a deduction |
| “FLRW symmetry is the unique allowed ansatz” | Not settled by one sketch |
| “Migkas directional \(H\) proves a CoE” | False: systematics + bulk-flow degeneracy + CMB budget |
| “Alnes killed every off-center idea” | Overstated: killed *spherical large-contrast void* off-centering at \(\gtrsim 15\,\mathrm{Mpc}\); smaller contrast needs its own accounting |

---

## 9. Chemistry analogies woven tight

| Cosmology idea | Chem / applied-physics cousin |
|----------------|-------------------------------|
| CP / FLRW background | Continuum / mean-field limit |
| Homogeneity | Intensive, translationally invariant free-energy density |
| Isotropy about us | No angular dependence on one lab bench’s instruments |
| “Large scales” | Dilute-limit / coarse-graining length |
| Bulk flow | Advection current vs thermal motion |
| Cluster \(L_X\)–\(T\) directional \(A\) | Calibration curve whose intercept depends on which rack you grab |
| Off-center in a spherical blast | Lab frame offset from explosion center in a Sedov–Taylor problem |
| Gauss / shell theorem \(M(<r)\) | Electrostatics / gravity inside spherical charge / mass distributions |
| Averaging / backreaction | Failure of \(\langle n^2\rangle \approx \langle n\rangle^2\) |

---

## 10. What would pressure or falsify the CP (as used in FLRW)

Falsification in cosmology is rarely a single dagger. Think in terms of **sustained, multi-probe, systematics-controlled directional structure** that cannot be absorbed into peculiar velocity, calibration, or local structure.

Examples of pressure: confirmed late-universe expansion dipole coherent across independent tracers; matter dipole (radio/quasar) irreconcilable with the CMB kinematic frame after exhaustive systematics; directional failure of FLRW consistency relations; reproducible CMB isotropy violation beyond accepted foregrounds.

Examples that do **not** by themselves falsify CP: local bulk flows inside \(\sim 100\,\mathrm{Mpc}\); philosophical preference for a CoE story; a single \(R^2\) between a toy model and one digitised figure.

---

## 11. Exercises (physical intuition + project thread)

**Exercise 1 — Two sentences.**
Write one sentence a physical chemist would accept that defines *isotropy about us*, and one that defines *homogeneity of space*. Then write a third sentence that states a measurement which constrains the first much more directly than the second.

**Exercise 2 — Project thread: Uniform-ball autopsy.**
Using Gauss’s law / shell theorem only, explain why a **uniform-density** sphere is a bad pedagogical CoE if the goal is useful large-scale sky anisotropy in an effective \(g\) or \(H\) for an off-center observer. Contrast with centrally concentrated decreasing \(\rho\). Introduce symbols \(D\), \(r\), \(R_{\mathrm{CMB}}\) in one labeled sketch (words or figure).

**Exercise 3 — Project thread: Rhetoric vs physics + dipole budget.**
(a) Separate “CMB is isotropic about us” from “there is no center anywhere.” (b) Assume \((\delta T/T)_{\mathrm{geom}}\approx(D/R_{\mathrm{CMB}})(\Delta H/H)\) with \(\alpha=1\). For \(\Delta H/H=0.06\) and max geometric dipole \(1.2\times 10^{-3}\), compute max \(D/R_{\mathrm{CMB}}\). If someone proposes \(D/R_{\mathrm{CMB}}=0.3\), by what factor do they overshoot? What Alnes-style lesson survives even when contrast is only \(\sim 6\%\)?

---

## 12. Reading list (verified citations)

### Textbooks / notes

- Barbara Ryden, *Introduction to Cosmology*, 2nd ed. (Cambridge Univ. Press, 2017).
- Steven Weinberg, *Gravitation and Cosmology* (Wiley, 1972); *Cosmology* (Oxford Univ. Press, 2008).
- P. J. E. Peebles, *Principles of Physical Cosmology* (Princeton Univ. Press, 1993).
- David W. Hogg, “Distance measures in cosmology,” arXiv:astro-ph/9905116 (1999; revised 2000).

### Late-universe isotropy / dipoles

- K. Migkas et al., *A&A* **636**, A15 (2020); *A&A* **649**, A151 (2021); arXiv:2103.13904.
- N. J. Secrest et al., *ApJL* **908**, L51 (2021); *ApJL* **937**, L31 (2022).
- A. J. Fennelly, *MNRAS* **181**, 121–130 (1977).

### Off-center spherical models / CMB / kSZ / blast waves

- H. Alnes & M. Amarzguioui, *Phys. Rev. D* **74**, 103520 (2006); arXiv:astro-ph/0607334.
- H. Kodama, K. Saito & A. Ishibashi, arXiv:1004.3089 (2010); K. Saito, A. Ishibashi & H. Kodama, arXiv:1002.3855 (2010).
- J. García-Bellido & T. Haugbølle, *JCAP* **09**, 016 (2008); arXiv:0807.1326.
- E. Bertschinger, *ApJ* **268**, 17–29 (1983).
- J. Smoller & B. Temple, *Commun. Math. Phys.* **210**, 275–308 (2000); arXiv:astro-ph/9812063.

### Averaging (forward reference)

- D. L. Wiltshire, *Phys. Rev. Lett.* **99**, 251101 (2007); *Phys. Rev. D* **80**, 123512 (2009).
- K. Bolejko & M. Korzyński, *Int. J. Mod. Phys. D* **26**, 1730011 (2017); arXiv:1612.08222.

### Local shared-machine resources

- `/workspace/cosmology/figures/coe_geometry_cmb_dipole_schematic.png`
- `/workspace/cosmology/figures/photon_hits_per_arc_vs_angle.png`
- `/workspace/cosmology/papers/gr_ch2_ltb_2026-09-20/CMB_DIPOLE_ESTIMATE.md`
- `/workspace/cosmology/papers/gr_ch2_ltb_2026-09-20/NOTES_FOR_CHRIS.md`
- `/workspace/cosmology/migkas/fig6_digitized.csv`

---

## 13. Takeaway

The cosmological principle is a **symmetry assumption** that buys you FLRW. Isotropy about us is what most maps constrain. Homogeneity is stronger and harder. “No center” talk is not a substitute for either. The decreasing-\(\rho\) off-center CoE sketch — symbols \(D\), \(r\), \(R_{\mathrm{CMB}}\), Gauss/shell, dipole budget — is the **spine** that makes those distinctions concrete across Modules 2–9. Late-universe anisotropy claims (Migkas, Secrest) are serious empirical pressure *if* systematics fail to explain them — and they are not automatic proof of a CoE. Spherical off-center models with large contrast are in deep tension with the CMB dipole (Alnes–Kodama) and with kSZ. Smaller contrast and \(R_{\mathrm{CMB}}\gg D\) change the arithmetic but do not grant a free pass to \(D/R\sim\mathcal{O}(0.1\)–\(0.5)\). Hold those distinctions and you are ready for the rest of the spine.
