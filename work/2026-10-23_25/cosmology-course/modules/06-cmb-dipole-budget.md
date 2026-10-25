# Module 6 — CMB isotropy and the dipole budget

**Audience.** Chris: chemistry → applied physics.

**Stance.** No cheerleading of ΛCDM (Lambda cold dark matter) or of the Center-of-Expansion (CoE) sketch. Same evidentiary standard both ways.

### Acronyms

| Acronym | Meaning |
|---------|---------|
| **CP** | Cosmological principle |
| **CoE** | Center of Expansion |
| **FLRW** | Friedmann–Lemaître–Robertson–Walker |
| **CMB** | Cosmic microwave background |
| **LTB** | Lemaître–Tolman–Bondi |
| **kSZ** | Kinematic Sunyaev–Zeldovich effect |
| **GR** | General relativity |

**Spine symbols.** \(D\) = CoE→us; \(r\) = vision shell; \(R_{\mathrm{CMB}}\) = last-scattering sphere in the CoE-centered picture; working assumption \(R_{\mathrm{CMB}}\gg D\); decreasing \(\rho\) (not uniform ball); Gauss/shell.

---

## 1. Goal

State what **cosmic microwave background (CMB)** isotropy constrains. Teach the geometric dipole budget \(D/R_{\mathrm{CMB}}\lesssim 0.02\) for Migkas-sized contrast. Separate the heating-dipole anti-CoE sketch from optional lensing. Use shared figures.

---

## 2. Chem / applied-physics bridge

Blast waves / Sedov–Taylor: a spherical explosion with **centrally concentrated, decreasing** density (not a uniform ball). Uniform \(\rho\) gives \(g\propto -r\) inside (harmonic oscillator). Decreasing \(\rho\) plus Gauss / shell theorem \(M(<r)\) gives directional structure for an off-center observer. Outward homologous motion → geometric + kinematic **heating dipole hot on the anti-CoE side** in the simple sketch; lensing is optional extra physics.

---

## 3. What CMB isotropy says

The CMB is isotropic about us at the \(\sim 10^{-5}\) level in temperature after removing the \(\sim 10^{-3}\) dipole. The dipole is usually interpreted as our peculiar velocity relative to the CMB rest frame (\(\sim 370\,\mathrm{km\,s^{-1}}\); \(\delta T/T\sim 1.23\times 10^{-3}\)).

That is the strongest isotropy-about-us statement in cosmology. It is **not** automatically the statement “there is no center anywhere” (rhetoric vs the **cosmological principle (CP)** / **Friedmann–Lemaître–Robertson–Walker (FLRW)** physics — Module 1).

---

## 4. Geometric dipole sketch (project spine)

Along the axis, path lengths from an off-center observer to a CoE-centered shell of radius \(R_{\mathrm{CMB}}\) differ: roughly \(R+D\) toward the CoE and \(R-D\) away. Fractional path asymmetry \(\sim D/R_{\mathrm{CMB}}\).

Order-of-magnitude geometric temperature dipole (local estimate notes):

\[
\left(\frac{\delta T}{T}\right)_{\mathrm{geom}} \approx \alpha \,\frac{D}{R_{\mathrm{CMB}}}\,\left(\frac{\Delta H}{H}\right)_{\mathrm{eff}},\qquad \alpha\sim 0.5\text{–}2.
\]

**Migkas-sized contrast.** Take \((\Delta H/H)_{\mathrm{eff}}\sim 0.06\) (toy Amp\(\times\sin\) scale) to \(\sim 0.09\) (Migkas 2021 \(\sim 9\%\) \(H_0\)). For \(\alpha\sim 1\) and \((\delta T/T)_{\mathrm{geom}}\lesssim 1.2\times 10^{-3}\):

\[
\frac{D}{R_{\mathrm{CMB}}} \lesssim 0.02.
\]

If \(R_{\mathrm{CMB}}\sim 14\,\mathrm{Gpc}\), then \(D\lesssim\sim 0.3\,\mathrm{Gpc}\). Larger than classic Alnes \(\sim 15\,\mathrm{Mpc}\) because contrast is only \(\sim 6\)–\(9\%\), not order unity. Still incompatible with \(D/R\sim\mathcal{O}(0.1\)–\(0.5)\) without cancellation.

Details: `/workspace/cosmology/papers/gr_ch2_ltb_2026-09-20/CMB_DIPOLE_ESTIMATE.md`.

**Figures.**

- `/workspace/cosmology/figures/coe_geometry_cmb_dipole_schematic.png`
- `/workspace/cosmology/figures/photon_hits_per_arc_vs_angle.png`

---

## 5. Alnes–Kodama scope (preview; Module 7 deepens)

For **Gpc-scale Lemaître–Tolman–Bondi (LTB)** voids with **order-unity** contrast, **Alnes & Amarzguioui (2006)** force an off-center observer within \(\sim 15\,\mathrm{Mpc}\) of the center or the CMB dipole overshoots. **Kodama, Saito & Ishibashi (2010):** same class of bound, sometimes tens of Mpc — still \(\ll\) Gpc.

CMB constrains **spherical off-center + large contrast** hard. It does not automatically settle every CoE story at Migkas-sized contrast with \(R_{\mathrm{CMB}}\gg D\) — but the budget above still kills large \(D/R\).

**Kinematic Sunyaev–Zeldovich (kSZ)** pressure (García-Bellido & Haugbølle 2008) is flagged here and taught in Module 7: deep voids hurt even for centered observers.

---

## 6. Readings

- H. Alnes & M. Amarzguioui, *Phys. Rev. D* **74**, 103520 (2006); arXiv:astro-ph/0607334.
- H. Kodama, K. Saito & A. Ishibashi, arXiv:1004.3089 (2010).
- J. García-Bellido & T. Haugbølle, *JCAP* **09**, 016 (2008) — kSZ (flag now).
- E. Bertschinger, *ApJ* **268**, 17 (1983) — Sedov → cosmological self-similar blast (background).

---

## 7. Exercises

1. Why does a **uniform-density** ball fail as a pedagogical CoE for sky anisotropy? (Gauss / shell theorem.)

2. **Project thread: Dipole budget arithmetic.** For \(\Delta H/H=0.06\), \(\alpha=1\), max geom dipole \(1.2\times 10^{-3}\), compute max \(D/R_{\mathrm{CMB}}\). If someone proposes \(D/R_{\mathrm{CMB}}=0.3\), by what factor do they overshoot? How does allowed \(D\) compare to classic Alnes \(\sim 15\,\mathrm{Mpc}\) when contrast is only Migkas-sized?

3. Separate: “CMB is isotropic about us” vs “there is no center anywhere.” Which is closer to data?

---

## 8. Takeaway

The CMB dipole budget is a **gate** on \((D,(\Delta H/H)_{\mathrm{eff}})\), not a slogan. Module 7 places that gate next to LTB, Alnes–Kodama, and kSZ.
