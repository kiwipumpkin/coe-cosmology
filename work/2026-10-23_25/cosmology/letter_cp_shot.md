# Letter draft — A centered blast-wave Universe (short)

**Working title:** *Homology is isotropic; decreasing density is not: a first-principles angular match to Migkas et al. (2020)*

**Format note:** ApJL / A&A Letters / Nature Astronomy Correspondence length. Punchy, not a full article. Ontology explicit. No AI / Navier–Stokes meta commentary.

---

## Abstract (letter-style)

The cosmological principle (CP) asserts that there is no preferred center and that a suitably defined Hubble expansion is isotropic about every comoving observer. We take that claim seriously as an empirical target, not a theorem. Starting from a Sedov–Taylor (ST) blast-wave cosmos with a real geometric center of expansion (CoE) beyond the horizon, we show that pure homology \(\vec{v}=H\vec{r}\) yields \(H_{\mathrm{obs}}=H\) for *any* observer offset. Anisotropy requires homology breaking. With Gauss gravity \(g=-GM(<r)/r^2\) and a radially *decreasing* density \(\rho(r)\) (Phase A eq. 3.13), an off-center observer at distance \(D\) acquires a directional \(H_{\mathrm{obs}}(\hat{n})\). Mapping \(A\propto H^{-2}\) and comparing to digitised Migkas et al. (2020) Fig. 6 \(A/A_{\mathrm{all}}\) versus Galactic longitude—with only a phase (axis) rotation, not a fitted \(A\cos\varphi\) dipole template—gives \(R^2\simeq 0.54\) for the \(A\) series after best alignment (\(R^2\simeq 0.40\) before). That is not yet a kill shot on the CP, but it is a first-principles angular match of the kind the CP must survive. Open problems (CMB, light travel, Migkas cone-scan pipeline) are stated without evasion.

---

Dear Editor,

The cosmological principle is the quietest and most consequential assumption in modern cosmology: no preferred center; large-scale isotropy about us as typical observers. Late-time directional claims—especially Migkas et al. (2020, A&A 636, A15) on X-ray cluster \(L_X\)–\(T\) normalizations—have put pressure on that assumption at the several-to-tens-of-percent level. Those claims are contested (bulk flow, selection, scatter, simulation nulls). Contested or not, they define a concrete angular pattern against which a centered alternative can be tested *without* smuggling in a dipole template.

**Ontology.** We adopt a blunt working picture. The Universe has a real CoE. That center lies beyond the observable horizon. We sit off-center at Euclidean (Newtonian-pass) distance \(D\). Expansion is modeled as an ST blast wave about the CoE, with interior homology \(\vec{v}=H(t)\,\vec{r}\) plus central gravity from a decreasing density \(\rho(r)\). This is not “ΛCDM plus a small anisotropic patch.” The scientific goal is either to locate the CoE (direction and scale) or to kill the idea.

**Homology theorem (why geometry alone fails).** Place the CoE at the origin, the observer at \(\vec{D}\), a source at \(\vec{r}\), and \(\vec{d}=\vec{r}-\vec{D}\). Homologous ST gives \(\vec{v}_{\mathrm{rel}}=H\vec{d}\), hence
\[
H_{\mathrm{obs}}\equiv v_{\mathrm{los}}/d = H
\]
independent of direction. Pure offset does **not** anisotropize Hubble’s law. Drafts that claim otherwise usually omit subtracting the observer’s motion—a frame error that manufactures a fake dipole. Any real sky anisotropy must break homology.

**Gauss + decreasing \(\rho\) restores anisotropy.** Under spherical symmetry,
\[
\vec{g}=-\frac{GM(<r)}{r^2}\,\hat{r},\qquad M(<r)=\int_0^r 4\pi s^2\rho(s)\,ds.
\]
Uniform density is the rejected special case: \(M\propto r^3\Rightarrow g\propto -r\Rightarrow\) monopole-only \(H_{\mathrm{obs}}\). A decreasing profile—baseline Phase A eq. (3.13), the incomplete-gamma cumulative of \(\rho\propto e^{-r/r_s}\)—makes \(M(<r)\) grow slower than \(r^3\). Near observer scale \(r\sim D\), the field retains tidal \(1/r^2\) structure. Peculiar velocities from this field, differenced between source and observer and projected along \(\hat{n}\), yield
\[
H_{\mathrm{obs}}(\hat{n})=H+\delta v/d.
\]
At \(d\ll D\) the leading sky pattern is a quadrupole \(P_2(\mu)\) with \(\mu=\hat{n}\cdot\hat{D}\). The natural dimensionless amplitude is
\[
\kappa_{\mathrm{natural}}=\frac{GM_{\mathrm{enc}}(<D)}{D^3 H^2}=\frac{\Omega_m}{2}=0.15
\]
for \(\Omega_m=0.3\). We use \(\kappa=0.15\), \(\alpha\sim\mathcal{O}(1)\), and \(d/D=0.1\) as the binding baseline—not a free MCMC dipole amplitude.

**Migkas angular comparison (not a dipole fit).** Migkas Fig. 6 reports \(A/A_{\mathrm{all}}\) and \(B/B_{\mathrm{all}}\) versus Galactic longitude in \(\Delta\ell=90^\circ\) bins. We digitise those series and map the ST+Gauss \(H_{\mathrm{obs}}\) to \(A\) via \(A\times H^2\) conservation, i.e. \(A/A_{\mathrm{all}}\propto 1/(H_{\mathrm{obs}}/H)^2\), renormalised to mean unity. The comparison is **first-principles angular**: model \(H_{\mathrm{obs}}(\hat{n})\) through the same longitude binning, with a single phase \(\phi_0\) that aligns the preferred axis—not a fitted \(A\cos\varphi\) template.

Quantitatively (affine \(y=ax+b\) after amplitude scaling; see `migkas_r2_fit_summary.json`): for series **A**, \(R^2=0.402\) at \(\phi_0=0^\circ\), \(R^2=0.502\) at a Migkas-motivated prior \(\phi_0=280^\circ\), and \(R^2=0.536\) at best phase \(\phi_0\simeq 103.5^\circ\). Rotations improve the fit; the scatter plot of Migkas versus model is the primary diagnostic Chris prefers. Series **B** (slope) remains near-flat in both data and expectation (\(R^2\lesssim 0.06\))—a null that any selection or ICM add-on must not break.

**What this does to the CP.** The CP predicts \(A/A_{\mathrm{all}}\approx 1\) after local kinematic corrections. A decreasing-\(\rho\) centered blast wave predicts a coherent angular pattern at \(\kappa_{\mathrm{natural}}\). Matching Migkas’s *longitude* curve at \(R^2\sim 0.5\) with physics fixed by Gauss and \(\Omega_m/2\)—rather than by fitting a dipole template—is precisely the kind of late-time smoking gun the CP cannot shrug off if the morphology survives a fair observational filter. We are not claiming the CP is dead. We are claiming it is under first-principles fire.

**Honest open problems.**

1. **Morphology.** Full-sky Legendre content at \(d\ll D\) is quadrupole-dominated (\(C_2\sim 2\alpha\kappa\), \(C_1\approx 0\)). Migkas’s 1D \(A(\ell)\) looks more monopole+dipole-like (one deep trough near \(\ell\sim 250^\circ\)–\(300^\circ\)). Whether the Migkas \(\theta\simeq 75^\circ\) cone-scan pipeline converts a \(P_2\) sky into that 1D shape is an open forward-modelling question—not license to insert an \(A\cos\varphi\) by hand.

2. **CMB.** Off-center observers in centered cosmologies generically imprint CMB multipoles (Alnes & Amarzguioui; Saadeh et al. 2016 Bianchi limits). A CoE beyond the horizon without a last-scattering compatibility mechanism is the main adversary of this program. Late-time angular success does not waive that test.

3. **Light travel.** The present pass is instantaneous Newtonian \(v_{\mathrm{los}}/d\). Redshift composition, luminosity distance, and evolving positions remain to be done.

4. **Pipeline.** Longitude strips are not Migkas’s 2D cone maps. The next required step is an identical cone-scan forward map of \(H_{\mathrm{obs}}(\hat{n})\) with axis at \((280^\circ,-20^\circ)\) or joint \((303^\circ,-27^\circ)\), then retune \(\alpha,\kappa\) only against that filter.

**Bottom line.** Homology is isotropic. Decreasing density under Gauss’s law is not. With \(\kappa_{\mathrm{natural}}=0.15\), a phase-aligned first-principles \(A(\ell)\) already tracks Migkas Fig. 6 at \(R^2\simeq 0.54\). That is enough to put the cosmological principle on notice; it is not enough to retire it. The cone-scan, light-travel, and CMB tests decide whether this letter becomes a footnote or a funeral oration.

Respectfully submitted,

[Authors]

---

## Recommended figures (letter: 3–4 max)

**Fig. 1 — Geometry.**  
File: `/workspace/cosmology/figures/geometry.png`  
Caption: Off-center observer in a centered ST cosmos. CoE at the origin; observer at \((D,0)\); sample sources on a circle of radius \(R\sim d\). Preferred axis \(\hat{D}\), CoE-radius \(r\), Earth-centric angle \(\varphi\). Pure homology \(\vec{v}=H\vec{r}\) → isotropic \(H_{\mathrm{obs}}\); anisotropy requires Gauss gravity with decreasing \(\rho\).

**Fig. 2 — Migkas vs model (money plot; letter-ready).**  
Preferred: `/workspace/cosmology/figures/migkas_r2_letter_panels.png` (compact 2×2).  
Full diagnostic: `migkas_r2_fit_panels.png`. Overlay companion: `migkas_fig6_model_overlay.png`.  
Caption: Digitised Migkas et al. (2020) Fig. 6 \(A/A_{\mathrm{all}}\) vs Galactic longitude with ST+Gauss decreasing-\(\rho\) prediction (\(\kappa=0.15\), \(\alpha=1\), \(d/D=0.1\)), mapped \(A\propto H^{-2}\). Scatter Migkas vs model + plain \(R^2\). No \(A\cos\varphi\) dipole template. Best \(A\)-series \(R^2=0.536\) at \(\phi_0\simeq 103.5^\circ\) (\(R^2(0)=0.402\); prior \(\phi_0=280^\circ\): \(R^2=0.502\)). \(B\) null \(R^2\simeq 0.058\).

**Fig. 3 — Phase alignment.**  
File: `/workspace/cosmology/figures/migkas_r2_vs_rotation.png`  
Caption: \(R^2\) vs model phase \(\phi_0\). Rotations improve the \(A\) fit; \(\phi_0=280^\circ\) Migkas prior marked; \(B\) stays a weak null. Axis alignment of a first-principles sky — not a free dipole amplitude fit.

**Fig. 4 (optional).**  
`density_mass.png` or `sky_mollweide.png` — omit from short letter if space-tight (save for longer paper).

**Numbers:** `/workspace/cosmology/figures/migkas_r2_fit_summary.json`.
