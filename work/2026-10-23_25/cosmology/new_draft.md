# Anisotropic Hubble Expansion from a Central-Mass Sedov–Taylor Model of the Universe (New Draft)

**Author working hypothesis:** Chris — Cosmology agent program  
**Draft date:** 2026-09-17 (America/Chicago); **language revised** same day (spherical decreasing-\(\rho\) LIVE clarification).  
**Status:** Complete new prose draft superseding the unfinished 2020s outline/extract. Incorporates cleaned Phase A equations, decreasing-\(\rho\) baseline (3.13), Gauss gravity, time-integrated peculiar velocities, \(\kappa_{\mathrm{natural}}\), and multipole honesty. Old draft.txt left untouched. **Binding 2026-09-17:** decreasing \(\rho(r)\) (eq. 3.13) + Gauss \(g=-Gm(r)/r^2\); uniform ball rejected. **Spherical decreasing-\(\rho\) is LIVE** (amplitude GO; Migkas-like angular curve exists); remaining caveat is multipole morphology vs a *pure dipole*, not a veto of spherical variable-density cosmology.  
**Tone:** Skeptical graduate thesis / paper draft. Not a discovery claim.

---

## Abstract

The cosmological principle asserts that the Universe has no preferred center and that a suitably defined Hubble expansion is isotropic about every comoving observer. Late-time probes—especially directional analyses of X-ray galaxy-cluster scaling relations (Migkas et al. 2020, 2021)—have reported apparent anisotropies at the several-to-tens-of-percent level in effective Hubble parameters or luminosity distances. Those claims remain contested: bulk flows, selection, scatter, and simulation-based null tests can all mimic directional signals. This draft develops, and then stress-tests, a concrete alternative ontology: the Universe expands as a Sedov–Taylor (ST) blast wave about a **real geometric center of expansion (CoE)** that lies beyond the observable horizon, and we are an **off-center observer**.

Starting from classical ST similarity and a homologous interior velocity field \(\vec{v}=H(t)\,\vec{r}\), we prove a elementary but decisive kinematic result: pure homology yields an **isotropic** apparent Hubble parameter \(H_{\mathrm{obs}}=H\) for *any* observer offset. Anisotropy therefore cannot come from geometry alone; it requires a homology-breaking mechanism. The working mechanism studied here follows from **Gauss’s law / Newton’s shell theorem** together with a **centrally concentrated, radially decreasing density** \(\rho(r)\) (baseline ontology—not a uniform ball). Under spherical symmetry,
\[
\vec{g}=-\frac{G M(<r)}{r^2}\,\hat{r},
\]
with \(M(<r)=\int_0^r 4\pi s^2\rho(s)\,ds\). Chris identifies the relevant enclosure with the mass toward the CoE; at observer scale this sets \(\kappa_{\mathrm{natural}}=G M_{\mathrm{enc}}(<D)/(D^3 H^2)=\Omega_m/2=0.15\) (\(\Omega_m=0.3\)). The default profile is the Phase A extended law (eq. 3.13)—incomplete-gamma cumulative of \(\rho\propto e^{-r/r_s}\)—or any similar ST-motivated decreasing profile. Peculiar velocities are **time-integrated** along ST trajectories. **Why decreasing density matters:** uniform \(\rho\Rightarrow M(<r)\propto r^3\Rightarrow g\propto -r\Rightarrow\) monopole-only \(H_{\mathrm{obs}}\) (no anisotropy). Decreasing \(\rho\) makes \(M(<r)\) grow slower than \(r^3\), so near \(r\sim D\) the field retains tidal \(1/r^2\) structure and **restores** sky anisotropy (leading quadrupole at \(d\ll D\)).

For this baseline Gauss + decreasing-\(\rho\) model beyond the horizon, the predicted sky is **not** a clean Migkas dipole: leading anisotropy is a **quadrupole** \(C_2\approx 2\alpha\kappa\), full-sky Legendre \(C_1\approx 0\) at \(\mathcal{O}(d/D)\). Fitted Migkas-style \(A\sim 0.1\) can appear as **selection / footprint leakage** without a true orthogonal dipole. Relative acceleration versus angle—Chris’s intended Migkas-curve mechanism—is presented honestly, then confronted with this multipole structure.

Parameters \((D,M_{\mathrm{enc}},\xi,\alpha_{\mathrm{eff}})\) remain largely unconstrained. The CMB remains the main adversary of a real beyond-horizon center. This draft does **not** claim to have smoked the cosmological principle. It states the working hypothesis clearly, corrects errors in the unfinished draft, and defines which analyses are worth doing next and which are not.

---

## 1. Introduction

### 1.1 Hubble tension, isotropy, and the cosmological principle

The standard ΛCDM model assumes statistical homogeneity and isotropy on large scales—the cosmological principle (CP). In that framework, Hubble’s law \(v=Hd\) is a local consequence of FLRW expansion, not evidence for a geometric center. Dark energy (or an equivalent modification) is invoked to explain late-time acceleration. Parallel to the well-known \(H_0\) tension between early- and late-Universe determinations, a second class of tension concerns **directional** variations: reports that effective \(H_0\), luminosity distance, or cluster scaling normalizations depend on sky direction at levels far above the kinematic dipole expected from our motion relative to the CMB.

If those late-time anisotropies are cosmological rather than systematic, they challenge the CP more sharply than a mere recalibration of \(H_0\). The CP is not a theorem; it is an empirical claim that has been extraordinarily successful. Smoking it requires more than a preferred axis in one catalog. It requires a coherent model that (i) predicts the observed late-time pattern, (ii) survives bulk-flow degeneracies, and (iii) does not destroy the near-isotropy of the CMB.

### 1.2 Motivation for a real center

The working hypothesis of this research program is stronger than “ΛCDM plus a small anisotropic perturbation.” It asserts:

1. The Universe has a **real center of expansion (CoE)**—a preferred origin for a blast-wave-like expansion.
2. That CoE lies **beyond the observable horizon**, so we do not see it as a local object.
3. We are an **off-center observer** at comoving (or Euclidean, in the Newtonian pass) distance \(D\) from the CoE.
4. Expansion is usefully modeled, at least at the Newtonian skeleton level, by a **Sedov–Taylor blast wave** with a **radially decreasing** density \(\rho(r)\) about the CoE; gravity follows Gauss’s law \(g=-G M(<r)/r^2\) with \(M_{\mathrm{enc}}\) the mass toward the center (Chris). Uniform density is the rejected special case.

The scientific goal is to **locate that center** (direction and distance scale) and constrain enclosed mass / similarity structure—or to **kill the idea** if CMB isotropy and related probes make a real center untenable. This is not a side project on dipoles; it is an attempt to replace Copernican cosmology with a centered ontology if the data force that conclusion.

### 1.3 Scope of this draft

This document supersedes the unfinished draft (local `draft.txt`). It retains the title and the ST + central-mass motivation, but replaces heuristic claims with cleaned equations, honest multipole predictions, and explicit GO / CONDITIONAL / NO-GO gates (spherical decreasing-\(\rho\) LIVE; dipole-only Migkas MCMC CONDITIONAL). Sections 3–5 incorporate results from the living notes: Phase A derivation, absurdity check (with enclosed-mass addendum), multipole structure, and dipole-fraction / footprint analyses. Relativistic and CMB issues are stated as open problems, not solved.

What this draft is **not**: a peer-reviewed localization of the CoE; a GR metric; a CMB-compatible construction; or a claim that Migkas et al. have already falsified the CP.

---

## 2. Observational landscape

### 2.1 Migkas et al. (2020, 2021)

Migkas et al. (2020, A&A 636, A15) analyzed the X-ray luminosity–temperature (\(L_X\)–\(T\)) relation in a new galaxy-cluster sample and reported directional variation in the scaling normalization. For 313 MCXC clusters they found a \(\gtrsim 4\sigma\) anisotropy toward Galactic \((l,b)\sim(280^\circ,-20^\circ)\); joint analyses with other samples reached \(\sim 5\sigma\) toward \(\sim(303^\circ,-27^\circ)\). One interpretation is a \(\sim 16\pm 3\%\) luminosity-distance variation across the sky; alternatives include X-ray absorption, calibration, and selection systematics.

Migkas et al. (2021, A&A 649, A151) extended the analysis to ten cluster scaling relations and reported an apparent \(\sim 9\%\) \(H_0\) dipole toward \(\sim(280^\circ,-15^\circ)\), significant at \(>5\sigma\) in their Monte Carlo tests—**or**, equivalently, a bulk flow of order \(\sim 900\,\mathrm{km\,s^{-1}}\) out to \(\gtrsim 500\,\mathrm{Mpc}\). That degeneracy is central: a coherent peculiar velocity field can mimic anisotropic expansion until higher-redshift clusters (e.g. eROSITA) break the redshift scaling.

For this program, Migkas provides the best all-sky late-time anisotropy map against which to **aim** a CoE direction. It is not proof of a geometric center.

### 2.2 Related late-time dipoles and bulk flows

Related claims include cluster \(\sigma_v\) scaling dipoles roughly consistent in direction with Migkas but amplitude-sensitive to scatter (Pandya et al. 2024), and weaker \(\sim 3\%\) Tully–Fisher \(H_0\) variations (Boubel et al. 2025). Classical bulk-flow measurements (Watkins, Feldman & Hudson 2009) find \(\sim 400\,\mathrm{km\,s^{-1}}\) on \(\sim 50\,h^{-1}\,\mathrm{Mpc}\) scales—large relative to vanilla ΛCDM expectations in some analyses, but still well below the \(\sim 900\,\mathrm{km\,s^{-1}}\) Migkas-equivalent scale. Contested “dark flow” claims from kinetic Sunyaev–Zel’dovich stacking (Kashlinsky et al.) illustrate how easily preferred directions appear and how hard they are to establish.

SN Ia residual dipoles after CMB-frame correction are mixed and often consistent with local structure. Quasar and radio-source dipoles that exceed the kinematic expectation remain contested; if real, any CoE model must either align with those axes or explain misalignment.

### 2.3 Critiques and simulation nulls

FLAMINGO-based and related mock analyses of scaling relations as isotropy tracers indicate that isotropic ΛCDM lightcones plus realistic scatter can produce large *apparent* anisotropies. Reported significance of Migkas-like signals can drop to the \(\sim 3\sigma\) class, with much of the “signal” attributable to noise and scatter. In this draft, Migkas amplitudes are treated as an **upper target** for model building, not as gospel.

Methodologically one must always separate: (a) detection of directional variation in fitted normalizations; (b) cosmological anisotropic expansion; (c) bulk flow; (d) X-ray absorption, calibration, and selection.

### 2.4 CMB isotropy — the main adversary

Planck temperature and polarization are consistent with statistical isotropy at high accuracy (Planck Collaboration 2018 VII / A&A 641, A7). Known large-angle anomalies (hemispherical asymmetry, low-\(\ell\) alignments) do **not** license an \(\mathcal{O}(10\%)\) late-time Hubble dipole from a geometric center. Global anisotropic expansion is strongly disfavored by Bianchi analyses (Saadeh et al. 2016, Phys. Rev. Lett. 117, 131302): odds \(\sim 10^5:1\) against, with \((\sigma_V/H)_0 < 4.7\times 10^{-11}\) (95%).

Off-center observers in spherically symmetric inhomogeneous cosmologies (LTB-like) are tightly constrained by the CMB dipole: Alnes & Amarzguioui (2006) find offsets \(\lesssim\sim 15\,\mathrm{Mpc}\) for their void models; Kodama et al. (2010) update analytic dipole/quadrupole formulae. A CoE *beyond the horizon* is a stronger claim than a Gpc void with us near center. It requires a mechanism by which last scattering looks nearly isotropic while \(z\lesssim 0.3\) clusters show \(\mathcal{O}(10\%)\) effects. Without that mechanism, the CMB wins. Section 7 returns to this as an open problem, not a solved loophole.

---

## 3. Mathematical framework

### 3.1 Ontology and assumptions (Newtonian pass)

**Working ontology.** CoE at the coordinate origin. Observer at \(\vec{D}\) with \(D=\|\vec{D}\|\). Sources at \(\vec{r}\) with CoE-radius \(r=\|\vec{r}\|\) and observer–source separation \(\vec{d}=\vec{r}-\vec{D}\), \(d=\|\vec{d}\|\). Line of sight \(\hat{n}=\vec{d}/d\). Preferred axis \(\hat{D}=\vec{D}/D\) points from CoE to observer; the direction from observer toward the CoE is \(-\hat{D}\).

**Demo angle convention.** Toy calculations place sources on an observer-centered circle of radius \(R\) with Earth-centric angle \(\phi\) measured from \(+\hat{D}\) (anti-center): \(\hat{n}=(\cos\phi,\sin\phi)\). Thus \(\phi=0\) looks away from the CoE and \(\phi=\pi\) looks toward it. When prose says “toward the center,” it means \(\theta=\pi-\phi\) unless stated otherwise. The unfinished draft mixed these conventions; this draft does not.

**Phase A assumptions.**

1. Newtonian kinematics in Euclidean 3-space (planar demos are equatorial slices).
2. Instantaneous velocities at a fixed cosmic time \(t\) (no light travel yet).
3. Dust-like medium; classical ST similarity for the blast-front background.
4. Interior velocity **homologous** about the CoE: \(\vec{v}=H(t)\,\vec{r}\).
5. Observer and sources are comoving tracers, plus optional peculiar \(\vec{u}\) from gravity.
6. **Baseline mass: decreasing \(\rho(r)\).** Spherical symmetry \(\Rightarrow\) Gauss \(g=-G M(<r) r^{-2}\hat{r}\). Default profile: Phase A eq. (3.13) (\(\rho\propto e^{-r/r_s}\)) or similar ST-motivated decreasing law; \(M_{\mathrm{enc}}(<D)\) sets \(\kappa_{\mathrm{natural}}\). **Rejected special case:** uniform ball \(\rho=\mathrm{const}\Rightarrow M\propto r^3\Rightarrow g\propto -r\Rightarrow\) monopole only.
7. “Observed” distance for \(H_{\mathrm{obs}}\) is Euclidean \(d\), not luminosity distance.

Out of scope in this section: GR metrics, redshift composition, CMB last scattering, beyond-horizon matching.

### 3.2 Sedov–Taylor background (eqs 3.1–3.3)

Classical Sedov–Taylor blast front in uniform ambient density \(\rho_0\) with energy \(E\):

\[
R_{\mathrm{ST}}(t)=\xi\left(\frac{E}{\rho_0}\right)^{1/5}t^{2/5}. \tag{3.1}
\]

Front speed:

\[
v_{\mathrm{front}}(t)=\frac{dR_{\mathrm{ST}}}{dt}=\frac{2}{5}\frac{R_{\mathrm{ST}}}{t}. \tag{3.2}
\]

Hubble-like rate from the front:

\[
H(t)\equiv\frac{\dot{R}_{\mathrm{ST}}}{R_{\mathrm{ST}}}=\frac{2}{5t}\qquad\Rightarrow\qquad v_{\mathrm{front}}=H(t)\,R_{\mathrm{ST}}. \tag{3.3}
\]

Interior working assumption:

\[
\vec{v}(\vec{r},t)=H(t)\,\vec{r}.
\]

**Caveat.** Cosmological identification of \((E,\rho_0,\xi)\) is undefined. Classical ST assumes a *static* uniform ambient medium; that is cosmologically false. Equations (3.1)–(3.3) are a kinematic analogy and similarity scaffold, not a derivation of FLRW.

Notation discipline: \(\xi\) is the ST similarity coefficient in (3.1); \(H=\dot{R}_{\mathrm{ST}}/R_{\mathrm{ST}}\); sheet parameter \(k\) is reserved for literal spreadsheet matching and should not be confused with \(\xi\). Blast-front radius is \(R_{\mathrm{ST}}\); toy circle radius is \(R\).

### 3.3 Geometry (eqs 3.4–3.5)

\[
r=\|\vec{r}\|=\sqrt{P_x^2+P_y^2+P_z^2}, \tag{3.4}
\]

\[
d=\|\vec{r}-\vec{D}\|=\sqrt{(P_x-D)^2+P_y^2+P_z^2}. \tag{3.5}
\]

Never use the same letter for CoE-radius and observer–source distance.

### 3.4 Homology isotropy theorem (eqs 3.6–3.8) — critical result

Observer velocity \(\vec{v}_O=H\vec{D}\); source velocity \(\vec{v}_S=H\vec{r}\). Relative velocity:

\[
\vec{v}_{\mathrm{rel}}=\vec{v}_S-\vec{v}_O=H(\vec{r}-\vec{D})=H\,\vec{d}. \tag{3.6}
\]

Line-of-sight recession:

\[
v_{\mathrm{los}}=\vec{v}_{\mathrm{rel}}\cdot\hat{n}=H\,d. \tag{3.7}
\]

Apparent Hubble parameter:

\[
H_{\mathrm{obs}}=\frac{v_{\mathrm{los}}}{d}=H(t). \tag{3.8}
\]

**Theorem (homology ⇒ isotropy).** For any observer offset \(D\), a purely homologous field \(\vec{v}=H\vec{r}\) yields \(H_{\mathrm{obs}}=H\) **independent of direction**. Geometry alone does not produce Migkas-like anisotropy. This is the classical reason that \(v=Hd\) does not imply a center.

The unfinished draft claimed after (3.8) that dependence on offset \(D\) “therefore introduces anisotropy.” That claim is **false** once observer motion is subtracted. Anisotropy requires breaking homology: gravity, non-homologous ST interior structure, selection on non-homologous shells, or relativistic / light-cone effects.

**Cautionary failure mode.** Using absolute CoE-frame speed \(H\vec{r}\cdot\hat{n}\) without subtracting \(\vec{v}_O\) falsely produces

\[
H_{\mathrm{obs}}^{\mathrm{(wrong)}}=H\left(1+\frac{D}{R}\cos\phi\right)
\]

on an observer-centered circle. That is a frame error, not physics. Any numerical sheet or demo that shows angular variation at \(M=0\) is broken.

### 3.5 Gravity from Gauss’s law with decreasing \(\rho(r)\) (eqs 3.9–3.12)

**Binding ontology.** Density is **not** uniform. The baseline is a **centrally concentrated** spherical profile with \(\rho(r)\) **decreasing outward** about the CoE (ST-motivated; Phase A extended mass). Uniform \(\rho\) is a rejected special case (§4.1).

**Gauss / shell theorem.** Under spherical symmetry the field at radius \(r\) depends only on enclosed mass:

\[
\vec{g}(\vec{r})=-\frac{G M(<r)}{r^2}\,\hat{r},\qquad
M(<r)=\int_0^r 4\pi s^2\rho(s)\,ds. \tag{3.9}
\]

Exterior shells contribute zero. This is not an “ad hoc point mass versus enclosed mass” dichotomy: \(1/r^2\) *is* the enclosed-mass law. What changes the phenomenology is the **radial run of \(\rho\)**, which determines how \(M(<r)\) scales with \(r\).

**Why decreasing density restores anisotropy.**

| Density | \(M(<r)\) | \(g(r)\) | Sky structure in \(H_{\mathrm{obs}}\) |
|---------|-----------|----------|----------------------------------------|
| Uniform \(\rho=\mathrm{const}\) | \(\propto r^3\) | \(\propto -r\) | **Monopole only** — anisotropy vanishes |
| Decreasing \(\rho(r)\) (baseline) | grows slower than \(r^3\); \(\to M_{\mathrm{enc}}\) at large \(r\) | \(\sim -GM_{\mathrm{enc}}/r^2\) near/outside concentration | **Tidal quadrupole** (leading) at \(d\ll D\); \(C_1\approx 0\) |

Uniform mean density fills every radius equally, so differential peculiar velocity between observer and source is purely proportional to separation \(\vec{d}\) and cancels into a direction-independent renormalization of \(H\). A decreasing profile concentrates mass toward the CoE: at observer distances \(r\sim D\), \(M(<r)\) is already close to the total enclosure \(M_{\mathrm{enc}}(<D)\), the force retains \(1/r^2\) tidal structure, and **sky anisotropy returns**. That is the physical reason the draft’s baseline is decreasing \(\rho\), not a uniform ball.

**Working enclosed mass and \(\kappa\).** Chris identifies \(M_{\mathrm{enc}}\) with the mass of everything between us and the CoE. Then

\[
\kappa_{\mathrm{natural}}=\frac{G M_{\mathrm{enc}}(<D)}{D^3 H^2}=\frac{\Omega_m}{2}=0.15\quad(\Omega_m=0.3)
\]

normalizes the Gauss field at observer scale. The default cumulative profile is Phase A (3.13) (§3.6); asymptotic exterior behavior matches \(g=-GM_{\mathrm{enc}}/r^2\).

**Peculiar velocity.** Acceleration is not velocity. Bolt-on placeholder:

\[
\vec{u}(\vec{r})=-\alpha\frac{G m(r)}{r^2 H}\,\hat{r}, \tag{3.10}
\]

with \(m(r)=M(<r)\) from the decreasing profile and \(\alpha\sim\mathcal{O}(1)\). Primary physics replaces (3.10) by **time integration** along ST trajectories (Sec. 3.7) under the same Gauss law.

Relative recession:

\[
\vec{v}_{\mathrm{rel}}=H\vec{d}+\bigl[\vec{u}(\vec{r})-\vec{u}(\vec{D})\bigr],
\]

\[
v_{\mathrm{los}}=H d+\delta v,\qquad \delta v=\bigl[\vec{u}(\vec{r})-\vec{u}(\vec{D})\bigr]\cdot\hat{n}. \tag{3.11}
\]

\[
H_{\mathrm{obs}}(\hat{n},d)=H+\frac{\delta v}{d}. \tag{3.12}
\]

**Why the unfinished draft’s (3.11)–(3.12) are rejected.** The draft wrote roughly \(v_{\mathrm{los}}\approx H_0 d-(GM/d^2)\cos\theta\) and \(H_{\mathrm{obs}}=H_0-(GM/d^3)\cos\theta\). Problems: (i) missing \(1/H\) (acceleration ≠ velocity); (ii) gravity depends on CoE-radius \(r\), not observer distance \(d\); (iii) omitted \(\vec{u}(\vec{D})\); (iv) unspecified angle projection; (v) hydro double-counting risk. The **Gauss + decreasing-\(\rho\)** force law is retained; the draft’s kinematic projection is what fails.

### 3.6 Baseline decreasing mass profile (eqs 3.13–3.14)

**Default ontology for \(\rho(r)\).** The working cumulative enclosed-mass profile (toy sheet / Phase A) is

\[
\frac{m(r)}{M}=1-e^{-X}\left(1+X+\frac{X^2}{2}\right),\qquad X=\frac{r}{r_s}. \tag{3.13}
\]

This is the incomplete-gamma cumulative of an **exponentially decreasing** density \(\rho\propto e^{-r/r_s}\) (specifically \(m(r)\propto\int_0^r s^2 e^{-s/r_s}\,ds\)). Scale radius \(r_s\) (sheet \(L\)) sets how rapidly mass concentrates toward the CoE. Total asymptotic mass \(M\) is identified with \(M_{\mathrm{enc}}\) for \(\kappa_{\mathrm{natural}}\) (sheet quirk \(M\approx 1.142\ldots\) ignored as a unit artifact).

\[
\vec{g}(r)=-\frac{G m(r)}{r^2}\,\hat{r}. \tag{3.14}
\]

Equation (3.14) is Gauss’s law with \(M(<r)=m(r)\). For \(r\gg r_s\), \(m(r)\to M_{\mathrm{enc}}\) and \(g\to -GM_{\mathrm{enc}}/r^2\)—the tidal limit used in multipole expansions. For \(r\lesssim r_s\), the soft core regularizes the origin. The profile is **phenomenological / ST-motivated**, not yet derived from a full ST+gravity similarity solution; any similarly decreasing law (e.g. other ST interior scalings) that yields \(M(<r)\) growing slower than \(r^3\) preserves the same qualitative anisotropy restoration.

Residual map (optional):

\[
\Psi(\hat{n},d)=H_{\mathrm{obs}}(\hat{n},d)-H. \tag{3.15}
\]

Pure ST: \(\Psi=0\). With gravity: \(\Psi=\delta v/d\).

### 3.7 Time-integrated gravity (primary, not bolt-on)

Chris’s binding correction: gravity should be **time-integrated from the blast**, not treated as an instantaneous peculiar bolt-on. The bolt-on is retained only as a shortcut / regression baseline.

Background ST homology: \(R_{\mathrm{ST}}\propto t^{2/5}\), \(H=2/(5t)\), \(r(t)=\chi\,t^{2/5}\). Extra central mass with \(g=-GM/r^2\) (still not a fully self-consistent ST+gravity similarity solution). Three treatments:

1. **Frozen path, no drag:** \(\Delta v=\int_{t_i}^{t_f} g(r_{\mathrm{bg}}(t))\,dt\).
2. **Hubble drag:** \(du/dt+Hu=g\) on the homologous background.
3. **Nonlinear ODE:** \(\ddot{r}=-(6/25)r/t^2-GM/r^2\) (ST radial acceleration preserving homology when \(M=0\), plus gravity).

Effective \(\alpha\) relative to bolt-on \(u=\alpha\,GM/(r^2 H)\):

| Treatment | \(\alpha_{\mathrm{eff}}\) (\(\eta=t_i/t_f\to 0\)) | \(\alpha_{\mathrm{eff}}\) (\(\eta=10^{-3}\)) |
|-----------|-----------------------------------------------|---------------------------------------------|
| Frozen \(\int g\,dt\) | 2 | \(\approx 1.50\) |
| Drag \(du/dt+Hu=g\) | \(2/3\) | \(\approx 0.66\) |
| Bolt-on Phase A | free \(\sim\mathcal{O}(1)\) | 1 (reference) |

**Consequence.** Integration rescales amplitudes by \(\mathcal{O}(1)\) only. It does **not** invent a new small-\(\kappa\) regime. Required \(\kappa\) for fixed sky variation scales as \(1/\alpha_{\mathrm{eff}}\). Null tests: \(M=0\) preserves homology and isotropic \(H_{\mathrm{obs}}\) in bolt-on and ODE implementations.

### 3.8 Enclosed mass and natural \(\kappa\)

**Binding reinterpretation.** \(M_{\mathrm{enc}}\) normalizing \(\kappa\) is the **enclosed cosmic mass** toward the CoE under the baseline **decreasing** \(\rho(r)\). By Gauss it sources \(g=-G m(r)/r^2\). It is not an extra microscopic point mass, and it is not uniform density:

\[
M_{\mathrm{enc}}(<D)=\frac{4\pi}{3}\,\Omega_m\,\rho_c\,D^3=\Omega_m\,M_{\mathrm{crit}}(<D).
\]

Dimensionless coupling:

\[
\kappa\equiv\frac{GM}{D^3 H^2}.
\]

Critical density \(\rho_c=3H^2/(8\pi G)\) implies \(G M_{\mathrm{crit}}(<D)/(D^3 H^2)=1/2\), hence

\[
\boxed{\kappa_{\mathrm{natural}}=\frac{G M_{\mathrm{enc}}(<D)}{D^3 H^2}=\frac{\Omega_m}{2}}.
\]

For \(\Omega_m=0.3\):

\[
\kappa_{\mathrm{natural}}=0.15.
\]

This \(\kappa_{\mathrm{natural}}\) normalizes the baseline decreasing-\(\rho\) Gauss model at \(r\sim D\). Required \(\kappa\sim 0.04\)–\(0.10\) for \(\sim 10\%\) half peak-to-peak variation (frozen / drag / ODE at \(d/D=0.1\)) is a fraction of \(\kappa_{\mathrm{natural}}\)—**fully natural**. “Absurd extra point mass” language is retired. Cluster comparisons are expected for cosmic \(M_{\mathrm{enc}}\) and are not a veto.

Time integration remains \(\mathcal{O}(1)\) and does not alter this conclusion.

---

## 4. Predicted sky structure

**Baseline throughout this section:** decreasing \(\rho(r)\) (eq. 3.13), Gauss \(g=-G m(r)/r^2\), \(\kappa=\kappa_{\mathrm{natural}}=0.15\), \(\alpha_{\mathrm{eff}}\sim\mathcal{O}(1)\). In the tidal / exterior regime \(m(r)\approx M_{\mathrm{enc}}\) so \(g\approx -GM_{\mathrm{enc}}/r^2\). Section 4.1 is the rejected uniform special case only.

### 4.1 Rejected special case: uniform \(\rho\Rightarrow g\propto -r\Rightarrow\) monopole only

If \(\rho=\mathrm{const}\), then \(M(<r)\propto r^3\) and Gauss gives \(g=-(4\pi G\rho/3)r\). Differential peculiar motion yields

\[
H_{\mathrm{obs}}=H(1-\alpha\kappa)
\]

independent of \(\hat{n}\)—**monopole only; anisotropy vanishes**. This is why uniform density is **not** the baseline. Decreasing \(\rho\) is required to restore angular structure (abstract / §3.5).

### 4.2 Baseline decreasing-\(\rho\) / Gauss model: leading quadrupole at \(d\ll D\)

With \(m(r)\approx M_{\mathrm{enc}}\) (or exact \(m(r)\) from 3.13) and \(\kappa_{\mathrm{natural}}\), the tidal limit \(\varepsilon=d/D\to 0\) yields

\[
\frac{H_{\mathrm{obs}}}{H}-1=\alpha\kappa\,(3\mu^2-1)=2\alpha\kappa\,P_2(\mu),
\]

where \(\mu=\hat{n}\cdot\hat{D}\). This is a **pure quadrupole**: \(C_2\to 2\alpha\kappa\), \(C_1\to 0\). Series in \(\varepsilon\):

\[
\mathrm{frac}=\alpha\kappa(3\mu^2-1)+\alpha\kappa\,\varepsilon\,(4.5\,\mu-7.5\,\mu^3)+\mathcal{O}(\varepsilon^2).
\]

The \(\mathcal{O}(\varepsilon)\) odd piece projects onto \(P_3\) (octupole): \(C_3\approx -3\alpha\kappa\,\varepsilon\). **Full-sky Legendre dipole \(C_1=0\) at this order.**

At \(\kappa_{\mathrm{natural}}=0.15\), \(\alpha=1\), \(d/D=0.1\): \(C_2\approx 0.30\), half peak-to-peak \(\approx 0.25\) (Migkas-*scale variation* natural under the baseline decreasing-\(\rho\) Gauss model), while dipole amplitude \(A\sim 10^{-7}\) (useless as a true Migkas dipole). Drag \(\alpha_{\mathrm{eff}}\approx 0.66\) scales amplitudes by \(\alpha_{\mathrm{eff}}\) without changing morphology. Finite \(r_s\) in (3.13) softens the core; multipole ordering at \(d\ll D\) remains quadrupole-led.

### 4.3 When a large true dipole appears

A substantial full-sky dipole under this baseline spherical model appears only when \(d/D\gtrsim 1\) (sources reach toward / past the CoE scale). That regime conflicts with the beyond-horizon ontology for Migkas cluster depths (hundreds of Mpc with \(D\) beyond the horizon \(\Rightarrow d/D\sim 0.01\)–\(0.1\)). Near-CoE singularities of a pure \(M=\mathrm{const}\) exterior law are softened by any finite core in \(m(r)\), but the beyond-horizon multipole conclusion is unchanged.

### 4.4 Fitted Migkas dipoles vs orthogonal multipoles

**Distinguish carefully:**

- **Full-sky Legendre \(C_1\):** orthogonal dipole coefficient. \(\approx 0\) at \(d\ll D\) for the baseline decreasing-\(\rho\) tidal field.
- **Fitted Migkas \(A\):** coefficient in \(A\cos\phi\) or \(M+A\mu\) on a sky mask. Can be large even when \(C_1\approx 0\).

On a CoE-facing hemisphere at \(d/D=0.1\), \(\kappa=0.15\), \(\alpha=1\): fake \(A_{\mathrm{moment}}\approx 0.11\) and \(A_{\mathrm{Migkas}}(M+A\mu)\approx 0.41\), while full-sky \(C_1\sim 10^{-8}\). Equatorial \(\phi\)-uniform Fourier leaks octupole into \(a_1\sim\mathcal{O}(d/D)\) (\(A_{\mathrm{leak}}\approx 0.375\,C_3\approx -0.017\) at \(d/D=0.1\))—nonzero fitted \(A\) without true \(C_1\), but still below typical Migkas \(A\sim 0.05\)–\(0.15\) targets on full sky.

Pure-dipole Migkas template on **full sky** at \(d/D\le 0.2\): \(R^2_{\mathrm{dip}}\lesssim 0.02\), while quadrupole template \(R^2_Q\gtrsim 0.92\). Same template at \(d/D=2\): \(R^2_{\mathrm{dip}}\approx 0.98\)—but then the CoE is not beyond survey depth.

**Softened language.** Do not say “a quadrupole has no dipole.” Say: full-sky Legendre \(C_1\approx 0\) at \(\mathcal{O}(d/D)\); fitted Migkas \(A\) can still be nonzero on partial sky or \(\phi\)-uniform measure.

### 4.5 Relative acceleration vs distance to CoE (Chris’s Migkas-curve mechanism)

Chris’s intended mechanism for the Migkas angular curve is that sources at different angles sit at different distances from the CoE and therefore experience different time-integrated Gauss accelerations from decreasing \(\rho\); the observer-centric \(H_{\mathrm{obs}}(\phi)=\Delta v_{\mathrm{los}}/d\) then varies with angle. That kinematics is exactly what Sections 3.5–3.7 and 4.2 expand under the primary force law.

Presented honestly: the mechanism **does** produce an angular curve of Migkas-relevant *variation amplitude* at natural \(\kappa\) on the baseline decreasing-\(\rho\) model. **Spherical decreasing-\(\rho\) is therefore viable** as a cosmology skeleton: amplitude GO, angular curve exists. Multipole honesty then requires stating that for beyond-horizon \(D\) and spherical \(M_{\mathrm{enc}}\), that curve is **quadrupole-led**, not a pure \(A\cos\phi\) dipole—so the fair test is a \(P_2\) comparison, not a dipole-only MCMC. (Uniform \(\rho\) would give no curve at all—§4.1.) Visual resemblance after Mollweide rotation, noise, and binning—as in the unfinished draft’s demos—is not a multipole identification.

---

## 5. Comparison with Migkas et al.

### 5.1 What matches visually

Unfinished-draft demos (recession curves vs angle; 3D surfaces; Mollweide maps rotated toward Migkas-like axes; cone binning; noise; Gaussian smoothing) can be made to *look* like published anisotropy maps. At \(\kappa\sim\kappa_{\mathrm{natural}}\) on the baseline decreasing-\(\rho\) model, half peak-to-peak sky variation of order \(10\%\)–\(25\%\) is available from the tidal quadrupole. Direction can be rotated by hand to \((l,b)\sim(280^\circ,-20^\circ)\) or joint \(\sim(303^\circ,-27^\circ)\). None of that constitutes a likelihood fit.

### 5.2 What fails as a *pure dipole* (not a veto of spherical decreasing-\(\rho\))

**Important:** the items below constrain matching Migkas as a clean full-sky *dipole*, not the viability of spherical variable-density cosmology. Spherical decreasing-\(\rho\) remains LIVE; the predicted sky is quadrupole-led. For the baseline decreasing-\(\rho\) Gauss model + beyond-horizon CoE:

1. Full-sky orthogonal anisotropy is a **leading quadrupole** (\(C_2\approx 2\alpha\kappa\)), not Migkas’s reported dipole. (Uniform \(\rho\) special case: monopole only—§4.1.)
2. Pure-dipole template fits on full sky at Migkas depths have \(R^2\approx 0\).
3. Forcing a true dipole coefficient of Migkas size at \(d\ll D\) requires absurd \(\kappa\sim 10^4\)–\(10^5\) even with enclosed-mass \(\kappa_{\mathrm{natural}}\) language.
4. Hemisphere / cone selection can manufacture fake fitted \(A\sim 0.1\)–\(0.4\) from a pure \(P_2\) sky under the baseline model—systematic channel, not a clean CP probe.

### 5.3 Verdict: spherical decreasing-\(\rho\) LIVE; CONDITIONAL vs Migkas as a pure dipole

| Question | Answer |
|----------|--------|
| Spherical decreasing-\(\rho\) cosmology viable? | **LIVE / GO** — amplitude natural; Migkas-like angular curve exists |
| Amplitude at natural \(\kappa\) (decreasing \(\rho\))? | **GO** — \(\kappa_{\mathrm{natural}}=0.15\); \(\sim 10\%\) variation natural |
| Predicted full-sky multipole at \(d\ll D\)? | **Quadrupole-led** (\(C_2\approx 2\alpha\kappa\); \(C_1\approx 0\)) |
| Time integration rescue to tiny \(\kappa\)? | **No** — \(\mathcal{O}(1)\) only |
| Clean full-sky Migkas *dipole* \((l,b,A)\) from baseline? | **CONDITIONAL NO** at \(d\ll D\) — wrong template; not a veto of the spherical model |
| Uniform \(\rho\) as working model? | **Rejected** — monopole only |
| Fake fitted \(A\) under selection? | **Yes** — quantified; not a smoking gun |
| Dipole-only MCMC on baseline spherical model? | **Do not run** — fits the wrong template |
| Right next test vs Migkas? | **\(P_2\) / quadrupole comparison** at \(\kappa\sim\kappa_{\mathrm{natural}}\) |
| Smoking CP via Migkas dipole MCMC on spherical? | **No claim** — over-broad “NO-GO for spherical” language retired; morphology caveat only |
| Non-spherical mass / lopsided blast? | **Open** — not yet computed |

**Language fix (binding).** Decreasing / variable \(\rho\) *fixes* the uniform-spherical monopole-only problem. **Spherical + decreasing \(\rho\) is LIVE.** The remaining caveat is *only* that at \(d\ll D\) the full-sky Legendre dipole \(C_1\approx 0\) and the leading multipole is quadrupole. “NO” applies to a clean full-sky Migkas \((l,b,A)\) *dipole* MCMC—**not** to spherical variable-density cosmology as a whole. Preferred summary: *Spherical decreasing-\(\rho\) is viable; predicted sky is quadrupole-led; CONDITIONAL vs Migkas as a pure dipole; do not run dipole-only MCMC; \(P_2\) comparison is the right test.* Selection / \(\phi\)-uniform leakage / \(d/D\gtrsim 1\) remain systematic or ontology channels, not licenses to claim a smoked CP. Open for non-spherical mass (dipole mass moment, lopsided blast).

Bulk-flow degeneracy (Migkas 2021) and FLAMINGO-style nulls further weaken any claim that matching a Migkas map equals locating a CoE.

---

## 6. Parameter estimation path

### 6.1 Parameters of interest

| Parameter | Meaning | Status |
|-----------|---------|--------|
| \((l,b)\) | CoE direction on the sky | Unconstrained; Migkas axis is a target, not a posterior |
| \(D\) | Observer–CoE distance | Unconstrained; beyond-horizon prior vs survey-depth conflict |
| \(M_{\mathrm{enc}}\), \(r_s\), \(\kappa\) | Enclosure, decrease scale, coupling | \(\kappa=\Omega_m/2\); baseline \(m(r)\) from (3.13) |
| \(\xi\) / \(H\) | ST similarity / expansion normalization | \(E,\rho_0\) cosmologically undefined; fit \(H\) directly for now |
| \(\alpha_{\mathrm{eff}}\) | Integrated-gravity efficiency | \(\sim 0.66\) (drag) to \(\sim 1.5\) (frozen); \(\mathcal{O}(1)\) |
| \(r_s\) | Decrease scale in (3.13) | Free; sets concentration |

### 6.2 What to do

1. Implement the cleaned forward model (3.11)–(3.12) with time-integrated \(\alpha_{\mathrm{eff}}\) on HEALPix shells at fixed \(d\) or \(z\)-proxy.
2. Regression: \(M=0\Rightarrow\Psi=0\) (homology theorem).
3. Prefer a **quadrupole (\(P_2\)) template** at \(\kappa\sim\kappa_{\mathrm{natural}}\) on the baseline decreasing-\(\rho\) forward model for Migkas-depth comparisons; report \(R^2\) / \(\chi^2\), not only best-fit amplitude.
4. Quantify footprint-induced dipole leakage from a pure \(P_2\) sky under Migkas-like masks (upper bound on fake \(A\)).
5. Optionally explore non-spherical / dipolar mass moments as a separate forward model.
6. Parallel CMB adversary calculation for the same \((l,b,D,\kappa)\) (Gate E).

### 6.3 What NOT to do yet

- **Do not** run a 3-parameter \((l,b,A)\) **dipole** MCMC on the baseline decreasing-\(\rho\) forward model. The model’s leading beyond-horizon pattern is quadrupole; a dipole MCMC fits the wrong template.
- Do not claim visual Mollweide agreement as parameter constraints.
- Do not quote toy-sheet units (\(G=k=M=1\), \(D=2\), \(R=1\)) as cosmological posteriors.
- Do not assume integration solves morphology or amplitude crises beyond \(\mathcal{O}(1)\).
- Do not eliminate dark energy or claim CP falsification until SN+BAO+CMB joint consistency exists.

### 6.4 Falsification criteria (program-level)

**NO-GO if:** (i) any CoE that fits late-time \(\mathcal{O}(10\%)\) anisotropy forces CMB multipoles beyond Planck / Saadeh limits with no viable screening; (ii) no causal GR construction for a beyond-horizon CoE that still sources the local field; (iii) amplitude needed for Migkas destroys BAO/SNIa isotropy; (iv) Migkas-level signal evaporates as scatter/systematics; (v) bulk flow + FLRW wins on evidence.

**GO only if:** a CoE posterior fits late-time anisotropy better than bulk-flow nulls, has an explicit CMB-compatibility mechanism or prediction, and forecasts distinctive eROSITA / Roman / CMB-S4 tests.

---

## 7. Relativistic / light-travel / CMB open problems

This section is intentionally unfinished as physics—not as an outline bullet, but as an honest list of gaps that currently block peer-reviewable center localization.

### 7.1 Light cones and luminosity distance

Newtonian \(v_{\mathrm{los}}/d\) is not an astronomical Hubble diagram. Real observables are redshift \(z\), luminosity distance \(d_L(\hat{n},z)\), and selection functions. Photon travel time means source and observer positions evolve along the null geodesic; “frozen \(t\)” maps are approximations. Aberration, gravitational redshift, and integrated Sachs–Wolfe–like effects have no Newtonian counterpart in the present skeleton.

### 7.2 ST is not a metric cosmology

Classical ST similarity assumes a static uniform ambient medium. Cosmic expansion of the medium itself contradicts that assumption. Mapping \((E,\rho_0,\xi)\) onto FLRW or LTB parameters is undefined. A serious paper needs an explicit dynamical system or metric (LTB-like, blast-wave metric, or other) and matching across horizons—not only \(V=kD\).

### 7.3 Beyond-horizon CoE construction

The ontology places the CoE outside the particle horizon so we do not “see” the center. No Newtonian or GR construction in eqs (3.1)–(3.14) shows how a center outside the particle horizon sets local \(H\) and \(m(r)\) while remaining causally consistent. This is asserted, not derived.

### 7.4 CMB adversary (make-or-break)

Off-center observers in centered cosmologies generically imprint CMB dipoles and higher multipoles (Alnes & Amarzguioui; Kodama et al.). Saadeh et al. (2016) bound global anisotropic expansion at extreme odds. The program’s credibility hinges on computing predicted CMB multipoles for the same \((l,b,D,M_{\mathrm{enc}})\) used in late-time fits. Rhetoric about a center “beyond the horizon” does not waive that calculation. **No such calculation exists yet in this program.**

### 7.5 Double-counting and back-reaction

Bolting or integrating \(GM/r^2\) onto an ST background risks double-counting blast mass-energy and ignores back-reaction on \(H(t)\). A self-consistent ST+gravity similarity solution (or GR analogue) is not in hand.

### 7.6 Redshift-dependent discriminants

Anisotropic expansion from a CoE and a bulk flow diverge with redshift. Predictions for eROSITA clusters and SN Ia anisotropic Hubble diagrams are required to break the Migkas 2021 degeneracy. Those forecasts are not yet delivered.

---

## 8. Discussion: can this smoke the cosmological principle?

### 8.1 What would smoking the CP require?

To smoke the cosmological principle with this program one would need, at minimum: (1) a reproducible late-time anisotropy not explained by bulk flow or systematics; (2) a centered dynamical model that predicts that anisotropy’s multipole structure and redshift dependence; (3) CMB compatibility or a sharp, tested mechanism for near-isotropy at \(z\sim 1100\); (4) evidence ratios favoring a CoE over FLRW + bulk flow. None of (2)–(4) are satisfied today. Item (1) remains contested.

### 8.2 Where the present model stands

**Strengths.** Clear ontology; cleaned homology theorem preventing a common frame error; natural \(\kappa\) under enclosed-mass; time-integrated gravity at \(\mathcal{O}(1)\) relative to bolt-on; quantitative multipole and footprint-leakage results; falsification gates stated up front.

**Caveats (not a spherical veto).** Baseline spherical decreasing-\(\rho\) is **LIVE** (amplitude GO; angular curve exists) but does **not** produce a clean full-sky Migkas *dipole* at beyond-horizon \(D\)—the sky is quadrupole-led. Migkas “matches” remain visual; parameters unconstrained; ST≠cosmology; CMB adversary uncomputed; Copernican / causal cost of a beyond-horizon center unpaid.

**Honest answer to the section title.** Not yet—and **not** because spherical decreasing-\(\rho\) is dead. Amplitude is no longer absurd; the morphology caveat is against treating Migkas as a pure dipole MCMC target. Smoking the CP remains the *goal* of the program, not an accomplished result of this draft. The right next step on the spherical channel is a **\(P_2\) comparison** at \(\kappa_{\mathrm{natural}}\) (do not run dipole-only MCMC). Paths that remain open: that quadrupole-template test; non-spherical mass; selection-aware reanalysis of Migkas; and—decisively—the CMB calculation that may kill the idea.

### 8.3 Relation to dark energy

The unfinished draft suggested the framework might offer an alternative to ΛCDM for accelerated expansion. This draft does **not** endorse that oversell. Acceleration, SN+BAO joint fits, and dark-energy elimination are out of scope until the anisotropy channel itself is viable and CMB-compatible.

---

## 9. Conclusions and next work

### 9.1 Conclusions

1. Working hypothesis: real CoE beyond the horizon; off-center observer; ST blast-wave expansion plus central gravity from enclosed cosmic mass.
2. Pure homologous ST ⇒ isotropic \(H_{\mathrm{obs}}\) for any offset (eqs 3.6–3.8). The unfinished draft’s claim that offset alone anisotropizes \(H_{\mathrm{obs}}\) is false.
3. Gravity must be time-integrated along ST trajectories; bolt-on is a shortcut. Integration changes amplitudes by \(\mathcal{O}(1)\) only.
4. Baseline: **decreasing \(\rho(r)\)** (eq. 3.13) + Gauss \(g=-G m(r)/r^2\); \(M_{\mathrm{enc}}\) toward CoE; \(\kappa_{\mathrm{natural}}=0.15\). Decreasing density restores anisotropy; uniform \(\rho\) is rejected (monopole only).
5. **Spherical decreasing-\(\rho\) is LIVE.** Beyond-horizon \(D\): leading **quadrupole**, full-sky \(C_1\approx 0\). Fitted Migkas dipoles can appear from selection/leakage/\(d\sim D\)—a morphology caveat, not a veto of the spherical model.
6. Relative acceleration vs angle produces a Migkas-relevant angular curve at natural \(\kappa\); multipole honesty says the prediction is quadrupole-led, not a pure dipole.
7. Comparison with Migkas: **CONDITIONAL** vs a pure-dipole template; spherical variable-density cosmology as a whole is **not** NO-GO. Do not run \((l,b,A)\) dipole-only MCMC; **\(P_2\) comparison is the right test.**
8. Parameters largely unconstrained. CMB is the main adversary. The CP is not smoked.

### 9.2 Next work (ordered)

1. Quadrupole (\(P_2\)) forward template at \(\kappa\sim\kappa_{\mathrm{natural}}\) with drag \(\alpha_{\mathrm{eff}}\); compare to Migkas maps/residuals without forcing a dipole.
2. Footprint leakage study: fake \(A\) from pure \(P_2\) under Migkas-like masks (partially started in living notes; fold into any comparison).
3. One-page quantitative CMB adversary memo for fiducial \((D,\kappa_{\mathrm{natural}})\) using Alnes/Kodama-style estimates as a lower bar, then a proper light-cone calculation.
4. Non-spherical / dipolar mass-moment forward model (only if 1–3 still leave a viable path).
5. Redshift-dependent forecasts vs bulk flow for eROSITA / SN Ia.
6. Only then: MCMC on a template that matches the model’s actual multipole content.

### 9.3 Assets

- Old unfinished draft (untouched): `/workspace/cosmology/draft.txt`
- This new draft: `/workspace/cosmology/new_draft.md` and Google Doc (this file)
- Living notes: `phase_a_derivation.md`, `absurdity_check.md`, `enclosed_mass_multipoles.md`, `dipole_fraction.md`, `literature_map.md`, `project_status_development_plan.md`

---

## References

1. Migkas, K., et al. (2020). Probing cosmic isotropy with a new X-ray galaxy cluster sample through the \(L_X\)–\(T\) relation. *Astronomy & Astrophysics*, 636, A15.
2. Migkas, K., et al. (2021). Cosmic anisotropy with ten cluster scaling relations. *Astronomy & Astrophysics*, 649, A151.
3. Sedov, L. I. (1946). Propagation of strong shock waves.
4. Taylor, G. I. (1950). The formation of a blast wave by a very intense explosion. *Proc. R. Soc. A*, 201(1065), 159–174.
5. Saadeh, D., Feeney, S. M., Pontzen, A., Peiris, H. V., & McEwen, J. D. (2016). How isotropic is the Universe? *Physical Review Letters*, 117, 131302.
6. Planck Collaboration (2020). Planck 2018 results. VII. Isotropy and statistics of the CMB. *Astronomy & Astrophysics*, 641, A7.
7. Alnes, H., & Amarzguioui, M. (2006). CMB anisotropies seen by an off-center observer in a spherically symmetric inhomogeneous universe. *Physical Review D* (astro-ph/0607334).
8. Kodama, H., et al. (2010). Analytic formulae for CMB anisotropies in general spherically symmetric spacetimes. arXiv:1004.3089.
9. Watkins, R., Feldman, H. A., & Hudson, M. J. (2009). Consistently large cosmic flows on scales of \(100\,h^{-1}\,\mathrm{Mpc}\). *Monthly Notices of the Royal Astronomical Society*.
10. Kashlinsky, A., et al. (2008+). Cluster kinetic SZ bulk-flow claims (contested systematics). *Astrophysical Journal*.
11. FLAMINGO / related A&A (2025) analyses of scaling relations as isotropy tracers — simulation reassessment of apparent cluster anisotropies.
12. Pandya, A., et al. (2024). Related cluster \(\sigma_v\) scaling dipole claims.
13. Boubel, P., et al. (2025). Tully–Fisher \(H_0\) variation constraints.
14. Living analysis notes (2026-09-17 CT): Phase A derivation; absurdity check; enclosed-mass multipoles; dipole fraction vs \(d/D\); literature map; project status & development plan.

---

*End of new draft. Supersedes unfinished draft prose claims where they conflict; retains ST + CoE ontology as working hypothesis under active skeptical test.*
