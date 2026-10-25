# Phase A — Clean derivation (eqs 3.1–3.14)

**Sedov–Taylor anisotropic Hubble / Chris — Cosmology agent**  
**Date:** 2026-09-17 (CT)  
**Scope:** Newtonian kinematics of an off-center observer in a centered Sedov–Taylor (ST) velocity field, then central gravity. Consistent notation; skeptical flags where the unfinished draft is hand-wavy or wrong.  
**Sources:** `/workspace/cosmology/draft.txt` §3; toy sheet parameters (k, G, M, D, R); Status doc Phase A.

---

## 0. Framing (binding)

Ontology: the Universe has a **real geometric center of expansion (CoE)**, argued to lie beyond the observable horizon. We are an **off-center observer** at distance \(D\) from the CoE. Expansion is modeled first as a **Sedov–Taylor blast-wave** (dust, self-similar) plus a **central gravitational mass** \(m(r)\). Goal of the larger program: locate the CoE (direction + \(D\)) or kill the idea (CMB adversary).

This note does **not** claim a GR metric, CMB compatibility, or a Migkas fit. It only cleans the Newtonian skeleton corresponding to draft eqs (3.1)–(3.14).

---

## 1. Notation glossary (canonical)

| Symbol | Meaning | Notes |
|--------|---------|-------|
| CoE, \(O\) | Center of expansion / “center of the Universe” | Origin of coordinates |
| \(\vec{D}\), \(D=\|\vec{D}\|\) | Observer (Earth) position relative to CoE | Preferred axis \(\hat{D}=\vec{D}/D\) |
| \(\vec{r}\), \(r=\|\vec{r}\|\) | Source position relative to CoE | **Not** observer–source distance |
| \(\vec{d}=\vec{r}-\vec{D}\), \(d=\|\vec{d}\|\) | Observer→source separation | Draft “distance from observer” |
| \(\hat{n}=\vec{d}/d\) | Line-of-sight unit vector (observer→source) | Sky direction |
| \(\theta\) | Angle between \(\hat{n}\) and \(+\hat{D}\) (toward CoE from observer? see §3) | **Must state convention**; draft uses Earth-centric angle on a circle of radius \(R\) |
| \(R_{\mathrm{ST}}(t)\) | ST blast-front radius | Draft used \(R(t)\) — here \(R_{\mathrm{ST}}\) to avoid clash with toy \(R\) |
| \(R\) | Toy: fixed observer–source distance on a sky circle | Spreadsheet / demos: \(R=1\) (code units) |
| \(\xi\) | ST similarity coefficient | Draft eq (3.1); sheet sometimes calls this \(k\) |
| \(k\) | Homologous expansion rate, \(k\equiv H(t)=\dot{R}_{\mathrm{ST}}/R_{\mathrm{ST}}\) | Sheet “Taylor Sedov” \(k\); **prefer \(H(t)\) or \(\xi\) carefully** |
| \(E\), \(\rho_0\) | Blast energy, ambient density | Cosmological meaning **undefined** (gap) |
| \(H(t)\), \(H_0\) | Expansion rate \(H=\dot{R}_{\mathrm{ST}}/R_{\mathrm{ST}}\); fiducial isotropic value | \(H=2/(5t)\) in classical ST |
| \(\vec{v}(\vec{r})\) | Absolute velocity in CoE frame | |
| \(\vec{v}_{\mathrm{rel}}\) | Velocity of source relative to observer | What enters recession |
| \(v_{\mathrm{los}}\) | Line-of-sight recession speed | \(\vec{v}_{\mathrm{rel}}\cdot\hat{n}\) |
| \(H_{\mathrm{obs}}(\hat{n},d)\) | \(v_{\mathrm{los}}/d\) | Apparent Hubble parameter |
| \(M\), \(m(r)\), \(r_s\) (or \(L\)) | Total / enclosed mass; scale radius | Sheet: \(L=0.5\), \(M\sim1\) |
| \(G\) | Newton constant | |
| \(\vec{g}(r)\), \(\vec{u}\) | Gravity acceleration; induced peculiar velocity | |
| \(\delta v\), \(\Psi\) | Peculiar LOS correction; residual \(H_{\mathrm{obs}}-H_0\) | Draft (3.15) |

**Rule:** Never use the same letter for CoE-radius and observer–source distance. Draft \(r\) vs \(d\) is kept; blast front is \(R_{\mathrm{ST}}\); toy circle radius is \(R\).

---

## 2. Assumptions (Phase A)

1. **Newtonian kinematics** in Euclidean 3-space (2D demos are equatorial slices).
2. **Instantaneous** velocities (no light travel time; positions frozen at one \(t\)).
3. **Dust-like** medium; classical ST similarity for the *background* expansion of the blast front.
4. Interior velocity field taken as **homologous** about the CoE: \(\vec{v}=H(t)\,\vec{r}\) (standard ST interior scaling near the self-similar solution; see caveats).
5. Observer and sources are **comoving tracers** of that field (plus optional peculiar \(\vec{u}\) from gravity).
6. Central mass is either a point mass \(M\) or an extended spherical \(m(r)\) about the CoE.
7. “Observed” distance for \(H_{\mathrm{obs}}\) is the Euclidean \(d\) (not luminosity distance).

**Explicitly out of scope here:** GR, redshift composition, CMB last-scattering, beyond-horizon matching, Migkas likelihood.

---

## 3. Geometric setup

Place CoE at the origin. Observer at \(\vec{D}=(D,0,0)\). Source at \(\vec{r}=(P_x,P_y,P_z)\).

**(3.4′)** Distance of source from CoE:
\[
r = \|\vec{r}\| = \sqrt{P_x^2+P_y^2+P_z^2}.
\]

**(3.5′)** Distance of source from observer:
\[
d = \|\vec{r}-\vec{D}\| = \sqrt{(P_x-D)^2+P_y^2+P_z^2}.
\]

**Sky angle (canonical):** define \(\theta\) by
\[
\cos\theta \equiv \hat{n}\cdot\hat{D},
\]
so \(\theta=0\) is the direction **from the observer toward the CoE** along \(+\hat{D}\) wait — careful:

- \(\hat{D}=\vec{D}/D\) points **from CoE to observer**.
- Direction **from observer toward CoE** is \(-\hat{D}\).
- Draft demos place sources on a circle \(\vec{r}=\vec{D}+R(\cos\phi,\sin\phi)\) with \(\phi\) an Earth-centric polar angle from \(+\hat{x}\) (same as \(\hat{D}\)). Thus \(\phi=0\) is **away from CoE** (anti-center), \(\phi=\pi\) is **toward CoE**.

**Convention adopted here (match demos):** \(\phi\) = Earth-centric angle from \(+\hat{D}\) (anti-center). Then
\[
\hat{n}=(\cos\phi,\sin\phi),\qquad \hat{n}\cdot\hat{D}=\cos\phi.
\]
When quoting “\(\theta\) toward the center,” use \(\theta=\pi-\phi\) or state explicitly. **Draft mixes \(\theta\) without a fixed convention — flagged as inconsistency.**

---

## 4. Sedov–Taylor background (cleaned 3.1–3.3)

Classical ST blast front in uniform density \(\rho_0\) with energy \(E\):

**(3.1)**
\[
R_{\mathrm{ST}}(t)=\xi\left(\frac{E}{\rho_0}\right)^{1/5}t^{2/5}.
\]

**(3.2)**
\[
v_{\mathrm{front}}(t)=\frac{dR_{\mathrm{ST}}}{dt}=\frac{2}{5}\frac{R_{\mathrm{ST}}}{t}.
\]

**(3.3)** Dimensional / kinematic analogy (not a derivation of FLRW):
\[
H(t)\equiv\frac{\dot{R}_{\mathrm{ST}}}{R_{\mathrm{ST}}}=\frac{2}{5t}\qquad\Rightarrow\qquad v_{\mathrm{front}}=H(t)\,R_{\mathrm{ST}}.
\]
Interior homologous model (Phase A working assumption):
\[
\vec{v}(\vec{r},t)=H(t)\,\vec{r}.
\]
Then relative to any comoving observer, \(v\propto d\) (Hubble-like). **Cosmological identification of \((E,\rho_0,\xi)\) remains undefined.**

---

## 5. Off-center observer: pure ST is isotropic (cleaned 3.6–3.8)

Observer velocity: \(\vec{v}_O=H\vec{D}\).  
Source velocity: \(\vec{v}_S=H\vec{r}\).

**(3.6′)** Relative velocity:
\[
\vec{v}_{\mathrm{rel}}=\vec{v}_S-\vec{v}_O=H(\vec{r}-\vec{D})=H\,\vec{d}.
\]

**(3.7′)** Line-of-sight recession:
\[
v_{\mathrm{los}}=\vec{v}_{\mathrm{rel}}\cdot\hat{n}=H\,d.
\]

**(3.8′)** Apparent Hubble parameter:
\[
H_{\mathrm{obs}}=\frac{v_{\mathrm{los}}}{d}=H(t).
\]

### Critical result (draft contradiction)

For **any** observer offset \(D\), a purely homologous field \(\vec{v}=H\vec{r}\) gives **\(H_{\mathrm{obs}}=H\) independent of direction**. Geometry alone does **not** produce Migkas-like anisotropy. This is the classical reason “\(v=Hd\) does not imply a center.”

The unfinished draft states after (3.8) that the expression “depends on … offset \(D\), and therefore introduces anisotropy.” That claim is **false** if (3.6)–(3.7) are applied correctly (subtract observer motion). Anisotropy requires **breaking homology** (radial ST structure, gravity, selection on non-homologous shells, or relativistic effects).

**Failure mode that falsely creates a dipole:** using absolute CoE-frame speed \(H\vec{r}\cdot\hat{n}\) without subtracting \(\vec{v}_O\). Then on a circle of radius \(R\) about the observer,
\[
v_{\mathrm{los}}^{\mathrm{(wrong)}}=H(D\cos\phi+R),\qquad
H_{\mathrm{obs}}^{\mathrm{(wrong)}}=H\left(1+\frac{D}{R}\cos\phi\right).
\]
Do **not** use this as physics; it is a frame error. (Useful only as a cautionary cross-check against buggy code.)

---

## 6. Central gravity (cleaned attempt at 3.9–3.12)

**(3.9′)** Spherical gravity from enclosed mass:
\[
\vec{g}(\vec{r})=-\frac{G m(r)}{r^2}\,\hat{r},\qquad \hat{r}=\vec{r}/r.
\]
Point mass: \(m(r)=M\).

### Peculiar velocity is underdetermined in Phase A

Acceleration is not velocity. Draft (3.10) writes
\[
u(r)\sim g(r)\,t_H=\frac{GM}{r^2 H},
\]
i.e. an impulse over a Hubble time with \(t_H\sim 1/H\). This is an **order-of-magnitude ansatz**, not a solution of the equations of motion inside an ST blast. Accept it only as a **placeholder** \(\vec{u}(\vec{r})=-u(r)\hat{r}\) with
\[
u(r)=\alpha\,\frac{G m(r)}{r^2 H},
\]
where \(\alpha\sim\mathcal{O}(1)\) is a free dimensionless fudge until a dynamical derivation exists.

**(3.10′)** Peculiar field (ansatz):
\[
\vec{u}(\vec{r})=-\alpha\frac{G m(r)}{r^2 H}\,\hat{r}.
\]

**(3.11′)** Relative peculiar + ST recession along LOS:
\[
\vec{v}_{\mathrm{rel}}=H\vec{d}+\bigl[\vec{u}(\vec{r})-\vec{u}(\vec{D})\bigr],
\]
\[
v_{\mathrm{los}}=H d+\delta v,\qquad
\delta v=\bigl[\vec{u}(\vec{r})-\vec{u}(\vec{D})\bigr]\cdot\hat{n}.
\]

**(3.12′)** Apparent Hubble:
\[
H_{\mathrm{obs}}(\hat{n},d)=H+\frac{\delta v}{d}.
\]

### Why draft (3.11)–(3.12) are rejected as written

Draft:
\[
v_{\mathrm{los}}(\theta)\approx H_0 d-\frac{GM}{d^2}\cos\theta,\qquad
H_{\mathrm{obs}}(\theta)=H_0-\frac{GM}{d^3}\cos\theta.
\]

Problems:
1. **Dimensions:** \(GM/d^2\) is acceleration, not velocity; the draft dropped the \(1/H\) from (3.10).
2. **Wrong radius:** gravity depends on CoE-radius \(r\), not observer distance \(d\).
3. **Missing differential peculiar:** observer also feels \(\vec{u}(\vec{D})\); only the difference projects.
4. **Angle geometry** unspecified; \(\cos\theta\) is not automatically the projection of \(-\hat{r}\).
5. **Double-counting risk:** ST already responds to mass-energy in a real blast; bolting \(GM/r^2\) onto a dust ST solution is not self-consistent hydrodynamics.

**Small-\(R/D\) dipole sketch (for intuition only):** with \(r\approx D\), \(\vec{u}(\vec{r})\approx\vec{u}(\vec{D})+(\vec{d}\cdot\nabla)\vec{u}\), the leading sky pattern in \(\delta v/d\) is dipolar in \(\hat{n}\cdot\hat{D}\) with amplitude set by \(\partial_r u\) and tidal field \(\sim Gm(D)/D^3\). Amplitude scaling \(\sim \alpha Gm(D)/(D^3 H)\) relative to \(H\) — i.e. a dimensionless \(\varepsilon\sim \alpha (GM)/(D^3 H^2)\). Fitting \(\varepsilon\sim\mathcal{O}(0.1)\) to Migkas without absurd \((M,D)\) is a **Phase C** question, not assumed here.

---

## 7. Extended mass profile (cleaned 3.13–3.14)

Draft / sheet profile (cumulative mass fraction):

**(3.13′)**
\[
\frac{m(r)}{M}=1-e^{-X}\left(1+X+\frac{X^2}{2}\right),\qquad X=\frac{r}{r_s}.
\]
(Sheet used \(L\) for \(r_s\), \(L=0.5\) in code units; normalization \(M\approx1.142\ldots\) appears as a sheet quirk — treat \(M\) as the free total-mass parameter.)

**(3.14′)**
\[
\vec{g}(r)=-\frac{G m(r)}{r^2}\,\hat{r}.
\]
Same peculiar ansatz as (3.10′) with \(m(r)\) replacing \(M\).

This profile is the incomplete-gamma cumulative associated with density \(\rho\propto e^{-r/r_s}\) (specifically \(m(r)\propto\int_0^r s^2 e^{-s/r_s}ds\)). **No derivation from ST hydro is given in the draft** — it is a convenient smooth core + asymptotic point-mass envelope.

---

## 8. Residual map (optional 3.15)

**(3.15′)**
\[
\Psi(\hat{n},d)=H_{\mathrm{obs}}(\hat{n},d)-H.
\]
For pure ST, \(\Psi=0\). With gravity ansatz, \(\Psi=\delta v/d\).

---

## 9. Clean ordered equation list (Phase A canonical)

| Tag | Equation | Status |
|-----|----------|--------|
| **3.1** | \(R_{\mathrm{ST}}(t)=\xi(E/\rho_0)^{1/5}t^{2/5}\) | Classical ST; \(E,\rho_0\) cosmologically undefined |
| **3.2** | \(v_{\mathrm{front}}=(2/5)R_{\mathrm{ST}}/t\) | OK |
| **3.3** | \(H(t)=2/(5t)\), \(v_{\mathrm{front}}=H R_{\mathrm{ST}}\) | OK as ST algebra |
| **3.4** | \(r=\|\vec{r}\|\) | Geometry |
| **3.5** | \(d=\|\vec{r}-\vec{D}\|\) | Geometry |
| **3.6** | \(\vec{v}=H\vec{r}\); \(\vec{v}_{\mathrm{rel}}=H\vec{d}\) | Homologous ST |
| **3.7** | \(v_{\mathrm{los}}=H d\) (pure ST) | OK |
| **3.8** | \(H_{\mathrm{obs}}=H\) (pure ST; **isotropic**) | **Corrects draft prose** |
| **3.9** | \(\vec{g}=-G m(r)r^{-2}\hat{r}\) | OK Newtonian |
| **3.10** | \(\vec{u}=-\alpha[G m(r)/(r^2 H)]\hat{r}\) | **Ansatz**; \(\alpha\) free |
| **3.11** | \(v_{\mathrm{los}}=H d+[\vec{u}(\vec{r})-\vec{u}(\vec{D})]\cdot\hat{n}\) | Replaces draft (3.11) |
| **3.12** | \(H_{\mathrm{obs}}=H+\delta v/d\) | Replaces draft (3.12) |
| **3.13** | \(m(r)/M=1-e^{-X}(1+X+X^2/2)\), \(X=r/r_s\) | Phenomenological |
| **3.14** | \(g=G m(r)/r^2\) (magnitude) | OK |
| **3.15** | \(\Psi=H_{\mathrm{obs}}-H\) | Optional residual |

Draft tags 3.1–3.14 are **retained**; contents of **3.8, 3.11, 3.12** are corrected. Do not renumber yet (keeps Status/lit cross-refs stable).

---

## 10. Toy sheet / demo correspondence

Sheet / demo parameters (code units): \(k\sim1\), \(G\sim1\), \(M\sim1\), \(D=2\), \(R=1\), \(r_s=L=0.5\).

Recommended cross-check (Phase A acceptance test):
1. **Gravity off** (\(\alpha=0\) or \(M=0\)): recession vs \(\phi\) on the observer-centered circle of radius \(R\) must be **flat** at \(v_{\mathrm{los}}=k R\) (or \(H R\)). If the sheet shows angular variation with \(M=0\), the sheet uses the frame-error formula — fix it.
2. **Gravity on**: implement (3.11′) with \(\vec{u}\) from \(m(r)\); compare shape (not absolute scale) to draft `recession_v` column only after documenting \(\alpha\) and whether \(\vec{u}(\vec{D})\) is subtracted.
3. Demographics of draft numerical curve (min near \(\phi\sim110^\circ\)–\(120^\circ\), max near \(180^\circ\)): interpret only after (1)–(2) agree.

---

## 11. Gaps / inconsistencies found

1. **Homology ⇒ isotropy:** Draft §3.2 claims \(D\) alone anisotropizes \(H_{\mathrm{obs}}\). False for \(\vec{v}=H\vec{r}\). Anisotropy must come from gravity (ansatz), non-homologous ST interior, or GR/lightcone effects.
2. **Eqs 3.11–3.12:** Dimensionally inconsistent; use \(d\) instead of \(r\); omit observer peculiar; drop \(1/H\).
3. **\(k\) vs \(\xi\) vs \(H\):** Sheet \(k\), ST \(\xi\), Hubble \(H\) conflated. Canonical: \(\xi\) in (3.1); \(H=\dot{R}_{\mathrm{ST}}/R_{\mathrm{ST}}\); reserve \(k\) only when matching the sheet literally.
4. **\(\theta\) convention:** Earth-centric \(\phi\) from anti-center vs angle toward CoE never fixed in prose; demos use \(\phi\) from \(+\hat{x}\).
5. **\(R(t)\) vs toy \(R\):** Same letter, different meanings — resolved here as \(R_{\mathrm{ST}}\) vs \(R\).
6. **ST ≠ cosmology:** \(E,\rho_0,\xi\) have no FLRW map; blast-wave ambient density contradicts cosmic expansion of the medium itself.
7. **Gravity bolt-on:** Not a self-consistent ST+gravity similarity solution; \(\alpha\) and \(t_H\) are free.
8. **Beyond-horizon CoE:** Asserted in framing; **no Newtonian or GR construction** in eqs 3.1–3.14 shows how a center outside the particle horizon sets local \(H\) and \(m(r)\).
9. **Migkas “match”:** Visual Mollweide rotations only; not constrained by (3.12′).
10. **Mass profile (3.13):** Phenomenological; sheet \(M=1.142\ldots\) unexplained in draft prose.
11. **2D vs 3D:** Formalism is 3D; most demos are planar rings extruded in azimuth — quadrupole content may be misrepresented.
12. **Newtonian-only:** No redshift, no DL vs \(d\), no aberration, no integrated Sachs–Wolfe analog — insufficient for CMB or SN Ia claims.

---

## 12. Newtonian limitations (skeptical summary)

- Instantaneous Euclidean \(v_{\mathrm{los}}/d\) is **not** an astronomical Hubble diagram (needs \(z\), \(d_L\), selection).
- ST similarity assumes a **static uniform** ambient medium — cosmologically false.
- Central gravity as peculiar impulse over \(1/H\) ignores cosmic expansion work, continuity, and back-reaction on \(H(t)\).
- Off-center observers in centered GR models (LTB etc.) are **CMB-dipole constrained** (Alnes & Amarzguioui; Kodama et al.); Phase A cannot waive that.
- Saadeh et al. (2016) anisotropic-expansion limits are outside Newtonian scope but bound the program’s viability.

---

## 13. What is still undefined (must resolve before MCMC)

| Item | Needed for |
|------|------------|
| Map \((\xi,E,\rho_0)\to H_0\) or drop ST energy params and fit \(H\) directly | Normalization |
| Dynamical \(\vec{u}\) (or fit \(\alpha\), or replace with tidal template) | Amplitude of \(\Psi\) |
| Angle convention + HEALPix \(\hat{n}\) relative to CoE | Direction fit \((l,b)\) |
| Relation \(d\) vs redshift / \(d_L(\hat{n})\) | Migkas DL interpretation |
| Prior ranges for \((D, M, r_s, \alpha)\) in physical units | MCMC |
| CMB prediction pipeline | Go/no-go |

---

## 14. Recommended next micro-step (toward MCMC / Migkas)

**Micro-step A→C (minimal forward model):**
1. Implement (3.11′)–(3.12′) on a HEALPix grid for fixed shell \(d\) (or fixed \(z\) proxy), with parameters \(\{D, M, r_s, \alpha, H\}\) and CoE direction \((l,b)\).
2. Confirm \(M=0\Rightarrow\Psi=0\) (regression test against homology theorem).
3. Project dipole amplitude \(A_1(d)=\sqrt{(3/4\pi)}\,|a_{10}|\) (align CoE to \(z\)-axis) vs \(D,M\).
4. Ask: is there any **non-absurd** \((GM)/(D^3 H^2)\) giving \(\sim5\)–\(15\%\) variation in \(H_{\mathrm{obs}}\) or \(d_L\) on Migkas’s effective depth — **before** full MCMC?
5. If yes, run a **3-parameter** MCMC on \((l,b,A)\) matching a Migkas dipole template; if no, document amplitude failure (informs go/no-go).

Parallel (do not skip): one-page **CMB adversary memo** (Status action #2).

---

## 15. Deliverable checklist

- [x] Consistent notation glossary  
- [x] Cleaned eqs 3.1–3.14 (with corrected 3.8, 3.11, 3.12)  
- [x] Assumptions + derivation steps  
- [x] Gaps / Newtonian limits flagged  
- [x] Next micro-step toward Migkas/MCMC  
- [ ] Spreadsheet numerical agreement test (pending code/sheet fix if \(M=0\) is not flat)

---

*End of Phase A note.*
