# CoE letter — every calculation

**Source:** Chris Alameddin, *Centre-of-Expansion shell cosmology: from Migkas anisotropy to a flat CMB* (PDF 29 Sep 2026 / `letter.tex` 1 Oct 2026), plus supporting `kill_function_from_geometry.md`, `_draft_from_chris.md`, and cycle scripts under `cmb_reconcile/`.

**Purpose:** One place for every display/inline equation, fitted number, shell rule, and multipole/kill diagnostic quoted or used in the letter flow. Not a claim that CMB isotropy at \(\sim 10^{-5}\) has been achieved.

---

## 0. Stance (problem)

CMB isotropy (Planck):

\[
\frac{\Delta T}{T} \sim 10^{-5}
\]

*Plain language:* After the dipole, the last-scattering sky is flat to about one part in \(10^{5}\).

Method chain (Migkas fixes geometry; CMB is an external check):

\[
\text{Migkas data} \;\rightarrow\; \text{fixed CoE+shell} \;\rightarrow\; \text{independent CMB stress test}
\]

*Plain language:* Do not refit the shell to the CMB; predict a multipole from the Migkas geometry and ask whether it can be cancelled.

---

## 1. Best-fit geometry and kinematic picture

### 1.1 Local Hubble velocity

\[
\vec{V}_H = \vec{V}_g - \vec{V}_O
\]

*Plain language:* Observed directional Hubble flow is free-fall toward the CoE minus our own motion.

### 1.2 Vision radius and shell rule

Vision / diagnostic emission radius:

\[
R \approx 14.3\,\mathrm{Gpc}
\]

(cycle tables use \(R = 14.30048287260589\,\mathrm{Gpc}\).)

Shell edges relative to observer distance \(D\) (headline Newtonian \(\Omega_s=0.049\) rule used throughout cycles 4–6):

\[
R_{\mathrm{in}} = D - k_{\mathrm{in}} R,\qquad
R_{\mathrm{out}} = D + k_{\mathrm{out}} R
\]

with

\[
k_{\mathrm{in}} \approx 1.028,\qquad
k_{\mathrm{out}} \approx 1.475
\]

(precise stored values: \(k_{\mathrm{in}}=1.0276268339399446\), \(k_{\mathrm{out}}=1.4754794866889613\)).

Shell-rule lower bound (observer mid-shell; vision sphere must not poke through the inner edge toward the CoE in a way that breaks the fit):

\[
D \gtrsim k_{\mathrm{in}} R \approx 14.7\,\mathrm{Gpc}
\]

*Plain language:* \(D\) must sit above \(\sim 14.7\,\mathrm{Gpc}\) so US remains between \(R_{\mathrm{in}}\) and \(R_{\mathrm{out}}\).

### 1.3 Preferred relativistic headline (letter)

\[
D \approx 32.9\,\mathrm{Gpc},\quad
R_{\mathrm{in}} \approx 18.1\,\mathrm{Gpc},\quad
R_{\mathrm{out}} \approx 49.7\,\mathrm{Gpc}
\]

*Plain language:* Relativistic LTB fit preferred in the letter’s geometry section (observer mid-shell).

### 1.4 Newtonian \(\Omega_s=0.049\) solved geometry (letter §4 / Table 1 column)

\[
\begin{aligned}
D &\approx 49.17\,\mathrm{Gpc},\\
R_{\mathrm{in}} &\approx 34.48\,\mathrm{Gpc},\\
R_{\mathrm{out}} &\approx 70.27\,\mathrm{Gpc},\\
\Omega_s &= 0.049,\\
\ell_0 &= 13^{\circ}.
\end{aligned}
\]

High-precision stored solve (`newt_ltbgeom_A_vs_ell.json`):

| Quantity | Value |
|----------|------:|
| \(D\) | \(49.17456190433319\,\mathrm{Gpc}\) |
| \(R_{\mathrm{in}}\) | \(34.47900196614479\,\mathrm{Gpc}\) |
| \(R_{\mathrm{out}}\) | \(70.27463103261002\,\mathrm{Gpc}\) |
| \(\Omega_s\) | \(0.049\) |
| \(\ell_0\) | \(13.0^{\circ}\) |
| \(R^2\) vs Migkas | \(0.893\) |
| \(\langle H\rangle\) | \(71.25\,\mathrm{km\,s^{-1}\,Mpc^{-1}}\) |
| \(R_s\) (vision) | \(14.300\,\mathrm{Gpc}\) |
| \(m_{cc}\) (core / \(M_{\mathrm{crit}}\)) | \(0.001714\) |
| age | \(13.886\,\mathrm{Gyr}\) |

*Plain language:* Full Newtonian shell+core numbers that match Migkas \(A(\ell)\) at \(\Omega_s=0.049\).

### 1.5 Mean Hubble conversion (cycle convention)

\[
\langle H\rangle = 71.6\,\mathrm{km\,s^{-1}\,Mpc^{-1}}
\]

(used in cycle 4/5/6 CMB diagnostics; letter fit store uses \(\sim 71.25\) for the Newtonian solve above).

---

## 2. Table 1 — Newtonian and relativistic fits to Migkas \(A(\ell)\)

| | Newt. \(\Omega_s\) free | Newt. \(\Omega_s=0.049\) | Rel. LTB \(\Omega_s=0.049\) |
|--|--:|--:|--:|
| \(D\) (Gpc) | \(93.7\) | \(49.2\) | \(32.9\) |
| \(\Omega_s\) | \(1.01\) | \(0.049\) | \(0.049\) |
| \(\ell_0\) (deg) | \(12.8\) | \(13.0\) | \(13.5\) |
| \(R^2\) vs Migkas | \(0.909\) | \(0.893\) | \(0.885\) |
| \(R^2\) vs \(A_F\) | \(0.987\) | \(0.977\) | \(0.963\) |

*Plain language:* Three headline fits; relativistic column is the letter’s preferred \(D\approx 32.9\,\mathrm{Gpc}\).

Cycle-1 convention note on the same three \(D\) points: model \(\langle H\rangle \approx 71.79 / 71.61 / 71.63\) (free / omb / rel).

Redshift slices used for consistency (Fig. 3): \(0.03\)–\(0.07\), \(z>0.07\), and all-\(z\); model recomputed at each slice’s median distance with fit parameters held fixed. Also shown: \(D\)–core-mass degeneracy.

---

## 3. Einstein and Friedmann equations (CoE+shell / LTB)

LTB origin is the **CoE**, not US. Observer sits at comoving distance \(D\). Sky anisotropy is off-centre sampling of radial structure \(\beta(r)\), **not** a direction-dependent scale factor \(a(t,\theta)\).

### 3.1 Einstein equation

\[
G_{\mu\nu} = 8\pi G\, T_{\mu\nu}.
\]

*Plain language:* Matter–geometry coupling for the CoE-centred metric.

### 3.2 LTB line element

\[
\mathrm{d}s^{2}
=
-\mathrm{d}t^{2}
+
\frac{(R')^{2}}{1+2E(r)}\,\mathrm{d}r^{2}
+
R(t,r)^{2}\,\mathrm{d}\Omega^{2},
\]

with areal radius \(R(t,r)\), \(R'\equiv\partial_r R\), \(\dot{R}\equiv\partial_t R\).

*Plain language:* Concentric dust shells about the CoE; no angular dependence in the metric coefficients.

### 3.3 Friedmann-like shell constraint

\[
\dot{R}^{2} = 2E(r) + \frac{2GM(r)}{R}.
\]

*Plain language:* Each shell’s expansion energy plus self-gravity.

Sign of \(E(r)\):

- \(E=0\): marginally bound dust;
- \(E>0\): expands forever;
- \(E<0\): can recollapse.

### 3.4 Interior mass and radial multiplier \(\beta(r)\)

\[
M(r) = \frac{4\pi}{3}\, R^{3}\, \nu\, \beta(r).
\]

*Plain language:* \(\nu\) sets the background density scale; \(\beta(r)\) is a radial mass multiplier (not a scale factor). \(\beta=1\) recovers critical-density normalisation. In the CoE+shell geometry \(\beta\) is small in the cavity, rises through the thick shell, and returns toward the exterior background beyond \(R_{\mathrm{out}}\).

### 3.5 Matter + radiation + curvature on the same radial scaffolding

\[
H_a^{2}
\propto
\frac{\nu_m \beta_m(r)}{a^{3}}
+
\frac{\nu_r \beta_r(r)}{a^{4}}
-
\frac{k(r)}{a^{2}},
\]

with \(H_a=\dot{a}/a\) and \(a(t,r)\) the local expansion factor of each shell.

Asymptotics:

\[
a \sim t^{1/2}\quad\text{(early, radiation)},\qquad
a \sim t^{2/3}\quad\text{(late, dust)}.
\]

*Plain language:* Standard early/late power laws on each shell; curvature \(k(r)\) follows from \(E(r)\) with no angular metric dependence.

### 3.6 Weak-field matter response (logic-path context)

\[
\mathbf{a} \simeq -\nabla\Phi,
\qquad
\Delta\mathbf{v} \sim -\int \nabla\Phi\,\mathrm{d}t.
\]

*Plain language:* Matter can accumulate large peculiar velocities over cosmic time.

Schematic ISW / Rees–Sciama contribution:

\[
\frac{\Delta T}{T}
\propto
\frac{2}{c^{2}}
\int \dot{\Phi}\,\mathrm{d}\eta.
\]

*Plain language:* Photons see a changing potential differently from accumulated free-fall of matter.

---

## 4. Sky orientation

Sequence (fit parameters fixed; only placement varies):

1. Untilted, no-roll.
2. \(+15^{\circ}\) tilt.
3. Tilt + radial gradient.
4. Tilt + gradient + \(+80.8^{\circ}\) roll about the CoE axis.

Predicted MAX/MIN bands compared to Migkas Fig. 8.

Fitted CoE longitude phases from Table 1: \(\ell_0 \approx 12.8^{\circ}\)–\(13.5^{\circ}\) (Migkas \(A(\ell)\) phase).

*Plain language:* Place the CoE axis on the Galactic sky so model MAX/MIN track Migkas Fig. 8.

---

## 5. Time evolution (donut universe) and volumes

CoE fixed; thick shells expand outward faster than light in a coordinate sense; the sphere of events visible to us grows at \(c\). Recombination is a spherical surface of **fixed physical radius** around the CoE. Observer remains mid-shell between inner/outer edges (S1 and S2).

Hubble-bubble volume about US:

\[
V_{\mathrm{HB}} = 0
\quad\text{until the observer passes } D_{\mathrm{recomb}};
\quad\text{afterward the visible sphere grows at } c.
\]

Useful volume fractions (storyboard language):

\[
\frac{V_{\mathrm{HB}}}{V_{\mathrm{ann}}},\qquad
\frac{V_{\mathrm{HB}}}{V_{\mathrm{S1}}},\qquad
\frac{V_{\mathrm{S2}}}{V_{\mathrm{S1}}}.
\]

*Plain language:* Track how the Hubble bubble fills the shell annulus relative to inner/outer spheres as the geometry evolves.

---

## 6. CMB stress test — Doppler, static, kill function, multipoles

### 6.1 Working stress-test geometry (fixed \(R\))

- \(R \approx 14.3\,\mathrm{Gpc}\)
- \(\Omega_s = 0.049\), \(R_{\mathrm{in}}=D-1.028\,R\), \(R_{\mathrm{out}}=D+1.475\,R\)
- Compact point core; \(M_c\) re-solved at each \(D\) so cluster quadrupole stays at

\[
a_2 = 0.084
\quad\text{at}\quad
R_s = 0.31\,\mathrm{Gpc}
\]

- Preferred diagnostic distances:

\[
D = 14.7,\quad
22.28,\quad
49\,\mathrm{Gpc}
\]

(with static-null \(D_{\mathrm{null}} = 22.27586\,\mathrm{Gpc}\) from `cycle4/partB.json`).

- Angle convention: \(\theta=0^{\circ}\) toward CoE, \(180^{\circ}\) anti-CoE.

### 6.2 Kinematic Doppler (linear convention)

\[
\left(\frac{\Delta T}{T}\right)_{\mathrm{Doppler}}
=
-\,\hat{\mathbf{n}}\cdot\frac{\mathbf{V}_{\mathrm{emit}}-\mathbf{V}_{\mathrm{obs}}}{c}.
\]

*Plain language:* Line-of-sight relative velocity of emitters vs observer, in linear \(v/c\).

### 6.3 Static gravitational term

\[
\left(\frac{\Delta T}{T}\right)_{\mathrm{static}}
=
\frac{\Phi_{\mathrm{emit}}-\Phi_{\mathrm{obs}}}{c^{2}}.
\]

By the shell theorem, on the \(R\)-sphere this depends only on

\[
M_{\mathrm{eff}} = M_{\mathrm{interior}} - M_{\mathrm{cavity}}.
\]

Code form used in `lib4.static_pattern` (point-mass equivalent):

\[
\frac{GM}{c^{2}D}
=
\frac{M_e}{2}\left(\frac{HD}{c}\right)^{2}.
\]

*Plain language:* Static pattern can be nulled when core mass equals cavity mass (\(D\approx 22.28\,\mathrm{Gpc}\)); Doppler cannot.

### 6.4 Dimensionless \(HR/c\)

\[
\frac{HR}{c} \approx 3.42
\quad(\langle H\rangle=71.6,\; R=14.3\,\mathrm{Gpc}).
\]

*Plain language:* Linear \(v/c\) is formally outside \(v\ll c\); fractional Doppler \(/ (HR/c)\) still leaves \(|a_2|\) many orders above \(10^{-5}\).

### 6.5 Fixed-\(R\) stress-test numbers (`dtt.json`)

| \(D\) (Gpc) | Doppler \(a_2\) | Doppler ptp | Static \(a_2\) | \(M_c/M_{\mathrm{crit}}\) | notes |
|--:|--:|--:|--:|--:|--|
| \(14.7\) | \(-0.0387\) | \(0.616\) | \(-0.011\) | \(0.00187\) | \(v_{r,\min}<0\) (turnaround near axis) |
| \(22.28\) | \(-0.411\) | \(3.03\) | \(\approx 0\) | \(0.00193\) | static null |
| \(49\) | \(-0.385\) | \(0.795\) | \(+0.087\) | \(0.00196\) | cavity outweighs core |

Target for \(\ell\ge 2\):

\[
|a_\ell| \lesssim 10^{-5}.
\]

*Plain language:* Fixed-\(R\) Doppler \(a_2\sim 0.04\)–\(0.4\) is thousands to tens of thousands times too large.

### 6.6 Kill function definition

After monopole (mean) subtraction on the Doppler curve:

\[
\mathrm{kill}(\theta)
=
-\Biggl[
\Bigl(\tfrac{\Delta T}{T}\Bigr)_{\mathrm{Doppler}}(\theta)
-
\Bigl\langle\tfrac{\Delta T}{T}\Bigr\rangle_{\mathrm{Doppler}}
\Biggr].
\]

On the full axisymmetric great circle \(\phi\in[0^{\circ},360^{\circ}]\):

\[
\mathrm{kill}(\phi)
=
-\,\Delta T/T_{\mathrm{Doppler}}\big|_{\mathrm{demeaned}},
\qquad
\mathrm{kill}(\phi)=\mathrm{kill}(360^{\circ}-\phi).
\]

Required kill quadrupole:

\[
a_2^{\mathrm{kill}} = -a_2^{\mathrm{Doppler}}.
\]

| \(D\) (Gpc) | \(a_2^{\mathrm{kill}}\) | Peak \(\lvert 100\times\mathrm{kill}\rvert\) (% \({\rm d}E/E\)) |
|--:|--:|--:|
| \(14.7\) | \(+0.0387\) | \(\approx 39\%\) |
| \(22.28\) | \(+0.411\) | \(\approx 263\%\) (turnaround spike; use quadrupole) |
| \(49\) | \(+0.385\) | \(\approx 49\%\) |

*Plain language:* The kill function is the electromagnetic \(\Delta T/T\) pattern (kinematic–Doppler character) needed to cancel the gravitational multipole: opposite sign, same angular shape.

Letter wording: the kill function is minus the Doppler residual as a function of angle from the CoE.

### 6.7 Recombination-track variant (cycle 4)

Instead of fixed \(R=14.3\,\mathrm{Gpc}\):

\[
R_{\mathrm{see}} \approx c\,\Delta t \approx 4.25\,\mathrm{Gpc},
\qquad
f = T_0/3000 \approx 9.08\times 10^{-4}.
\]

Kill \(a_2\) plateaus near \(\sim 0.11\)–\(0.115\) for \(D_{\mathrm{now}}\gtrsim 8\,\mathrm{Gpc}\) — still \(\sim 10^{4}\times\) over \(10^{-5}\).

### 6.8 Legendre multipoles about the CoE axis

Axisymmetric pattern \(f(\theta)\) expanded as

\[
a_\ell
=
\frac{2\ell+1}{2}
\int_{-1}^{1}
f(\mu)\, P_\ell(\mu)\, \mathrm{d}\mu,
\qquad
\mu=\cos\theta.
\]

*Plain language:* Standard Legendre coefficients used for Doppler, static, kill, and partner residuals through \(\ell\le 4\).

---

## 7. Two-emitter construction and partner multipoles

### 7.1 Shared bulk velocity → flat kinematic sky

If opposite plates (or shells) share the observer’s bulk Cartesian velocity \(\mathbf{V}\), relative Doppler vanishes and the kinematic sky is flat **before** gravitational multipoles are added.

Special-relativistic observer boost through an isotropic rest-frame field (plates):

\[
T_{\mathrm{obs}}(\hat{\mathbf{n}})
=
\frac{T_{\mathrm{rest}}}{\gamma\,(1-\boldsymbol\beta\cdot\hat{\mathbf{n}})},
\qquad
\gamma=(1-\beta^{2})^{-1/2}.
\]

Toy: \(\beta_{\mathrm{obs}}=0.01\) → demeaned peak \(|\Delta T/T|\approx 0.020\) (classic smooth dipole).

### 7.2 Concentric-shell Doppler toy

Optional Hubble-like recession (first-order, same linear convention):

\[
\boldsymbol\beta_{\mathrm{em}}
=
\beta_{\mathrm{rec}}\,\frac{r_{\mathrm{hit}}}{R_{\mathrm{out}}}\,\hat{\mathbf{e}}_r,
\qquad
\frac{\Delta T}{T}
\approx
-\,\hat{\mathbf{n}}\cdot(\boldsymbol\beta_{\mathrm{em}}-\boldsymbol\beta_{\mathrm{obs}}).
\]

Comoving observer at radius \(D\): \(\beta_{\mathrm{obs,radial}}=\beta_{\mathrm{rec}}(D/R_{\mathrm{out}})\).

Example geometry: \(R_{\mathrm{in}}=1\), \(R_{\mathrm{out}}=3\), \(D=1.8\) → inner-disk half-angle

\[
\alpha = \arcsin(R_{\mathrm{in}}/D) \approx 33.7^{\circ}.
\]

### 7.3 Shape correlation of radial/shear residuals vs kill

Unit-normalized Pearson correlations (`vs_kill_term_0_360.json`, toy \(\beta_0=0.02\)):

| \(D\) (Gpc) | A′ same radial speed | B \(V\propto r\) | C \(V\propto 1/r\) | D mild \(\kappa\) | kill \(a_2\) |
|--:|--:|--:|--:|--:|--:|
| \(14.7\) | \(0.71\) | \(0.57\) | \(0.60\) | \(0.69\) | \(\approx 0.039\) |
| \(22.28\) | \(0.55\) | \(0.60\) | \(0.68\) | \(0.56\) | \(\approx 0.427\) |
| \(49\) | \(0.63\) | \(0.67\) | \(0.73\) | \(0.63\) | \(\approx 0.384\) |

Amplitude scale factors to match kill peaks: \(\sim 8\)–\(28\). Shape correlation \(\lesssim 0.73\) is morphological interest, not identity.

Cases that spoil shared-\(\mathbf{V}\) flatness (peak \(|\Delta T/T|\) at \(\beta_0=0.02\)):

| Case | Peak \(|\Delta T/T|\) |
|--|--:|
| A′ same radial speed (not same vector) | \(\approx 0.035\) |
| B \(V\propto r\) | \(\approx 0.048\) |
| C \(V\propto 1/r\) | \(\approx 0.028\) |
| D mild shear \(\kappa=0.1\) | \(\approx 0.037\) |
| Plates, observer \(\beta\) only | \(\approx 0.020\) |

### 7.4 Multipole-partner cancellation (`simple_cmb_two_quad_cancel.json`)

Fit a low-order opposite partner to kill / CoE Doppler (monopoles removed). Target: \(10^{-5}\).

**At \(D=49\,\mathrm{Gpc}\):**

| Partner | Fraction of demeaned peak cancelled | Residual peak | Under \(10^{-5}\)? |
|--|--:|--:|--|
| Dipole + \(P_2\) | \(\approx 61\%\) | \(\approx 0.191\) | No |
| \(\ell\le 4\) | \(\approx 85\%\) | \(\approx 0.076\) | No |

**At \(D=14.7\) and \(22.28\):** \(\ell\le 4\) cancels \(\sim 22\%\)–\(39\%\); residuals still \(O(0.25\)–\(1.6)\).

Letter statement: exploratory two-emitter and multipole-partner constructions follow kill-function morphology at large \(D\) and remove most demeaned power through \(\ell\le 4\), but do **not** yet reach the CMB bound of order \(10^{-5}\).

---

## 8. Apparent acceleration (letter discussion)

Off-centre sampling illusion:

- No dark-energy component required.
- No actual outward acceleration of the material shell required.
- A contracting CoE geometry would share the **average** appearance of acceleration but differ in **angular** signature.

*Plain language:* What looks like late-time acceleration is a geometric sampling artefact at distance \(D\) from the CoE.

---

## 9. What we do next (open light-cone calculation)

1. Evolve \(D(t)\), \(R_{\mathrm{in}}(t)\), \(R_{\mathrm{out}}(t)\) from the present back through recombination (observer mid-shell).
2. Integrate null geodesics on that history to obtain \(T(\theta)\) on our sky (gravitational + kinematic together).
3. Ask whether shared bulk velocity of emitting surfaces plus the gravitational multipole can keep

\[
\ell \ge 2
\quad\text{below}\quad
\sim 10^{-5}
\]

while low-redshift \(A(\ell)\) remains Migkas-compatible.

4. If a residual remains, identify the smallest extra degree of freedom (recombination timing, shell-thickness history, or a mild partner multipole) that removes it without reintroducing a fixed-curve conflict across redshift slices.

---

## 10. Appendix — shell expansion asymptotics and pointing

LTB shell equation integrates to early radiation / late matter asymptotics (Fig. A1): flat matter+radiation track; dust \(E=0\); \(E>0\) expand-forever; \(E<0\) recollapsing.

CoE is a fixed celestial direction (fitted CoE coordinates). From Central Texas / Plano: rises NE, sets NW on a sidereal schedule. Figs. A2–A3 and Table A1: Sun–CoE timing through 2026 and nine Sun–CoE event crossings for 2027 (rise–rise, culmination–culmination, set–set diagonal).

---

## 11. Quick reference — all letter-quoted numbers

| Symbol / quantity | Value | Role |
|--|--|--|
| CMB isotropy | \(\sim 10^{-5}\) | External bound |
| \(R\) | \(\approx 14.3\,\mathrm{Gpc}\) | Vision / diagnostic radius |
| Shell rule | \(D\gtrsim 14.7\,\mathrm{Gpc}\) | Mid-shell bound |
| \(k_{\mathrm{in}}, k_{\mathrm{out}}\) | \(\approx 1.028,\; 1.475\) | Shell thickness rule |
| Rel. \(D, R_{\mathrm{in}}, R_{\mathrm{out}}\) | \(32.9,\; 18.1,\; 49.7\,\mathrm{Gpc}\) | Preferred relativistic headline |
| Newt. \(\Omega_s=0.049\) \(D,R_{\mathrm{in}},R_{\mathrm{out}}\) | \(49.17,\; 34.48,\; 70.27\,\mathrm{Gpc}\) | Solved Newtonian geometry |
| \(\Omega_s\) | \(0.049\) (or free \(1.01\)) | Shell density parameter |
| Free-\(\Omega_s\) \(D\) | \(93.7\,\mathrm{Gpc}\) | Table 1 col. 1 |
| \(\ell_0\) | \(12.8^{\circ}\)–\(13.5^{\circ}\) | CoE longitude phase |
| \(R^2\) vs Migkas | \(0.909 / 0.893 / 0.885\) | Fit quality (3 cols) |
| \(R^2\) vs \(A_F\) | \(0.987 / 0.977 / 0.963\) | Fit quality vs Fourier form |
| Tilt | \(+15^{\circ}\) | Sky orientation |
| Roll | \(+80.8^{\circ}\) | Sky orientation |
| Cluster \(a_2\) | \(0.084\) at \(R_s=0.31\,\mathrm{Gpc}\) | Core-mass anchor |
| \(\langle H\rangle\) | \(71.6\) (cycles) / \(\sim 71.25\) (newt solve) | Velocity unit |
| \(D_{\mathrm{null}}\) | \(22.27586\,\mathrm{Gpc}\) | Static gravitational null |
| Kill \(a_2\) (fixed \(R\)) | \(\sim 0.039\)–\(0.41\) | Required cancel quadrupole |
| Partner cancel at \(D=49\) | \(\sim 61\%\) (dip+\(P_2\)), \(\sim 85\%\) (\(\ell\le 4\)) | Exploratory, not \(10^{-5}\) |

---

## 12. Source map

| Letter topic | Primary files |
|--|--|
| Letter equations & Table 1 | `letter.tex`, attached CoE PDF, `_draft_from_chris.md` |
| Kill / Doppler / static | `kill_function_from_geometry.md`, `cycle4/dtt/`, `cycle4/lib4.py` |
| Shell rule / Newtonian solve | `letter_brief/figs/newt_ltbgeom_A_vs_ell.json`, `cycle4/README.md` |
| Partner cancel / two emitters | `cycle6_two_emitters/`, `simple_cmb_geometry.md` |
| Logic path | `logic_path_CMB_flatness.md` |
| Donut volumes | `cycle7_recombination/donut_storyboard_notes.md` |

