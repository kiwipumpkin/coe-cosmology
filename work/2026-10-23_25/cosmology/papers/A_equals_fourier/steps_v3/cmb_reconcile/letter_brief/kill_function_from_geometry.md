# Kill function from CoE geometry — how it was calculated, and how plate / shell Dopplers relate

**Status:** exploratory stress-test note for Chris's Centre-of-Expansion (CoE) / donut-universe geometry.  
**Not a claim** that CMB isotropy at \(\sim 10^{-5}\) has been achieved. Partner cancellation removes a large fraction of demeaned power at large \(D\), but residuals remain \(\gg 10^{-5}\).

Primary sources on disk:

- Logic path: `../logic_path_CMB_flatness.md`
- Fixed-\(R\) kill pipeline: `../cycle4/dtt/` (`compute_dtt.py`, `dtt.json`, `plot_kill_term*.py`)
- Shared helpers: `../cycle4/lib4.py`
- Recombination-track variant: `../cycle4/dtt_rec_track/`
- Two-emitter / plate–shell toys: `../cycle6_two_emitters/`
- Partner-cancel summary: `../cycle6_two_emitters/simple_cmb_geometry.md`

---

## 1. Intro

Chris's CoE model places us (US) **off-centre** inside a thick ordinary-matter shell, at distance \(D\) from the Centre of Expansion. Local Hubble anisotropy is interpreted as a free-fall velocity difference relative to the CoE, fitted to Migkas directional \(H\) data. That fit **fixes** CoE direction and shell parameters; the CMB is then an **external** stress test, not a fit target.

The same geometry that fits Migkas predicts a gravitational multipole / Doppler residual on the last-scattering sky. Left alone, that residual is far above the observed CMB anisotropy scale \(\sim 10^{-5}\). The **required kill function** is defined as the exact opposite of that residual along a great circle about the CoE axis — the fractional \(\Delta T/T\) (or \({\rm d}E/E\)) pattern that would cancel it point-by-point after monopole removal.

Cycle 6 then asks a related geometry question: if US sits between **two infinite plates** or **two concentric shells** that share a bulk velocity with the observer, when is the sky flat, and when do plate/shell Doppler residuals morphologically resemble the kill function?

---

## 2. Geometry \(\rightarrow\) Doppler residual

### Working geometry (headline fixed-\(R\) stress test)

- Vision / emission radius \(R \approx 14.3\,\mathrm{Gpc}\) (observer-centred sphere used for today's-field diagnostics).
- Thick shell: \(\Omega_s = 0.049\),
  - \(R_{\rm in} = D - 1.028\,R\),
  - \(R_{\rm out} = D + 1.475\,R\).
- Compact core (point) with mass \(M_c\) re-solved at each \(D\) so the cluster quadrupole stays at \(a_2 = 0.084\) at \(R_s = 0.31\,\mathrm{Gpc}\).
- Mean conversion \(\langle H \rangle = 71.6\,\mathrm{km\,s^{-1}\,Mpc^{-1}}\).
- Preferred fixed-\(R\) distances:
  - \(D = 14.7\,\mathrm{Gpc}\) (shell-minimum side),
  - \(D = 22.28\,\mathrm{Gpc}\) (static-null; \(D_{\rm null} = 22.27586\) from `../cycle4/partB.json`),
  - \(D = 49\,\mathrm{Gpc}\) (headline-like).
- Angle convention on plots: \(\theta = 0^\circ\) toward CoE, \(180^\circ\) anti-CoE.

### What enters \(\Delta T/T\)

On the emission sphere, cycle 4/lib4 split the Newtonian today's-field residual into:

1. **Kinematic Doppler** (linear convention, same as cycle 1)

\[
\left(\frac{\Delta T}{T}\right)_{\rm Doppler}
=
-\,\hat{\mathbf n}\cdot\frac{\mathbf V_{\rm emit}-\mathbf V_{\rm obs}}{c}.
\]

2. **Static gravitational**

\[
\left(\frac{\Delta T}{T}\right)_{\rm static}
=
\frac{\Phi_{\rm emit}-\Phi_{\rm obs}}{c^2}.
\]

By the shell theorem, the static pattern on the \(R\)-sphere depends only on
\(M_{\rm eff} = M_{\rm interior} - M_{\rm cavity}\). Static can be nulled at \(D \approx 22.28\,\mathrm{Gpc}\) (core mass = cavity mass). **Doppler cannot** be removed by any interior profile that keeps the Migkas cluster pattern (cycle 4 no-go: \(V/r\) constant over \(D-R\ldots D+R\) would kill shear at Earth).

### Stress-test numbers (fixed \(R=14.3\,\mathrm{Gpc}\), from `dtt.json`)

- \(D=14.7\): Doppler \(a_2 \approx -0.0387\); peak-to-peak \(\approx 0.616\); \(M_c/M_{\rm crit}\approx 0.00187\). Note: \(v_{r,{\rm min}}<0\) (turnaround near axis).
- \(D=22.28\): Doppler \(a_2 \approx -0.411\); peak-to-peak \(\approx 3.03\); static \(a_2 \approx 0\).
- \(D=49\): Doppler \(a_2 \approx -0.385\); peak-to-peak \(\approx 0.795\); static \(a_2 \approx +0.087\).

Target CMB scale for \(\ell\ge 2\): \(|a_\ell|\lesssim 10^{-5}\). Fixed-\(R\) Doppler \(a_2\) is \(\sim 0.04\)–\(0.4\) — thousands to tens of thousands times too large. Caveat: \(HR/c \approx 3.42\), so linear \(v/c\) is outside its formal range; a fractional estimate \(\mathrm{Doppler}/(HR/c)\) still leaves \(|a_2|\sim 0.01\)–\(0.12 \gg 10^{-5}\).

---

## 3. Definition of the kill function

After mean (monopole) subtraction on the tabulated Doppler curve,

\[
\mathrm{kill}(\theta)
=
-\Bigl[
\bigl(\tfrac{\Delta T}{T}\bigr)_{\rm Doppler}(\theta)
-
\bigl\langle\tfrac{\Delta T}{T}\bigr\rangle_{\rm Doppler}
\Bigr].
\]

Equivalently on the full axisymmetric great circle \(\phi\in[0^\circ,360^\circ]\) about the CoE axis (mirror of the \(\theta=0\)–\(180^\circ\) half):

\[
\mathrm{kill}(\phi)
=
-\,\Delta T/T_{\rm Doppler}\big|_{\rm demeaned},\qquad
\mathrm{kill}(\phi)=\mathrm{kill}(360^\circ-\phi).
\]

- Required kill quadrupole is then \(a_2^{\rm kill} = -a_2^{\rm Doppler}\):
  - \(D=14.7\): \(a_2^{\rm kill} \approx +0.0387\)
  - \(D=22.28\): \(a_2^{\rm kill} \approx +0.411\)
  - \(D=49\): \(a_2^{\rm kill} \approx +0.385\)
- Photon-energy form used in plots: \(100\times\mathrm{kill}\) as percent \({\rm d}E/E\) (positive = energy gain). Peak absolute required corrections on the fixed-\(R\) circle:
  - \(D=14.7\): \(\approx 39\%\)
  - \(D=22.28\): \(\approx 263\%\) (inflated by turnaround spike; quadrupole is the robust metric)
  - \(D=49\): \(\approx 49\%\)

A successful cancellation would leave a post-cancel residual \(\approx 0\) point-by-point at the \(\sim 10^{-5}\) level. **That has not been demonstrated**; see §8–§10.

---

## 4. How the kill was computed (steps)

Pipeline in `../cycle4/dtt/`:

1. **Build mass model** (`lib4.build` / `lib4.run`): compact point core + \(\Omega_s=0.049\) shell; Chris Newtonian integrator (`prof_core.integrate`).
2. **Re-solve \(M_c\)** at each \(D\) so cluster \(a_2=0.084\) at \(R_s=0.31\,\mathrm{Gpc}\) (`compute_dtt.mcc_for`).
3. **Sample the \(R=14.3\,\mathrm{Gpc}\) sphere**: positions \(\mathbf X\) and velocities \(\mathbf V\) for observer + emitters; form look directions \(\hat{\mathbf n}\).
4. **Doppler** \(\Delta T/T = -\hat{\mathbf n}\cdot(\mathbf V_{\rm emit}-\mathbf V_{\rm obs})/c\) with velocity unit from \(\langle H\rangle\) conversion; **static** from \(M_{\rm eff}\).
5. **Legendre** \(a_\ell\) about the CoE axis (`lib4.legendre`); reverse code \(\theta\) so plots have \(\theta=0\) toward CoE.
6. **Kill** \(= -\)(Doppler − mean); extend \(0^\circ\)–\(180^\circ\) to \(0^\circ\)–\(360^\circ\) by axial reflection.
7. Write `dtt.json`; plot kill and energy-percent curves.

Fixed-\(R\) vs recombination track:

- Fixed \(R\approx 14.3\,\mathrm{Gpc}\) is the headline stress test used for kill \(a_2\sim 0.04\)–\(0.4\).
- `../cycle4/dtt_rec_track/` instead uses \(R_{\rm see}\approx c\,\Delta t \approx 4.25\,\mathrm{Gpc}\) with \(f=T_0/3000\approx 9.08\times 10^{-4}\). Kill \(a_2\) plateaus near \(\sim 0.11\)–\(0.115\) for \(D_{\rm now}\gtrsim 8\,\mathrm{Gpc}\) — still \(\sim 10^{4}\times\) over \(10^{-5}\). Smaller \(R_{\rm see}\) softens amplitude vs fixed \(R\) but does **not** pass the CMB test while keeping the cluster fit.

---

## 5. Plots / files produced

### Cycle 4 fixed-\(R\) kill

| File | Role |
|------|------|
| `../cycle4/dtt/compute_dtt.py` | Integrator → Doppler/static curves |
| `../cycle4/dtt/dtt.json` | Tabulated \(\theta\), \(\Delta T/T\), Legendre coeffs |
| `../cycle4/dtt/dtt_vs_theta.png` | Doppler / static / total vs \(\theta\) |
| `../cycle4/dtt/kill_term_vs_theta.png` | Kill vs \(\theta\) (\(0\)–\(180^\circ\)) |
| `../cycle4/dtt/kill_term_0_360.png` | Kill + Doppler residual on full great circle |
| `../cycle4/dtt/kill_energy_pct_0_360.png` | \(100\times\mathrm{kill}\) with \(\pm 10^{-5}\) band |
| `../cycle4/dtt/plot_kill_term_0_360.py` | Script that defines `kill = -demean(doppler)` |

### Cycle 4 recombination track

| File | Role |
|------|------|
| `../cycle4/dtt_rec_track/track.json` | Kill \(a_2\) along \(D_{\rm now}\) track |
| `../cycle4/dtt_rec_track/kill_term_0_360_track.png` | Kill circle at selected track \(D\) |
| `../cycle4/dtt_rec_track/kill_a2_vs_Dnow_Rsee.png` | Kill \(a_2\) vs \(D_{\rm now}\), \(R_{\rm see}\) |

### Cycle 6 two-emitters and comparison to kill

| File | Role |
|------|------|
| `../cycle6_two_emitters/compute_two_emitters.py` | Plate / shell \(T_{\rm obs}\) + Doppler toys |
| `../cycle6_two_emitters/two_emitters.json` | Rest-frame and Doppler skies |
| `../cycle6_two_emitters/plates_T_vs_theta.png`, `shells_T_vs_theta.png` | Great-circle \(T_{\rm obs}\) |
| `../cycle6_two_emitters/comoving_sameV_explore.json` | Shared-\(\mathbf V\) vs radial/shear cases |
| `../cycle6_two_emitters/vs_kill_term_0_360.json` | Shape correlations of cycle-6 residuals vs kill |
| `../cycle6_two_emitters/vs_kill_term_0_360.png`, `kill_vs_plates_shells.png` | Overlay / morphology plots |
| `../cycle6_two_emitters/simple_cmb_two_quad_cancel.json` | Dipole+\(P_2\) / \(\ell\le 4\) partner cancel |
| `../cycle6_two_emitters/simple_cmb_geometry.md` | Human-readable cancel summary |
| `figs/kill_vs_cycle6_shapes.jpg`, `figs/cancel_partner_residual.png` | Letter-brief figures |

---

## 6. Infinite-plate Doppler

Idealization in `compute_two_emitters.py`:

- Infinite parallel plates at \(z=\pm d\); observer at \(z\in(-d,d)\).
- Every upward ray hits the \(+\)plate; every downward ray hits the \(-\)plate (hemispheres fixed — **offset never changes solid angles**).
- Rest-frame blackbody: if both plates share \(T\), sky is **exactly isotropic** (mid or offset).

With observer boost \(\boldsymbol\beta = \beta\,\hat{\mathbf z}\) through an otherwise isotropic field, exact SR Doppler:

\[
T_{\rm obs}(\hat{\mathbf n})
=
\frac{T_{\rm rest}}{\gamma\,(1-\boldsymbol\beta\cdot\hat{\mathbf n})},
\qquad
\gamma=(1-\beta^2)^{-1/2}.
\]

For equal plate temperatures this is the classic **smooth dipole**. Toy amplitude used: \(\beta_{\rm obs}=0.01\) → demeaned peak \(|\Delta T/T|\approx 0.020\) (`comoving_sameV_explore.json`, “plates: observer \(\beta\) only”).

**Relation to kill:** plate observer-boost Doppler is a pure dipole. The CoE kill function from cycle 4 has a strong quadrupole and higher multipoles (radial expansion / shear geometry). A plate dipole alone is **not** the kill morphology; it is the simplest controlled Doppler on an otherwise flat two-emitter sky.

---

## 7. Concentric-shell Doppler

Idealization:

- Concentric shells \(R_{\rm in}\), \(R_{\rm out}\); observer at \(D\) with \(R_{\rm in}<D<R_{\rm out}\).
- Toy ratios in cycle 6: \(R_{\rm in}=1\), \(R_{\rm out}=3\), \(D=1.8\) → inner-disk half-angle \(\alpha=\arcsin(R_{\rm in}/D)\approx 33.7^\circ\).
- Rays with \(|\theta|<\alpha\) (toward CoE) hit the inner shell; others hit the outer shell.
- Rest-frame equal \(T\): sky again **exactly isotropic** — geometry only chooses *which* shell is hit.

Optional Hubble-like recession on emitters (first-order, same linear convention as cycle 4):

\[
\boldsymbol\beta_{\rm em}
=
\beta_{\rm rec}\,\frac{r_{\rm hit}}{R_{\rm out}}\,\hat{\mathbf e}_r,
\qquad
\left(\frac{\Delta T}{T}\right)
\approx
-\,\hat{\mathbf n}\cdot(\boldsymbol\beta_{\rm em}-\boldsymbol\beta_{\rm obs}).
\]

Comoving observer at radius \(D\) sets \(\beta_{\rm obs,radial}=\beta_{\rm rec}(D/R_{\rm out})\). Thin-shell models also jump at the limb \(\alpha\) where hit radius (hence \(\beta\)) is discontinuous.

**Relation to kill:** radial velocity fields on curved surfaces generate smooth angle dependence plus possible limb features — closer in *kind* to the CoE Doppler residual than a pure plate dipole, but still a toy (discontinuous \(r\), no GR, no last-scattering hypersurface).

---

## 8. Comparison to the kill (what matched / what didn't)

Unit-normalized shapes of cycle-6 residuals were compared to cycle-4 kill \(=-\)Doppler (`vs_kill_term_0_360.json`, toy \(\beta_0=0.02\)):

Cases (demeaned residuals left when emitters are **not** sharing one Cartesian bulk \(\mathbf V\) with US):

- **A′** same radial *speed* \(\beta_0\,\hat{\mathbf e}_r\) (directions differ) — not a shared vector
- **B** \(V\propto r\) homologous, US matched at \(D\)
- **C** \(V\propto 1/r\), US matched at \(D\)
- **D** mild radial gradient / shear \(\kappa\approx 0.1\) about \(D\)

Shape correlations with kill (Pearson on demeaned curves), fixed-\(R\) kill at three \(D\):

| \(D\) (Gpc) | A′ corr | B \(V\propto r\) | C \(V\propto 1/r\) | D mild \(\kappa\) | kill \(a_2\) |
|-------------|---------|------------------|--------------------|-------------------|--------------|
| 14.7 | 0.71 | 0.57 | 0.60 | 0.69 | \(\approx 0.039\) |
| 22.28 | 0.55 | 0.60 | 0.68 | 0.56 | \(\approx 0.427\) |
| 49 | 0.63 | 0.67 | 0.73 | 0.63 | \(\approx 0.384\) |

- **Matched (qualitative):** morphology of radial / shear residuals tracks kill better at larger \(D\); letter fig `figs/kill_vs_cycle6_shapes.jpg` shows unit-normalized curves converging in shape.
- **Not matched (quantitative):** cycle-6 toy peak amplitudes at \(\beta_0=0.02\) are \(\sim 0.03\)–\(0.05\), while kill peaks are \(\sim 0.39\)–\(2.6\). Amplitude scale factors to match kill peaks are \(\sim 8\)–\(28\). Shape correlation \(\lesssim 0.73\) is interesting, not identity. None of these residuals are at \(10^{-5}\).

### Partner cancel (exploratory)

Fit a low-order opposite partner to kill / CoE Doppler (`simple_cmb_two_quad_cancel.json`):

- Dipole-only, dipole+\(P_2\), or \(\ell\le 4\) Legendre partner; monopoles removed.
- At \(D=49\,\mathrm{Gpc}\):
  - dip+\(P_2\): \(\approx 61\%\) of demeaned peak cancelled; residual peak \(\approx 0.19\)
  - \(\ell\le 4\): \(\approx 85\%\) cancelled; residual peak \(\approx 0.076\)
- At \(D=14.7\) and \(22.28\): cancellations \(\sim 22\%\)–\(39\%\) for \(\ell\le 4\); residuals still \(O(0.25\)–\(1.6)\).
- **Under \(10^{-5}\)?** False for all listed \(D\) and partners.

Positive takeaway (as in `simple_cmb_geometry.md`): opposite partners *do* cancel; reaching \(10^{-5}\) needs a partner that tracks kill past \(P_2\), or a recombination geometry where the kill itself is already tiny (\(\varepsilon_{\rm rec}=D_{\rm rec}/R_{\rm see}\) small helps in principle — see logic path §§11, 15).

---

## 9. Shared-\(\mathbf V\) flat case

The pedagogically clean limit (`comoving_sameV_explore.json`):

- **Plates + no relative \(\beta\)**: flat (\(T_{\rm obs}=T\) everywhere).
- **Shells + US share one Cartesian bulk velocity vector** with both emitters: exactly flat — relative Doppler vanishes in a common rest frame.
- **Donut / letter picture:** US between two emitters that coast with the same \(\mathbf V\) ⇒ flat sky *before* gravity's multipole is added (`letter.tex` / `figs/two_emitters_comoving.png`).

What **spoils** flatness (same JSON):

- Same *radial speed* but not the same vector (A′): peak \(|\Delta T/T|\approx 0.035\) \(\sim 3500\times 10^{-5}\)
- \(V\propto r\) homologous (B): \(\approx 0.048\) — shear / limb jump
- \(V\propto 1/r\) (C): \(\approx 0.028\)
- Mild shear \(\kappa=0.1\) (D): \(\approx 0.037\)
- Plates with observer \(\beta\) only: classic dipole \(\approx 0.020\)

So: **shared bulk \(\mathbf V\)** is the isotropic two-emitter baseline; **radial expansion / differential \(V\)** reintroduces a kill-like anisotropy — the same conceptual distinction drawn in `logic_path_CMB_flatness.md` §§6–7 and §§12–13 (rigid translation vs radial expansion; parallel-plate limit vs curvature).

---

## 10. Open limits (honest)

- Fixed-\(R\) and recombination-track kills both remain \(\gg 10^{-5}\) while the Migkas cluster quadrupole is held fixed.
- Static gravitational term can be nulled (\(D\approx 22.28\,\mathrm{Gpc}\) or compensated profiles); **Doppler from anisotropic emitter motion cannot** by interior profile alone (cycle 4 verdict).
- Cycle-6 plate/shell models are SR / thin-surface toys: no full metric, no past light cone to recombination, no SW/ISW, discontinuous thin-shell \(\beta(r)\).
- Partner cancel at \(\ell\le 4\) removes \(\sim 85\%\) at \(D=49\,\mathrm{Gpc}\) — encouraging morphology, **not** a CMB pass.
- Linear Doppler with \(HR/c\sim O(1)\)–few is formally outside \(v\ll c\); fractional rescaling does not change the order-of-magnitude failure.
- The decisive open calculation remains the one in `logic_path_CMB_flatness.md` §16: null geodesics on the evolving metric to the recombination hypersurface, with \(D(t)\), \(R_{\rm in}(t)\), \(R_{\rm out}(t)\), and \(\varepsilon_{\rm rec}\).

**Bottom line:** the kill function is a well-defined, reproducible output of Chris's CoE geometry under the fixed-\(R\) (and track) stress tests — namely \(-\)Doppler after mean subtraction on the CoE great circle. Infinite plates and concentric shells clarify *when* two-emitter skies are flat (shared \(\mathbf V\), equal rest-frame \(T\)) and show that radial/shear Dopplers can resemble the kill in shape at large \(D\). They do not yet cancel it to \(10^{-5}\).

