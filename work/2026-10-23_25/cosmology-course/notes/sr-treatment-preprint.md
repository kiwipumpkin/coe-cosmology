# Special-relativistic Doppler and central-mass correction to the centre-of-expansion shell geometry

## Motivation

A centre-of-expansion (CoE) cosmology places the observer off the preferred centre of a spherically expanding shell. The cosmic microwave background (CMB) last-scattering surface then sits at unequal distances from that centre, so a pure geometric shell model of the Migkas $A/A_{\rm all}$ anisotropy produces a fractional temperature (or amplitude) swing of order several percent. Observed CMB anisotropies after the kinematic dipole are of order $10^{-5}$. The question addressed here is whether special relativity (SR) applied to the radial motion of the emitter and the observer, together with gravitational redshift from a central mass at the preferred origin, can cancel most of that geometric anisotropy, and what residual pattern remains once that cancellation is carried through.

This section replaces an earlier “kill-function” description. It states the model, the fitted parameters, the residual after the SR-plus-mass correction, and a Fourier analysis of that residual about the shell axis. Numbers are taken from the existing fits in `cmb-cancel.json`, `cmb-cancel-story.json`, `cmb-cancel-minus-m1to12.json`, and `migkas-geometry-fit.json`; no new optimisation is reported here.

## Geometry input

The geometric target is the normalised Migkas shell curve from the weighted least-squares (WLS) fit in `migkas-geometry-fit.json`:

- visible-universe radius $D = 58.2\,\mathrm{Gpc}$,
- shell radius $R = 14.3\,\mathrm{Gpc}$,
- shell axis $\ell_0 = 15.16^\circ$ in Galactic longitude,
- weighted $R^2 = 0.55$–$0.58$, $\chi^2 = 18.75$ on 71 bins.

Let $A/A_{\rm all}(\ell)$ be that full-precision WLS curve on a $1^\circ$ grid $\ell = 0,\ldots,359$. The normalised geometry is

$$
G(\ell)=\frac{A/A_{\rm all}(\ell)}{\langle A/A_{\rm all}\rangle},
$$

where $\langle\cdot\rangle$ denotes the mean over the $360$ longitude samples. The root-mean-square (RMS) of $G-1$ is

$$
\mathrm{RMS}(G-1)=5.75\times10^{-2}.
$$

A multiplicative correction $S(\ell)$ that would cancel $G$ exactly, up to an overall scale, is

$$
S_{\rm req}=\frac{1/G}{\langle 1/G\rangle}.
$$

The model below constructs a physical $S_{\rm model}$ from SR Doppler and gravitational redshift and compares it to $S_{\rm req}$.

## Definitions and conventions

Units: $c=1$, and the observer’s radial offset from the preferred centre $O$ in the reference configuration is $d_1=1$. Angles are measured from the fitted shell axis,

$$
\theta=\ell-\phi,
$$

with $\phi$ free in the fit (it returns to $\ell_0$).

- Preferred centre $O$ at $(0,0)$; central mass sits at $O$.
- Observer in the reference configuration at $(d_1,0)=(1,0)$.
- Emitter (CMB last-scattering parcel in the direction $\theta$) at reference distance $q$ from the observer:

$$
\mathbf{x}_{e0}=(1+q\cos\theta,\ q\sin\theta),\qquad
r_{e0}=\sqrt{1+q^2+2q\cos\theta}.
$$

Thus $\theta=0$ points away from $O$ (anti-centre) and $\theta=180^\circ$ toward $O$.

- A common expansion scale $s$ maps reference radii to present radii: $r_o=s$, $r_e=s\,r_{e0}$.
- Velocities are radial from $O$, outward positive. Their magnitudes follow from a Hubble-like term plus the work against a central mass parameter $\mu$ (with dimensions of length in these units, i.e.\ $GM/c^2$ measured in units of $d_1$):

$$
v^2=\beta_1^2 r_0^2+2\mu\left(\frac{1}{r}-\frac{1}{r_0}\right),
$$

applied separately at the observer ($r_0=1$, $r=r_o=s$) and at the emitter ($r_0=r_{e0}$, $r=r_e=s\,r_{e0}$). The dimensionless parameter $\beta_1$ is therefore the free-expansion speed at unit radius in the reference configuration (the product of the expansion rate and $d_1/c$).

- Velocity vectors:

$$
\boldsymbol{\beta}_o=v_o\,\hat{\mathbf{x}},\qquad
\boldsymbol{\beta}_e=v_e\,\mathbf{x}_e/r_e.
$$

- Photon flight from emitter to observer: with $\mathbf{d}=\mathbf{x}_o-\mathbf{x}_e$ and the observer moving along $+\hat{\mathbf{x}}$, the positive root $t_f$ of

$$
(v_o^2-1)\,t_f^2+2\,d_x\,v_o\,t_f+|\mathbf{d}|^2=0
$$

gives the propagation direction

$$
\hat{\mathbf{n}}=\frac{\mathbf{d}+v_o t_f\,\hat{\mathbf{x}}}{t_f}.
$$

## Special-relativistic and gravitational factors

The relativistic Doppler factor between emitter and observer is

$$
D_{\rm motion}=\frac{\gamma_o\bigl(1-\boldsymbol{\beta}_o\cdot\hat{\mathbf{n}}\bigr)}{\gamma_e\bigl(1-\boldsymbol{\beta}_e\cdot\hat{\mathbf{n}}\bigr)},
\qquad
\gamma=\frac{1}{\sqrt{1-v^2}}.
$$

The gravitational redshift between the emitter and observer radii, for the same central mass $\mu$, is

$$
D_{\rm grav}=\sqrt{\frac{1-2\mu/r_e}{1-2\mu/r_o}}.
$$

The model correction is their product,

$$
S_{\rm model}=D_{\rm motion}\,D_{\rm grav}.
$$

The residual after applying this correction to the geometry and re-normalising is

$$
R(\ell)=\frac{G\,S_{\rm model}}{\langle G\,S_{\rm model}\rangle}-1.
$$

The objective minimised in the existing fits is $\mathrm{RMS}(R)$. Physical guards reject parameter points with $v^2\notin(0,1)$ or $2\mu/r\ge1$ anywhere on the sky.

## Fitted parameters (run with $s\le3$)

Three existing runs are recorded in `cmb-cancel.json`. Run 1 allows only free expansion ($v=kr$, $\mu=0$, $s=1$) and collapses to a degenerate limit $\beta_1\to1$, $q\to0$. Run 3 allows $s$ up to $50$ and drives $s$ to that bound; because the intercept geometry is scale-free, $s$ enters only through $1/s$ in the speed formula, and the $s\to\infty$ limit is finite. The run used below is run 2, with $1\le s\le3$:

| Parameter | Value | Meaning in the code |
|-----------|-------|---------------------|
| $\beta_1$ | $0.5866$ | Free-expansion speed at unit reference radius; enters $v^2=\beta_1^2 r_0^2+\cdots$ |
| $q$ | $0.2303$ | Emitter distance from the observer in units of $d_1$ |
| $\mu$ | $0.09332$ | Central-mass length $GM/c^2$ in units of $d_1$; enters both $D_{\rm grav}$ and the $2\mu(1/r-1/r_0)$ term in $v^2$ |
| $s$ | $3.000$ (upper bound) | Common expansion scale of all radii relative to the reference configuration |
| $\phi$ | $15.16^\circ$ | Axis of the Doppler/mass pattern; coincides with the shell axis $\ell_0$ |

At this point the present observer speed is $v_o=0.469\,c$; emitter speeds range from $0.205\,c$ to $0.648\,c$; the largest $2\mu/r$ on the sky is $0.081$.

## Results after SR plus central mass

| Quantity | Value |
|----------|-------|
| Starting geometry RMS, $\mathrm{RMS}(G-1)$ | $5.75\times10^{-2}$ |
| Residual RMS after SR + mass, $\mathrm{RMS}(R)$ | $2.37\times10^{-3}$ |
| Maximum $\|R\|$ | $3.76\times10^{-3}$ |
| Fraction of starting RMS removed (amplitude) | $95.9\%$ |
| Fraction of starting variance removed | $99.83\%$ |

Equivalently, the starting anisotropy is about $5{,}750$ times the CMB $10^{-5}$ level; after SR + mass it is about $237$ times that level.

## Solving the residual with Fourier analysis

The residual $R(\ell)$ is expanded in Fourier modes about the shell axis $\phi=15.16^\circ$:

$$
R(\ell)=a_0+\sum_{m=1}^{M}\Bigl[c_m\cos\bigl(m(\ell-\phi)\bigr)+s_m\sin\bigl(m(\ell-\phi)\bigr)\Bigr].
$$

Every sine coefficient satisfies $|s_m|\sim10^{-8}$–$10^{-12}$, so the residual is a pure cosine series about $\ell_0$, analogous to the zonal harmonics $J_n$ of the Earth’s gravitational field. Writing $c_m$ for the cosine amplitudes about $\phi$, the existing decomposition through $M=12$ gives leftover RMS

$$
\mathrm{RMS}\Bigl(R-\sum_{m=1}^{12}c_m\cos\bigl(m(\ell-\phi)\bigr)\Bigr)=1.83\times10^{-6}
$$

on the full-precision geometry curve. On the four-decimal Migkas CSV the same procedure leaves $2.79\times10^{-5}$, which matches the rounding difference between that CSV and the full-precision curve ($\mathrm{RMS}(G_{\rm round}-G_{\rm full})=2.85\times10^{-5}$).

This Fourier step is an analysis of the residual that remains after the physical SR-plus-mass correction. It does not replace that correction. The coefficients $c_m$ are the amplitudes that a complete physical model of the leftover pattern would have to reproduce.

### Origin of each mode

To first order in small fractional contrasts,

$$
R\approx\ln G-\ln(1/S_{\rm model})=\ln G-\ln(1/D_{\rm motion})-\ln(1/D_{\rm grav}),
$$

so each cosine coefficient of the residual is the difference of the corresponding coefficients of $\ln G$ and of $\ln(1/S_{\rm model})$. Measured amplitudes from `cmb-cancel-story.json`:

| $m$ | Period | $c_m$ | $\%$ of residual variance | $(\ln G)_m$ | $(\ln 1/S)_m$ | Origin |
|-----|--------|-------|---------------------------|-------------|---------------|--------|
| 1 | $360^\circ$ | $+1.80\times10^{-4}$ | $0.29\%$ | $-5.68\times10^{-3}$ | $-5.86\times10^{-3}$ | Nearly cancelled by SR + mass |
| 2 | $180^\circ$ | $-2.73\times10^{-4}$ | $0.66\%$ | $+7.81\times10^{-2}$ | $+7.83\times10^{-2}$ | Nearly cancelled by SR + mass |
| 3 | $120^\circ$ | $-2.00\times10^{-3}$ | $35.6\%$ | $-2.08\times10^{-2}$ | $-1.88\times10^{-2}$ | Shell spectrum vs slower-falling Doppler |
| 4 | $90^\circ$ | $-1.85\times10^{-3}$ | $30.3\%$ | $+2.42\times10^{-3}$ | $+4.27\times10^{-3}$ | Shell spectrum vs slower-falling Doppler |
| 5 | $72^\circ$ | $+1.73\times10^{-3}$ | $26.6\%$ | $+1.96\times10^{-4}$ | $-1.53\times10^{-3}$ | Shell spectrum vs slower-falling Doppler |
| 6 | $60^\circ$ | $-8.08\times10^{-4}$ | $5.8\%$ | $-1.13\times10^{-4}$ | $+6.97\times10^{-4}$ | Transition |
| $\ge7$ | $\le51.4^\circ$ | $\lesssim2.8\times10^{-4}$ | $\lesssim0.7\%$ each | $\sim0$ in $\ln G$ | geometric Doppler tail | Point-source Doppler tail |

Modes $m=1$ and $m=2$ of $\ln G$ are matched almost exactly by $\ln(1/S_{\rm model})$, which is why the SR-plus-mass correction removes $95.9\%$ of the starting RMS. Modes $m=3$, $4$, and $5$ together carry about $92\%$ of the leftover variance: the shell geometry’s multipole spectrum falls faster with $m$ than the relativistic Doppler pattern of a point-mass-plus-Hubble flow, so those three coefficients do not cancel. Modes $m\ge7$ are almost absent from $\ln G$ and are the geometric high-$m$ tail of the Doppler factor alone.

## Caveats

1. **Parameter $s$ sits on its upper bound.** In run 2, $s=3$ is the imposed ceiling. Run 3 with $s\le50$ reaches $s=50$ with a slightly lower RMS ($2.15\times10^{-3}$). The intercept geometry is scale-free, so $s$ is only weakly constrained; the physical content of the fit is carried by $\beta_1$, $q$, $\mu$, and $\phi$.
2. **The Fourier series is diagnostic, not a dynamical model.** Solving $R$ with cosines about $\ell_0$ identifies which multipoles a complete CoE-plus-SR theory must still explain. It is not a claim that those multipoles have been physically cancelled.
3. **Residual vs.\ CMB.** After SR + mass, $\mathrm{RMS}(R)=2.37\times10^{-3}$ is still $\sim237$ times the CMB $10^{-5}$ level. Reaching the CMB amplitude with the Fourier content of $R$ requires accounting for the $m=3$–$5$ mismatch (and smaller higher modes) by additional physics, not by the two-parameter SR-plus-mass model alone.
4. **Friedmann–Lemaître–Robertson–Walker (FLRW) comparison.** This calculation does not assert that the Cosmological Principle (CP) fails, nor that FLRW is ruled out. It quantifies how much of one specific geometric anisotropy is cancelled by SR Doppler and a central mass under CoE kinematics.
5. **Exploratory four-clump sketch.** An optional tidal model with four equal masses about the same axis (existing `clump-fit-residual.json`) reduces the post-SR residual RMS from $2.37\times10^{-3}$ to $8.34\times10^{-4}$ (about $88\%$ of the post-SR variance). That sketch is not part of the SR treatment proper and is recorded only as an existence check that higher multipoles can be sourced by off-axis mass.

## Figures supplied by the author

These three figures go with this section, ahead of the existing figures listed below.

1. **Photon intercept geometry** (`notes/sr-fig-photon-intercept.png`). Exact photon paths between two objects that move radially away from the centre $O$ at frozen speeds, with $c=1$. Object 1 starts at $(1,0)$ and object 2 at about $(1,0.91)$. A photon sent from 1 reaches 2 at $t=1.169$, and the return photon reaches 1 at $t=2.337$. The figure shows the light-travel-time root $t_f$ and the arrival direction $\hat{\mathbf n}$ that enter $D_{\rm motion}$.
2. **Doppler temperature pattern versus speed** (`notes/sr-fig-temperature-anisotropy.png`). The fractional blackbody-temperature anisotropy $(T_{\rm obs}-\bar T)/\bar T$, in percent, against emitter angle $\theta$ ($0^\circ$ outward, $180^\circ$ inward), for $v(d_1)=0.01c$, $0.203c$ and $0.278c$. At $0.01c$ the pattern is flat. At $0.203c$ it runs from about $-0.45\%$ outward to $+0.26\%$ near $110^\circ$ and $250^\circ$. At $0.278c$ it runs from about $-1.2\%$ to $+0.73\%$. In every case there is a shallow local minimum near $180^\circ$, which shows how the SR sky pattern grows with expansion speed.
3. **Exploratory four-clump tidal sketch** (`notes/sr-fig-4clumps.jpg`). The post-SR residual compared with the tidal pull of four equal masses of $7.7\times10^{51}$ kg at $24.5$ Gpc, at $\ell=75^\circ$, $162^\circ$, $228^\circ$ and $315^\circ$. The masses explain 88% of the post-SR variance. This figure belongs with caveat 5 and is not part of the SR treatment itself.

## Suggested figures (existing paths)

| Figure | Path | Content |
|--------|------|---------|
| Equations summary | `notes/cmb-cancel-eqs.png` | Full equation list and run-2 / run-3 numbers |
| Preprint equation card | `notes/sr-treatment-eqs.png` | Key equations wrapped for the preprint |
| Curves and residual | `notes/cmb-cancel-curves.png` | $G$, $S_{\rm req}$, $S_{\rm model}$, and $R(\ell)$ |
| Mode amplitudes | `notes/cmb-cancel-fourier.png` | $\|A_m\|$ of $R$ for runs 1–3 |
| Step-by-step story | `notes/cmb-cancel-story.png` | Geometry → SR+mass → residual → mode sum → leftover |
| Mode origin | `notes/cmb-cancel-mode-origin.png` | $\|(\ln G)_m\|$, $\|(\ln 1/S)_m\|$, $\|c_m\|$ |
| Polar mode sheet | `notes/cmb-modes-sheet.png` | Each $m=1$–$12$ as a polar plot plus table |
| Waterfall | `notes/cmb-waterfall.png` | RMS in CMB units after each analysis step |
| Geometry fit | `notes/migkas-geometry-fit.png` | Migkas bins and WLS shell curve |

## Source files

- `notes/cmb-cancel.py`, `notes/cmb-cancel.json` — model, fits, Fourier of $R$
- `notes/cmb-cancel-story.py`, `notes/cmb-cancel-story.json` — mode-origin coefficients
- `notes/cmb-cancel-minus-m1to12.json` — leftover after solving $m=1$–$12$
- `notes/migkas-geometry-fit.json` — $D$, $R$, $\ell_0$, $\chi^2$, $R^2$
- `notes/clump-fit-residual.json` — optional four-clump sketch (exploratory only)
