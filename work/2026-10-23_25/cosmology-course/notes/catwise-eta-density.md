# CatWISE-constrained η-density (generalized exponential) + Migkas refit

**Date:** 2026-10-06 (America/Chicago)

**Acronyms:** CatWISE = CatWISE quasar dipole catalogue (Secrest et al. 2021); CoE = centre of expansion; Migkas = Migkas et al. cluster scaling $A/A_{\mathrm{all}}$ anisotropy.

## Density law

$$\rho(r)=\rho_0\,e^{-(r/L)^\eta},\qquad \rho_0=\rho(D)\,e^{(D/L)^\eta}$$

- $\eta=1$ recovers the pure exponential $\rho(r)=\rho_0 e^{-r/L}$.
- Finite at CoE: $\rho(0)=\rho_0$; $\rho\to 0$ as $r\to\infty$.
- Earth–CoE baseline for density plot: $D=58.2$ Gpc.
- Local normalisation: $\rho(D)=\Omega_m\rho_{\rm crit}$ with **$\Omega_m=0.3$**.
- $\rho_{\rm crit}(H_0=70)\approx 1.3599e+11\,M_\odot/\mathrm{Mpc}^3$.
- Local slope: $(\rho'/\rho)_D=-(\eta/L)(D/L)^{\eta-1}$.
- BH-density flag: $\rho(0)/\rho_{\rm crit}=1.0319$ (below stellar-BH reference; no floor invented).

## CatWISE 3% lock (only fixed constraint)

Pole-to-pole contrast at vision distance $d$:

$$\delta=\frac{\rho(D-d)-\rho(D+d)}{\rho(D)}=e^{-(D/L)^\eta\big((1-d/D)^\eta-1\big)}-e^{-(D/L)^\eta\big((1+d/D)^\eta-1\big)}\approx 0.031.$$

For $\eta=1$ this collapses to $\delta=2\sinh(d/L)$, so $L=d/\mathrm{arcsinh}(\delta/2)$ (independent of $D$).

Primary $d=R_{\rm vis}\approx14.3$ Gpc; secondary $d=3$ Gpc.

| case | $\eta$ | $d$ (Gpc) | $D$ (Gpc) | $L$ (Gpc) | $\delta$ check | $\rho(0)/\rho_{{\rm crit}}$ | $(\rho'/\rho)_D$ (Gpc$^{-1}$) |
|---|---:|---:|---:|---:|---:|---:|---:|
| eta=0.0500009_d=14.3005 | 0.05 | 14.3 | 58.2 | 0.849 | 0.03100 | 1.0319 | -0.00106135 |
| eta=1_d=14.3005 | 1 | 14.3 | 58.2 | 922.6 | 0.03100 | 0.31953 | -0.00108384 |

**η=1 recovery check:** $L=922.6$ Gpc vs $d/\mathrm{arcsinh}(\delta/2)=922.6$ Gpc (relative difference 0.0000%; success criterion ~1%).

## Migkas $A/A_{\rm all}$ refit (CatWISE $\delta$ lock)

Geometry only: $A/A_{\rm all}(\ell)=a_f(\ell-\ell_0)/\langle a_f\rangle$, Ob=0, same WLS as migkas-geometry-fit / CatWISE-exp. Free $\{D, m_c, \ell_0, \eta\}$; $L$ locked by $\delta=0.031$ (not freed independently).

| case | $d$ | $\eta$ | $L$ (locked) | $D$ | $L/D$ | $m_c$ | $\ell_0$ | $R^2$ unw. | $\chi^2$ | $\sigma$ vs null | $A_{\rm amp}$ |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| d14.3 | 14.3 | 0.05 | 2.905e+04 | 35 | 831 | 4.344 | 19.91 | 0.7453 | 9.278 | 4.99 | 0.1971 |
| d3 | 3 | 0.08266 | 0.003435 | 34.4 | 0.0001 | 2.854 | 19.77 | 0.7480 | 9.174 | 4.99 | 0.1971 |

### Comparison table (σ vs null $\chi^2_{\rm null}=45.137$, $N=71$)

| model | free | $R^2$ | $\chi^2$ | $\sigma$ nested |
|---|---|---:|---:|---:|
| ptshell WLS | 5 | 0.546 | 18.7 | 3.79 |
| CatWISE-$L$ exp ($d\approx 14.3$) | 3 | 0.716 | 10.6 | 5.12 |
| free-exp {$D,L,m_c,\ell_0$} | 4 | 0.759 | 8.5 | 5.06 |
| **CatWISE-$\eta$ (this work, $d\approx 14.3$)** | 4 | **0.745** | **9.3** | **4.99** |
| CatWISE-$\eta$ ($d=3$) | 4 | 0.748 | 9.2 | 4.99 |

**Best-fit (primary $d\approx14.3$ Gpc):** $\eta=0.05$, $D=35$ Gpc, $L=2.905e+04$ Gpc, $m_c=4.344$, $\ell_0=19.91^\circ$; $R^2=0.7453$, $\chi^2=9.278$, $\sigma_{\rm null}=4.99$.

**η boundary:** Under the CatWISE $\delta=0.031$ lock, $\chi^2$ improves monotonically as $\eta$ decreases ($\eta=1\to\chi^2\approx10.57$; $\eta=0.05\to\chi^2\approx9.28$). $\eta\lesssim 0.02$ fails the $L$-solve / dynamics numerically, so the reported best $\eta\approx 0.05$ is the lower viable edge, not an interior minimum.

**$d=3$ caveat:** locked $L=0.003435$ Gpc ($L/D=0.0001$) sits at a tiny-$L$ edge; treat as exploratory.

## Files

| Path | Content |
|---|---|
| `catwise-eta-density.png` | $\rho(r)$ with best-fit $\eta$ (+ $\eta=1$ overlay) |
| `catwise-eta-density.csv` | radial profile columns |
| `catwise-eta-density.md` | this note |
| `catwise-eta-migkas-fit-d14.3.png` / `.csv` | Migkas curve, primary $d\approx 14.3$ |
| `catwise-eta-migkas-fit-d3.png` / `.csv` | Migkas curve, $d=3$ (if run) |
| `catwise-eta-migkas-fit.json` | machine-readable fits + density params |
| `catwise-eta-density-migkas.py` | script |
| `prof_core.py` (bubble) | minimal patch: `kind=="eta"` → `np.exp(-(x**a["eta"]))`; LB/LOGV `eta` |
