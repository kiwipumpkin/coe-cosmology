# CatWISE-constrained exponential CoE density + Migkas refit

**Date:** 2026-10-06 (America/Chicago)

## Density law

$$\rho(r)=\rho(D)\,e^{-(r-D)/L}=\rho_0\,e^{-r/L},\qquad \rho_0=\rho(D)\,e^{D/L}$$

- Finite at $r=0$; smooth; $\rho\to 0$ as $r\to\infty$ (marked on plot).
- Earth–CoE for the density plot: $D=58.2$ Gpc (ptshell / vision consistency; not refit here).
- Local normalisation: $\rho(D)=\Omega_m\rho_{\rm crit}$ with **$\Omega_m=0.3$** (stated choice; not $\Omega_s=0.049$).
- $\rho_{\rm crit}(H_0=70)\approx 1.3599\times 10^{11}\,M_\odot/\mathrm{Mpc}^3$.
- $\rho(0)\ll\rho_{\rm BH}$ for both $d$ cases below.

## Scale length from CatWISE contrast

Pole-to-pole contrast $\delta=[\rho(D-d)-\rho(D+d)]/\rho(D)=2\sinh(d/L)$.

Target $\delta=0.031\approx 2\mathcal{D}_{\rm CatWISE}$ with $\mathcal{D}=0.01554$ (Secrest et al. 2021). Hemisphere-scale would use $\delta\approx\mathcal{D}=0.01554$ and roughly double $L$.

$$L=\frac{d}{\mathrm{arcsinh}(\delta/2)}\approx\frac{2d}{\delta}\quad(d\ll L)$$

Derivative at Earth (pure exp): $(\rho'/\rho)_D=-1/L$, so $(d\rho/dr)_D=-\rho(D)/L$.

| $d$ (Gpc) | role | $L$ (Gpc) | $L\approx 2d/\delta$ | $\rho(0)/\rho_{\rm crit}$ | $(d\rho/dr)_D/\rho_{\rm crit}$ (Gpc$^{-1}$) | $(\rho'/\rho)_D$ (Gpc$^{-1}$) |
|---:|---|---:|---:|---:|---:|---:|
| 14.30 | **PRIMARY** (vision $R$) | 922.6 | 922.6 | 0.31953 | $-3.252\times 10^{-4}$ | $-1.084\times 10^{-3}$ |
| 3.0 | secondary (quasar-depth) | 193.6 | 193.5 | 0.40524 | $-1.550\times 10^{-3}$ | $-5.166\times 10^{-3}$ |

At Earth: $\rho(D)/\rho_{\rm crit}=0.3$. Plot annotates the Earth tangent for the primary curve.

## Migkas $A/A_{\rm all}$ refit ($L$ fixed by CatWISE)

Geometry only: $A/A_{\rm all}(\ell)=a_f(\ell-\ell_0)/\langle a_f\rangle$, pure exp, $O_b=0$. Free $\{D,\,m_c,\,\ell_0\}$; $L$ locked. No Amp/phi sine. Same WLS as `migkas-geometry-fit.py`.

| case | $d$ | $L$ (fixed) | $D$ | $L/D$ | $m_c$ | $\ell_0$ | $R^2$ unw. | $\chi^2$ | $A_{\rm amp}$ | $\delta(R_{\rm vis})$ | $H t$ | age (Gyr) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| d14.3 | 14.30 | 922.6 | 46.0 | 20.062 | 4.531 | 20.9 | 0.7160 | 10.57 | 0.1962 | 0.0310 | 0.139 | 1.95 |
| d3 | 3.0 | 193.6 | 45.5 | 4.256 | 3.309 | 20.8 | 0.7184 | 10.49 | 0.1962 | 0.1479 | 0.461 | 6.44 |

Comparison: ptshell WLS $R^2\approx 0.546$ ($D\approx 58.2$); free-$L$ exp $R^2\approx 0.759$ ($L\approx 6.07$ Gpc, $\delta(14.3)\gg 3\%$).

**Read carefully:** with CatWISE-scale $L\gg D$ the density gradient alone cannot source Migkas-sized $A/A_{\rm all}$ variation. The optimizer recovers $R^2\sim 0.72$ by driving **large $m_c$** (deep CoE potential), which yields young ages ($\sim 2$–$6$ Gyr) and low $Ht$. That is a dynamical workaround, not a density-dipole solution. Free-$L$ exp fits Migkas with $L\sim 6$ Gpc but then $\delta(14.3)\sim 10$, far above CatWISE. Locking $L$ to $\delta\approx 0.031$ makes that trade-off explicit.

## Files

| Path | Content |
|---|---|
| `catwise-exp-density.png` | $\rho(r)$ with $D\pm d$, Earth tangent / $(d\rho/dr)_D$ annotated |
| `catwise-exp-density.csv` | radial $\rho/\rho_{\rm crit}$ columns |
| `catwise-exp-density.md` | this note |
| `catwise-exp-migkas-fit-d14.3.png` / `.csv` | Migkas curve, primary $d\approx 14.3$ |
| `catwise-exp-migkas-fit-d3.png` / `.csv` | Migkas curve, $d=3$ |
| `catwise-exp-migkas-fit.json` | machine-readable fits + density params |
| `catwise-exp-density-migkas.py` | script |
