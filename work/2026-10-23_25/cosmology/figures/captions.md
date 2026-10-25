# Figure captions — Sedov–Taylor centered-universe model

**Ontology (binding):** CoE beyond the horizon; observer off-center at distance \(D\); decreasing \(\\rho(r)\) (Phase A eq. 3.13 exponential); Gauss \(g=-GM(<r)/r^2\). Pure homologous ST \(\Rightarrow\) isotropic \(H_{\\mathrm{obs}}\); gravity restores anisotropy. At \(d\\ll D\) the leading sky pattern is a quadrupole \(P_2\); \(\\kappa_{\\mathrm{natural}}=\\Omega_m/2=0.15\). **Success criterion:** first-principles $H_{\\mathrm{obs}}(\hat n)$ / mapped $A(\ell)$ matches Migkas angular pattern — not a best-fit dipole.

| File | One-line caption |
|------|------------------|
| `geometry.png` | 2D schematic: CoE at the origin, observer at \((D,0)\), sample sources on a circle of radius \(R\\sim d\) about the observer, with \(\\hat{D}\), \(d\), \(r\), and \(\\varphi\) labeled. |
| `density_mass.png` | Phase A eq. (3.13): \(\\rho(r)/\\rho_0=e^{-r/r_s}\) and cumulative \(m(r)/M\), with an example mark at \(r=D\) (\(D/r_s=3\)). |
| `Hobs_vs_angle.png` | Migkas+2020 Fig.6 digitised $A/A_{\mathrm{all}}$ (+ $B$) vs Galactic longitude overlaid with ST+Gauss model ($\kappa=0.15$, axis $\ell=280^\circ$); no $A\cos\varphi$ dipole. |
| `multipoles.png` | Full-sky Legendre coefficients \(C_0\\ldots C_3\) at \(d/D=0.1\), \(\\kappa=0.15\): \(C_2\\approx 0.30\) dominant, \(C_1\\approx 0\). |
| `sky_mollweide.png` | Mollweide map of the leading \(P_2\) anisotropy with preferred axis \(\\hat{D}\) (quadrupole lobes, not a dipole). |
| `kappa_amplitude.png` | Half peak-to-peak anisotropy vs \(\\kappa\), with \(\\kappa_{\\mathrm{natural}}=0.15\) marked. |

**Parameters used unless noted:** \(\\alpha=1\), \(\\kappa_{\\mathrm{natural}}=0.15\), \(d/D=0.1\); tidal series \(H_{\\mathrm{obs}}/H-1=\\alpha\\kappa(3\\mu^2-1)+\\alpha\\kappa\\varepsilon(4.5\\mu-7.5\\mu^3)\).

**Regenerate:** `python3 make_figures.py`
