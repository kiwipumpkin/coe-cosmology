# CatWISE δ=3% density level set + Migkas refit

**Date:** 2026-10-06 (America/Chicago)

**Acronyms:** CatWISE = CatWISE quasar dipole (Secrest et al. 2021); CoE = centre of expansion; Migkas = Migkas et al. cluster $A/A_{\mathrm{all}}$ anisotropy.

## Constraint (only fixed)

$$\delta=\frac{\rho(D-d)-\rho(D+d)}{\rho(D)}=0.031\quad\text{at primary }d=R_{\rm vis}\approx 14.3\,\mathrm{Gpc}.$$

- $\rho(D)=\Omega_m\rho_{\rm crit}$ with $\Omega_m=0.3$.
- Density-plot baseline $D=58.2$ Gpc; Migkas fit frees $D$.
- Smooth entire way; $\rho(0)$ finite (flag if above stellar-BH density; no floor).
- Prefer $\ell_0$ nearer ptshell $\sim 15.2^\circ$ (Chris: $19.9^\circ$ too much): report free $\ell_0$ and $\ell_0=15.2^\circ$ fixed.

## Density family (level set at $D=58.2$ Gpc)

| law | formula | shape | $L$ (Gpc) | $\delta$ check | $\rho(0)/\rho_c$ | $(\rho'/\rho)_D$ | cliff index | gentle? |
|---|---|---|---:|---:|---:|---:|---:|---|
| exp | `ρ=ρ₀ e^{−r/L}` | — | 922.6 | 0.03100 | 0.3195 | -0.00108384 | 0 | yes |
| pl | `ρ=ρ₀ (1+r/L)^{−1.2}` | n=1.2 | 1049 | 0.03100 | 0.3201 | -0.00108367 | 0.0129 | yes |
| cored | `ρ=ρ₀/(1+(r/L)²)` | — | 321.9 | 0.03100 | 0.3098 | -0.00108776 | 0.23 | yes |
| gauss | `ρ=ρ₀ e^{−½(r/L)²}` | — | 231.5 | 0.03100 | 0.3096 | -0.00108591 | 0.246 | yes |
| soft | `ρ=ρ₀/(1+(r/L)^0.3)^{1.0}` | alpha=0.3, beta=1 | 5110 | 0.03100 | 0.3784 | -0.00106752 | 0.19 | yes |
| eta | `ρ=ρ₀ e^{−(r/L)^0.25}` | eta=0.25 | 1.525e+04 | 0.03100 | 0.3847 | -0.00106772 | 0.187 | yes |

**Cliff index** $=|(\rho'/\rho)_{D+d/2}-(\rho'/\rho)_{D-d/2}|$$/|(\rho'/\rho)_D|$: small ⇒ slope varies gently through Earth (no precipice). Under the δ lock all laws are near-linear through Earth (small index). The η≈0.05 edge's cliff sits in the interior, so it is not captured here and is excluded by fiat.

## Ranking under $\ell_0=15.2^\circ$ fixed

| rank | law | $R^2$ | $\chi^2$ | $n_{\rm free}$ | $\Delta\chi^2$ | $\sigma$ null (nested) | $\sqrt{\Delta\chi^2}$ | $D$ | $L$ | shape | $m_c$ | $\ell_0$ | cliff | beats 0.72? |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|---|
| 1 | eta | 0.7428 | 10.395 | 3 | 34.742 | 5.14 | 5.89 | 38.47 | 5.701e+04 | eta=0.25 | 4.365 | 15.20 | 0.287 | YES |
| 2 | soft | 0.7428 | 10.425 | 3 | 34.711 | 5.14 | 5.89 | 39.11 | 1.816e+04 | alpha=0.3, beta=1 | 4.362 | 15.20 | 0.278 | YES |
| 3 | pl | 0.7178 | 11.582 | 3 | 33.554 | 5.03 | 5.79 | 47.86 | 1060 | n=1.2 | 4.504 | 15.20 | 0.0129 | no |
| 4 | exp | 0.7173 | 11.614 | 2 | 33.522 | 5.32 | 5.79 | 48.26 | 922.6 | — | 4.51 | 15.20 | 0 | no |
| 5 | cored | 0.6933 | 12.829 | 2 | 32.307 | 5.21 | 5.68 | 63.33 | 335.4 | — | 4.657 | 15.20 | 0.21 | no |
| 6 | gauss | 0.6928 | 12.850 | 2 | 32.286 | 5.20 | 5.68 | 63.37 | 241.6 | — | 4.664 | 15.20 | 0.226 | no |

**Best under fixed15.2:** `eta` with $R^2=0.7428$, $\chi^2=10.395$, $\ell_0=15.20^\circ$, cliff=0.287.

## Ranking under free $\ell_0$

| rank | law | $R^2$ | $\chi^2$ | $n_{\rm free}$ | $\Delta\chi^2$ | $\sigma$ null (nested) | $\sqrt{\Delta\chi^2}$ | $D$ | $L$ | shape | $m_c$ | $\ell_0$ | cliff | beats 0.72? |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|---|
| 1 | eta | 0.7394 | 9.565 | 4 | 35.572 | 4.96 | 5.96 | 37.23 | 6.35e+04 | eta=0.25 | 4.385 | 20.08 | 0.297 | R²>0.72 but ℓ₀>16° |
| 2 | soft | 0.7383 | 9.593 | 4 | 35.544 | 4.96 | 5.96 | 37.25 | 2.112e+04 | alpha=0.3, beta=1 | 4.39 | 20.13 | 0.292 | R²>0.72 but ℓ₀>16° |
| 3 | pl | 0.7168 | 10.546 | 4 | 34.591 | 4.87 | 5.88 | 45.91 | 1061 | n=1.2 | 4.524 | 20.91 | 0.0129 | no |
| 4 | exp | 0.7157 | 10.569 | 3 | 34.567 | 5.12 | 5.88 | 45.80 | 922.6 | — | 4.532 | 20.97 | 0 | no |
| 5 | cored | 0.6936 | 11.596 | 3 | 33.540 | 5.03 | 5.79 | 59.50 | 325.4 | — | 4.677 | 21.70 | 0.225 | no |
| 6 | gauss | 0.6936 | 11.611 | 3 | 33.525 | 5.03 | 5.79 | 59.79 | 234.7 | — | 4.683 | 21.64 | 0.239 | no |

**Best under free:** `eta` with $R^2=0.7394$, $\chi^2=9.565$, $\ell_0=20.08^\circ$, cliff=0.297.

## σ vs null for every fit

Null: flat $A/A_{\rm all}=1$, $\chi^2_{\rm null}=45.137$, $N=71$. Nested: $\Delta\chi^2=\chi^2_{\rm null}-\chi^2$, $p=P(\chi^2_{\Delta\nu}>\Delta\chi^2)$, $\Delta\nu=n_{\rm free}$ (D, $m_c$, [$\ell_0$], [shape]); $\sigma=\Phi^{-1}(1-p)$ — same method as `catwise-sigma-null.py`. Fixing $\ell_0=15.2^\circ$ removes one free parameter.

| law | $\ell_0$ mode | free params | $n_{\rm free}$ | $\chi^2$ | $\Delta\chi^2$ | $p$ | $\sigma$ nested | $\sqrt{\Delta\chi^2}$ | $R^2$ |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| cored | fixed15.2 | D,mc | 2 | 12.829 | 32.307 | 9.650e-08 | 5.21 | 5.68 | 0.6933 |
| cored | free | D,mc,l0 | 3 | 11.596 | 33.540 | 2.477e-07 | 5.03 | 5.79 | 0.6936 |
| eta | fixed15.2 | D,mc,eta | 3 | 10.395 | 34.742 | 1.381e-07 | 5.14 | 5.89 | 0.7428 |
| eta | free | D,mc,l0,eta | 4 | 9.565 | 35.572 | 3.544e-07 | 4.96 | 5.96 | 0.7394 |
| exp | fixed15.2 | D,mc | 2 | 11.614 | 33.522 | 5.257e-08 | 5.32 | 5.79 | 0.7173 |
| exp | free | D,mc,l0 | 3 | 10.569 | 34.567 | 1.504e-07 | 5.12 | 5.88 | 0.7157 |
| gauss | fixed15.2 | D,mc | 2 | 12.850 | 32.286 | 9.753e-08 | 5.20 | 5.68 | 0.6928 |
| gauss | free | D,mc,l0 | 3 | 11.611 | 33.525 | 2.495e-07 | 5.03 | 5.79 | 0.6936 |
| pl | fixed15.2 | D,mc,n | 3 | 11.582 | 33.554 | 2.461e-07 | 5.03 | 5.79 | 0.7178 |
| pl | free | D,mc,l0,n | 4 | 10.546 | 34.591 | 5.638e-07 | 4.87 | 5.88 | 0.7168 |
| soft | fixed15.2 | D,mc,alpha | 3 | 10.425 | 34.711 | 1.402e-07 | 5.14 | 5.89 | 0.7428 |
| soft | free | D,mc,l0,alpha | 4 | 9.593 | 35.544 | 3.591e-07 | 4.96 | 5.96 | 0.7383 |

## Best R² per law (one line each; σ = nested vs null, Δν = n_free)

| law | best R² | σ | χ² | n_free | ℓ₀ | D (Gpc) | L (Gpc) | shape | gentle / caveat |
|---|---:|---:|---:|---:|---:|---:|---:|---|---|
| soft | 0.7428 | 5.14 | 10.425 | 3 | 15.2 | 39.1 | 1.816e+04 | alpha=0.3, beta=1 | gentle through Earth; α pinned at 0.3 → central cusp |
| eta | 0.7428 | 5.14 | 10.395 | 3 | 15.2 | 38.5 | 5.701e+04 | eta=0.25 | η pinned at 0.25 floor (edge, cusp) — not a winner |
| pl | 0.7178 | 5.03 | 11.582 | 3 | 15.2 | 47.9 | 1060 | n=1.2 | gentle; n pinned at 1.2 bound |
| exp | 0.7173 | 5.32 | 11.614 | 2 | 15.2 | 48.3 | 922.6 | — | gentle (near-linear through Earth); tiny central kink |
| cored | 0.6936 | 5.03 | 11.596 | 3 | 21.7 | 59.5 | 325.4 | — | gentle, smooth centre |
| gauss | 0.6936 | 5.03 | 11.611 | 3 | 21.6 | 59.8 | 234.7 | — | gentle, smooth centre |

## Edge caveat (important)

Every free shape parameter runs to its flattest bound (pl n→1.2, soft α→0.3, η→0.25). Under the δ lock, the Migkas fit keeps preferring a profile that is almost flat through Earth and steepens toward the CoE. That is the same pull that drove η→0.05. Soft α profile at ℓ₀=15.2° (β=1, {D, m_c} refit):

| α | R² | χ² | σ (ν=3) | D (Gpc) | centre |
|---:|---:|---:|---:|---:|---|
| 0.31 | 0.7422 | 10.443 | 5.13 | 39.1 | cusp: rho'(0) -> -inf |
| 0.50 | 0.7352 | 10.775 | 5.10 | 41.5 | cusp: rho'(0) -> -inf |
| 0.70 | 0.7281 | 11.110 | 5.07 | 44.1 | cusp: rho'(0) -> -inf |
| 1.00 | 0.7180 | 11.576 | 5.03 | 47.9 | kink: rho'(0) finite <0 |
| 1.50 | 0.7039 | 12.258 | 4.97 | 54.9 | rho'(0)=0 but rho''(0) diverges |
| 2.00 | 0.6933 | 12.829 | 4.91 | 63.3 | smooth: rho'(0)=0, rho'' finite |

R² climbs steadily as α falls; it passes 0.72 at α≈0.9. α=2 (smooth core) is identical to cored (0.693). So R²>0.72 at ℓ₀=15.2° comes only from cusped centres (α<1), not from a smooth core.

## Comparison priors

| model | $R^2$ | $\chi^2$ | $\ell_0$ | note |
|---|---:|---:|---:|---|
| ptshell WLS | 0.546 | 18.7 | 15.2 | geometry ptshell |
| CatWISE-$L$ exp | 0.716 | 10.6 | 20.9 | $\delta$ lock, $\eta=1$ |
| CatWISE-η edge | 0.745 | 9.3 | 19.9 | η≈0.05 cliff — Chris dislikes |
| free-exp (no δ lock) | 0.759 | 8.5 | — | over-produces contrast |

## Recommendation

- **Winner at ℓ₀=15.2°: soft core+tail** $\rho=\rho_0/[1+(r/L)^{\alpha}]$ (β=1), α=0.3, L=1.82e4 Gpc, D=39.1 Gpc, $m_c$=4.36. Fit: R²=0.743, χ²=10.43, σ=5.14 (ν=3). It beats CatWISE-exp (0.716 free / 0.717 at 15.2°) with ℓ₀≤16°. Through Earth it is as gentle as the others: $(\rho'/\rho)_D\approx-1.07\times10^{-3}$ Gpc⁻¹, and ρ(D±d) misses the straight line by <1%. The cost: α is pinned at its 0.3 bound, giving a central cusp with ρ(0)/ρ(D)=1.26 at D=58.2. That is much softer than the η≈0.05 cliff (ρ(0)/ρ(D)=3.4), but it is still an edge solution, not an interior optimum.
- **η≥0.25** gives the same R² (0.743, χ²=10.40), also pinned at its floor. It is treated as a duplicate of soft, not a separate winner.
- **Smooth-centre laws (cored, gauss, soft α=2)** cap at R²≈0.693. **pl / exp** sit at 0.717–0.718 (pl's n is pinned at 1.2). If Chris requires a strictly smooth centre, nothing beats 0.72; soft needs α≲0.9.
- The cliff index is measured over D±d/2 only. The δ lock makes every law nearly linear there, so the indices are all small. The laws differ in the interior (centre), not at Earth.

## Files

| Path | Content |
|---|---|
| `catwise-density-levelset.png` / `.csv` / `.md` | $\rho(r)$ overlay, all $\delta=3\%$ laws |
| `catwise-density-levelset-fit.json` | machine-readable fits |
| `catwise-density-levelset-migkas-*.png` | per-law / overlay Migkas curves |
| `catwise-density-levelset-migkas.py` | this script |
| `catwise-density-levelset-soft-alpha-profile.py` / `.json` | soft α profile at ℓ₀=15.2° |
| `catwise-density-levelset-run.log` | fit log |
| `prof_core.py` (bubble) | soft kind: `(1+|x|^α)^{-β}` |
