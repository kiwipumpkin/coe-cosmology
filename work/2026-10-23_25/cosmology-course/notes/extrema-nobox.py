"""Extremum angles and heights of A/A_all with NO boxcar and NO fitted sine
rotation: A/A_all = h(l)/<h>, h = H_obs/H at one input angle. Expanded to
first order in d/D, solved dA/dl = 0, written in D, L, d, M, G, H, alpha.
Checked against the exact no-boxcar curve in nobox_check.py. Nothing fitted."""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from nobox_check import exact_ext, K0
plt.rcParams.update({"mathtext.fontset": "cm", "font.family": "serif"})
FU = 0.94

def card_formulas(D, L, d, M, al=1.0, G=1.0, H=1.0):
    S = (1-FU)*(6*L**3*np.exp(D/L) - 6*L**3 - 6*D*L**2 - 3*D**2*L) + FU*D**3
    g = al*G*M/H**2
    Aall = 1 + g/2*(1/D**3 - 3/S)
    pk = lambda sg: 1 + g*(2/D**3 - 3/S + sg*3*d/D*(1/D**3 - (2*L + (1-FU)*D)/(2*L*S)))
    lmin = 90 + np.degrees(3*d/(4*D))
    return (0.0, 180.0, lmin, 360 - lmin), (pk(-1)/Aall, pk(+1)/Aall, (1 - g/D**3)/Aall, (1 - g/D**3)/Aall), Aall

lf, Af, Aall_f = card_formulas(1.0, 1/5.9, 0.1, K0)
ex, hbar, _ = exact_ext(5.9, 0.1, 1.0, K0)
le = [e[1] for e in ex]; Ae = [e[2] for e in ex]   # order: max0, min1, max180, min2
print("formula l", lf, "A", Af, "Aall", Aall_f)
print("exact   l", le, "A", Ae, "Aall", hbar)

Sx = r"(1-f_u)\left(6L^{3}e^{D/L}-6L^{3}-6DL^{2}-3D^{2}L\right)+f_uD^{3}"
iS = r"\dfrac{1}{" + Sx + r"}"
T = [
 (r"$\bf{Angles\ and\ heights\ of\ the\ maxima\ and\ minima\ of\ }A/A_{\rm all}\bf{,\ no\ boxcar}$", 16),
 (r"$D$: distance from the centre of expansion (CoE) to the observer;  $L$: scale length of the core density $e^{-r/L}$;", 12.5),
 (r"$d$: radius of the source shell around the observer;  $M$: mass inside $D$;  $G$: Newton's constant;  $H$: background expansion rate;", 12.5),
 (r"$\alpha$: peculiar-velocity scale;  $f_u=0.94$: share of the density at $D$ from the uniform floor;  $\ell$: angle from the anti-CoE direction.", 12.5),
 ("gap", 0.2),
 (r"$\bf{Input\ expression\ to\ first\ order\ in\ }d/D$   (your $H_{\rm obs}/H$ at one angle $\ell$)", 14),
 (r"$\dfrac{H_{\rm obs}}{H}\approx1-\dfrac{\alpha GM}{D^{3}H^{2}}+\dfrac{3\alpha GM}{H^{2}}\left(\dfrac{1}{D^{3}}-" + iS + r"\right)\cos^{2}\ell$", 15),
 (r"$\qquad+\dfrac{3\alpha GM\,d}{2DH^{2}}\left[3\left(\dfrac{1}{D^{3}}-" + iS + r"\right)\cos\ell\right.$", 15),
 (r"$\qquad\qquad\left.+\left(\dfrac{(1-f_u)\,D}{L\left[" + Sx + r"\right]}-\dfrac{5}{D^{3}}+\dfrac{5}{" + Sx + r"}\right)\cos^{3}\ell\right]$", 15),
 (r"Setting $dA/d\ell=0$:   $\sin\ell=0$  (the maxima),   or   $\cos\ell=-\dfrac{3d}{4D}$  (the minima).", 14),
 ("gap", 0.2),
 (r"$\bf{Sky\ mean}$ of $H_{\rm obs}/H$ over $\ell$", 14),
 (r"$A_{\rm all}\approx1+\dfrac{\alpha GM}{2H^{2}}\left(\dfrac{1}{D^{3}}-\dfrac{3}{" + Sx + r"}\right)$", 15),
 ("gap", 0.15),
 (r"$\bf{Maxima}$ at $\ell=0^\circ$ (upper sign) and $\ell=180^\circ$ (lower sign; toward the CoE)", 14),
 (r"$\left(\dfrac{A}{A_{\rm all}}\right)_{\max}\approx\dfrac{1}{A_{\rm all}}\left\{1+\dfrac{\alpha GM}{H^{2}}\left[\dfrac{2}{D^{3}}-\dfrac{3}{" + Sx + r"}\right.\right.$", 15),
 (r"$\qquad\qquad\qquad\left.\left.\mp\dfrac{3d}{D}\left(\dfrac{1}{D^{3}}-\dfrac{2L+(1-f_u)D}{2L\left[" + Sx + r"\right]}\right)\right]\right\}$", 15),
 ("gap", 0.15),
 (r"$\bf{Minima}$ at $\ell_{\min,1}$ (near $90^\circ$) and $\ell_{\min,2}$ (near $270^\circ$)", 14),
 (r"$\ell_{\min,1}\approx90^\circ+\dfrac{135^\circ}{\pi}\,\dfrac{d}{D},\qquad\ell_{\min,2}\approx270^\circ-\dfrac{135^\circ}{\pi}\,\dfrac{d}{D}\qquad\left(\dfrac{135^\circ}{\pi}=42.97^\circ\right)$", 15),
 (r"$\left(\dfrac{A}{A_{\rm all}}\right)_{\min,1}=\left(\dfrac{A}{A_{\rm all}}\right)_{\min,2}\approx\dfrac{1}{A_{\rm all}}\left(1-\dfrac{\alpha GM}{D^{3}H^{2}}\right)$", 15),
 ("gap", 0.3),
 (r"$\bf{Notes}$", 14),
 (r"$\bullet$  No boxcar and no sine rotation, so $A/A_{\rm all}=(H_{\rm obs}/H)\,/\,A_{\rm all}$ at a single angle $\ell$.", 12.5),
 (r"$\bullet$  The fitted rotation $1-\mathrm{Amp}_A\sin(\ell+\phi_A)$ is left out: it is the Migkas one-cycle phase factor, fitted afterwards, not part of the input expression.", 12.5),
 (r"$\bullet$  The density is $\rho(r)=\rho_{\rm core}\,e^{-r/L}+\rho_{\rm floor}$, and the bracket $" + Sx + r"$ equals $3M/4\pi\rho(D)$.", 12.5),
 (r"$\bullet$  So $3\alpha GM/(H^{2}\times\mathrm{bracket})=4\pi\alpha G\rho(D)/H^{2}$ (density at the observer), and $\alpha GM/D^{3}H^{2}=4\pi\alpha G\bar\rho(<D)/3H^{2}$ with $\bar\rho(<D)=3M/4\pi D^{3}$.", 12.5),
 (r"$\bullet$  $\alpha$ only multiplies $GM$, so the curve sees $\alpha M$, not $\alpha$ and $M$ separately.", 12.5),
 (r"$\bullet$  The maxima at $0^\circ$ and $180^\circ$ are exact: without the boxcar and rotation the curve is mirror-symmetric about the CoE axis.", 12.5),
 (r"$\bullet$  The minima shift by $42.97^\circ\times d/D$ from $90^\circ$ and $270^\circ$; at this order the shift does not depend on $L$, $M$ or $\alpha$.", 12.5),
 (r"$\bullet$  Baseline ($D/L$=5.9, $d/D$=0.1, $\alpha GM/D^{3}H^{2}$=0.2306) angles: %.1f$^\circ$, %.1f$^\circ$, %.1f$^\circ$, %.1f$^\circ$  (exact no-boxcar %.1f, %.1f, %.1f, %.1f)."
  % (lf[0], lf[2], lf[1], lf[3], le[0], le[1], le[2], le[3]), 12.5),
 (r"$\bullet$  Baseline heights (same order): %.4f, %.4f, %.4f, %.4f  (exact %.4f, %.4f, %.4f, %.4f);  $A_{\rm all}$ = %.4f  (exact %.4f)."
  % (Af[0], Af[2], Af[1], Af[3], Ae[0], Ae[1], Ae[2], Ae[3], Aall_f, hbar), 12.5),
 (r"$\bullet$  $D/L$ from 1 to 15: angles within 0.02$^\circ$, heights within 0.008. The merges seen with the boxcar do not occur.", 12.5),
 (r"$\bullet$  $\alpha$ and $M$: same accuracy (heights within 0.012) while $\alpha GM<D^{3}H^{2}$; above that $H_{\rm obs}$ goes negative at the minima and the formulas fail.", 12.5),
 (r"$\bullet$  $d/D$: heights within 0.01 up to $d/D\approx0.22$ and within 0.05 up to 0.5; minimum angles within 1$^\circ$ up to $d/D\approx0.48$.", 12.5),
]
Hh = {16: 0.55, 15: 0.95, 14: 0.5, 12.5: 0.36}
tot = sum(v if t == "gap" else Hh[v] for t, v in T)
fig = plt.figure(figsize=(13, tot + 0.4)); y = tot + 0.2
for t, v in T:
    if t == "gap": y -= v; continue
    fig.text(0.01, y/(tot + 0.4), t, fontsize=v, va="top"); y -= Hh[v]
fig.savefig("extrema-nobox.png", dpi=200, bbox_inches="tight", pad_inches=0.25)
from PIL import Image; im = Image.open("extrema-nobox.png"); print("px", im.size, "in", im.size[0]/200, im.size[1]/200)
