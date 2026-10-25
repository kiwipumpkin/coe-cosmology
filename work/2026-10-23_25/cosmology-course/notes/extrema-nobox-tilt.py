"""Extremum angles and heights of A/A_all = T (H_obs/H)/A_all with the sine tilt
T = 1 - Amp_A sin(l + phi_A) and NO boxcar. First order in d/D and Amp_A, with the
Amp_A*d/D cross terms kept; heights carry the half-shift term. Written in
D, L, d, M, G, H, alpha. Checked in nobox_tilt_check.py. Nothing fitted."""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from nobox_check import exact_ext, K0
from nobox_tilt_check import tilt_formulas, AMP, PH
plt.rcParams.update({"mathtext.fontset": "cm", "font.family": "serif"})

lf, Af, Aall_f, meanT, ds = tilt_formulas(1.0, 1/5.9, 0.1, K0)
ex, hbar, _ = exact_ext(5.9, 0.1, 1.0, K0, amp=AMP)
mx = [e for e in ex if e[0] == "max"]; mn = [e for e in ex if e[0] == "min"]
le = [mx[0][1], mx[1][1], mn[0][1], mn[1][1]]; Ae = [mx[0][2], mx[1][2], mn[0][2], mn[1][2]]
print("formula", np.round(lf, 2), np.round(Af, 4), "exact", np.round(le, 2), np.round(Ae, 4))

Sx = r"(1-f_u)\left(6L^{3}e^{D/L}-6L^{3}-6DL^{2}-3D^{2}L\right)+f_uD^{3}"
iS = r"\dfrac{1}{" + Sx + r"}"
mxb = r"\left(\dfrac{1}{D^{3}}-\dfrac{2L+(1-f_u)D}{2L\left[" + Sx + r"\right]}\right)"
T = [
 (r"$\bf{Angles\ and\ heights\ of\ the\ maxima\ and\ minima\ of\ }A/A_{\rm all}\bf{\ with\ the\ sine\ tilt,\ no\ boxcar}$", 16),
 (r"$D$: distance from the centre of expansion (CoE) to the observer;  $L$: scale length of the core density $e^{-r/L}$;", 12.5),
 (r"$d$: radius of the source shell around the observer;  $M$: mass inside $D$;  $G$: Newton's constant;  $H$: background expansion rate;", 12.5),
 (r"$\alpha$: peculiar-velocity scale;  $f_u=0.94$: share of the density at $D$ from the uniform floor;  $\ell$: angle from the anti-CoE direction;", 12.5),
 (r"Amp$_A=-0.06$ and $\phi_A=0.35$ rad: the tilt constants.", 12.5),
 ("gap", 0.2),
 (r"$\bf{Model}$:   $\dfrac{A}{A_{\rm all}}=\left[1-\mathrm{Amp}_A\sin(\ell+\phi_A)\right]\dfrac{H_{\rm obs}/H}{A_{\rm all}}$,   with $H_{\rm obs}/H$ to first order in $d/D$:", 14.6),
 (r"$\dfrac{H_{\rm obs}}{H}\approx1-\dfrac{\alpha GM}{D^{3}H^{2}}+\dfrac{3\alpha GM}{H^{2}}\left(\dfrac{1}{D^{3}}-" + iS + r"\right)\cos^{2}\ell$", 15),
 (r"$\qquad+\dfrac{3\alpha GM\,d}{2DH^{2}}\left[3\left(\dfrac{1}{D^{3}}-" + iS + r"\right)\cos\ell\right.$", 15),
 (r"$\qquad\qquad\left.+\left(\dfrac{(1-f_u)\,D}{L\left[" + Sx + r"\right]}-\dfrac{5}{D^{3}}+\dfrac{5}{" + Sx + r"}\right)\cos^{3}\ell\right]$", 15),
 ("gap", 0.1),
 (r"$\bf{Sky\ mean}$ of $H_{\rm obs}/H$ over $\ell$ (the model's normalisation)", 14),
 (r"$A_{\rm all}\approx1+\dfrac{\alpha GM}{2H^{2}}\left(\dfrac{1}{D^{3}}-\dfrac{3}{" + Sx + r"}\right)$", 15),
 ("gap", 0.15),
 (r"$\bf{Tilt\ shifts}$ (radians; upper sign $\delta_1$ or $\delta'_1$, lower sign $\delta_2$ or $\delta'_2$)", 14),
 (r"$\delta_{1,2}\approx-\mathrm{Amp}_A\cos\phi_A\times\left\{\dfrac{H^{2}}{\alpha GM}+\dfrac{2}{D^{3}}-\dfrac{3}{" + Sx + r"}\right.$", 15),
 (r"$\qquad\qquad\left.\mp\dfrac{3d}{D}" + mxb + r"\right\}$", 15),
 (r"$\qquad\quad\div\left\{\dfrac{6}{D^{3}}-\dfrac{6}{" + Sx + r"}\right.$", 15),
 (r"$\qquad\qquad\left.\mp\dfrac{9d}{2D}\left(\dfrac{4}{D^{3}}-\dfrac{4L+(1-f_u)D}{L\left[" + Sx + r"\right]}\right)\right\}$", 15),
 (r"$\delta'_{1,2}\approx-\mathrm{Amp}_A\left(\sin\phi_A\pm\dfrac{3d}{4D}\cos\phi_A\right)\left(\dfrac{H^{2}}{\alpha GM}-\dfrac{1}{D^{3}}\right)\div\left(\dfrac{6}{D^{3}}-\dfrac{6}{" + Sx + r"}\right)$", 15),
 ("gap", 0.15),
 (r"$\bf{Angles}$", 14),
 (r"$\ell_{\max,1}\approx\dfrac{180^\circ}{\pi}\,\delta_1,\qquad\ell_{\max,2}\approx180^\circ-\dfrac{180^\circ}{\pi}\,\delta_2$", 15),
 (r"$\ell_{\min,1}\approx90^\circ+\dfrac{135^\circ}{\pi}\,\dfrac{d}{D}+\dfrac{180^\circ}{\pi}\,\delta'_1,\qquad\ell_{\min,2}\approx270^\circ-\dfrac{135^\circ}{\pi}\,\dfrac{d}{D}-\dfrac{180^\circ}{\pi}\,\delta'_2$", 15),
 ("gap", 0.15),
 (r"$\bf{Heights}$ (upper sign max,1 / min,1; lower sign max,2 / min,2)", 14),
 (r"$\left(\dfrac{A}{A_{\rm all}}\right)_{\max,1,2}\approx\left(1\mp\mathrm{Amp}_A\sin\phi_A-\dfrac{1}{2}\mathrm{Amp}_A\cos\phi_A\,\delta_{1,2}\right)$", 15),
 (r"$\qquad\qquad\times\dfrac{1}{A_{\rm all}}\left\{1+\dfrac{\alpha GM}{H^{2}}\left[\dfrac{2}{D^{3}}-\dfrac{3}{" + Sx + r"}\right.\right.$", 15),
 (r"$\qquad\qquad\qquad\left.\left.\mp\dfrac{3d}{D}" + mxb + r"\right]\right\}$", 15),
 (r"$\left(\dfrac{A}{A_{\rm all}}\right)_{\min,1,2}\approx\left(1\mp\mathrm{Amp}_A\cos\phi_A+\dfrac{3d}{4D}\mathrm{Amp}_A\sin\phi_A+\dfrac{1}{2}\mathrm{Amp}_A\sin\phi_A\,\delta'_{1,2}\right)$", 15),
 (r"$\qquad\qquad\times\dfrac{1}{A_{\rm all}}\left(1-\dfrac{\alpha GM}{D^{3}H^{2}}\right)$", 15),
 ("gap", 0.3),
 (r"$\bf{Notes}$", 14),
 (r"$\bullet$  No boxcar. $\delta_1,\delta_2,\delta'_1,\delta'_2$ are in radians and $180^\circ/\pi$ turns them into degrees;  Amp$_A\cos\phi_A=-0.0564$,  Amp$_A\sin\phi_A=-0.0206$.", 12.5),
 (r"$\bullet$  First order in $d/D$ and in Amp$_A$, keeping the Amp$_A\times d/D$ cross terms (the $d$ parts of $\delta_{1,2}$ and $\delta'_{1,2}$, and $\frac{3d}{4D}$Amp$_A\sin\phi_A$).", 12.5),
 (r"$\bullet$  Those cross terms bring $\ell_{\max,2}$ and both minima to within 0.1$^\circ$ at baseline; $\ell_{\max,1}$ stays 1.1$^\circ$ high from $(d/D)^2$ terms this order drops.", 12.5),
 (r"$\bullet$  The $\frac{1}{2}$Amp$_A\,\delta$ terms in the heights: at the shifted peak the tilt factor rises by Amp$_A\times$shift and $H_{\rm obs}/H$ falls by half of that.", 12.5),
 (r"$\bullet$  They cut the height error at the maxima from 0.006 to 0.001 at baseline.", 12.5),
 (r"$\bullet$  Sky mean: the tilt's $\sin\ell$ part averages to zero; its product with the $\cos\ell$ part of $H_{\rm obs}/H$ changes the mean of the tilted curve by", 12.5),
 (r"$\qquad-\dfrac{9}{16}\,\mathrm{Amp}_A\sin\phi_A\,\dfrac{\alpha GM\,d}{DH^{2}}\left[\dfrac{L+(1-f_u)D}{L\left[" + Sx + r"\right]}-\dfrac{1}{D^{3}}\right]$", 13.5),
 (r"$\qquad$This is $%.1f\times10^{-5}$ at baseline (second order), and the model divides by the mean of $H_{\rm obs}/H$, so $A_{\rm all}$ above is unchanged." % ((meanT - Aall_f)*1e5), 12.5),
 (r"$\bullet$  The bracket $" + Sx + r"$ equals $3M/4\pi\rho(D)$, with $\rho(D)$ the density at the observer;  $\alpha$ only appears as $\alpha GM$.", 12.5),
 (r"$\bullet$  Baseline: $D/L$=5.9, $d/D$=0.1, $\alpha GM/D^{3}H^{2}$=0.2306.", 12.5),
 (r"$\bullet$  Baseline angles (max,1, max,2, min,1, min,2): %.1f$^\circ$, %.1f$^\circ$, %.1f$^\circ$, %.1f$^\circ$ (exact %.1f, %.1f, %.1f, %.1f)."
  % (*lf, *le), 12.5),
 (r"$\bullet$  Baseline heights, same order:  %.4f, %.4f, %.4f, %.4f  (exact %.4f, %.4f, %.4f, %.4f);  $A_{\rm all}$ = %.4f (exact %.4f)." % (*Af, *Ae, Aall_f, hbar), 12.5),
 (r"$\bullet$  $D/L$: max,1 within 2$^\circ$ for $D/L\geq4.5$ and within 0.7$^\circ$ for $D/L\geq8$; the other three within 0.7$^\circ$ for $D/L\geq4.3$; heights within 0.008.", 12.5),
 (r"$\bullet$  $\alpha$ and $M$: max,1 within 2$^\circ$ for $\alpha GM/D^{3}H^{2}\geq0.16$, the others within 0.4$^\circ$, heights within 0.012; fails for $\alpha GM>D^{3}H^{2}$ ($H_{\rm obs}<0$).", 12.5),
 (r"$\bullet$  $d/D$ (at $\alpha$=1.3): max,1 within 2$^\circ$ up to $d/D\approx0.17$, fails near 0.3 (the $0^\circ$ peak flattens); max,2 within 1$^\circ$ up to 0.4, minima up to 0.5.", 12.5),
 (r"$\bullet$  Below $D/L\approx2.36$ or $\alpha GM/D^{3}H^{2}\approx0.044$ the tilt wins: max,2 and min,1 merge and vanish, and the formulas no longer apply.", 12.5),
]
Hh = {16: 0.55, 15: 0.95, 13.5: 0.85, 14.6: 0.95, 14: 0.5, 12.5: 0.36}
tot = sum(v if t == "gap" else Hh[v] for t, v in T)
fig = plt.figure(figsize=(13, tot + 0.4)); y = tot + 0.2
for t, v in T:
    if t == "gap": y -= v; continue
    fig.text(0.01, y/(tot + 0.4), t, fontsize=(14 if v == 14.6 else v), va="top"); y -= Hh[v]
fig.savefig("extrema-nobox-tilt.png", dpi=200, bbox_inches="tight", pad_inches=0.25)
from PIL import Image; im = Image.open("extrema-nobox-tilt.png"); print("in", im.size[0]/200, im.size[1]/200)
