"""Same content and math as extrema-nobox-tilt.py, re-laid out as phone-shaped
pages (10 in wide, <=13 in tall, dpi 250, equations 18 pt, notes 15 pt)."""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from nobox_check import exact_ext, K0
from nobox_tilt_check import tilt_formulas, AMP
plt.rcParams.update({"mathtext.fontset": "cm", "font.family": "serif"})

lf, Af, Aall_f, meanT, ds = tilt_formulas(1.0, 1/5.9, 0.1, K0)
ex, hbar, _ = exact_ext(5.9, 0.1, 1.0, K0, amp=AMP)
mx = [e for e in ex if e[0] == "max"]; mn = [e for e in ex if e[0] == "min"]
le = [mx[0][1], mx[1][1], mn[0][1], mn[1][1]]; Ae = [mx[0][2], mx[1][2], mn[0][2], mn[1][2]]

Sx = r"(1-f_u)\left(6L^{3}e^{D/L}-6L^{3}-6DL^{2}-3D^{2}L\right)+f_uD^{3}"
iS = r"\dfrac{1}{" + Sx + r"}"
Q, QQ, QQQ = "\x01", "\x01\x01", "\x01\x01\x01"   # indent levels (0.55 in each)
EQ, NT, HD, TT = 18, 15, 17, 20   # equation, note, header, title sizes
B = r"$\bullet$  "

pages = {}
pages[1] = [
 (r"$\bf{Maxima\ and\ minima\ of\ }A/A_{\rm all}\bf{\ with\ the\ sine\ tilt,\ no\ boxcar}$", TT),
 (r"$\bf{Page\ 1:\ definitions\ and\ model}$", HD), ("gap", 0.1),
 (B + r"$D$: distance from the centre of expansion (CoE) to the observer", NT),
 (B + r"$L$: scale length of the core density $e^{-r/L}$", NT),
 (B + r"$d$: radius of the source shell around the observer", NT),
 (B + r"$M$: mass inside $D$;   $G$: Newton's constant;   $H$: background expansion rate", NT),
 (B + r"$\alpha$: peculiar-velocity scale", NT),
 (B + r"$f_u=0.94$: share of the density at $D$ from the uniform floor", NT),
 (B + r"$\ell$: angle from the anti-CoE direction", NT),
 (B + r"Amp$_A=-0.06$ and $\phi_A=0.35$ rad: the tilt constants", NT),
 (B + r"No boxcar.", NT),
 ("gap", 0.25),
 (r"$\bf{The\ bracket}$ that appears throughout, written out in full in every equation:", NT),
 (r"$" + Sx + r"$", EQ),
 (r"It equals $3M/4\pi\rho(D)$, with $\rho(D)$ the density at the observer.", NT),
 ("gap", 0.25),
 (r"$\bf{Model}$", HD),
 (r"$\dfrac{A}{A_{\rm all}}=\left[1-\mathrm{Amp}_A\sin(\ell+\phi_A)\right]\dfrac{H_{\rm obs}/H}{A_{\rm all}}$", EQ),
 ("gap", 0.15),
 (r"$\bf{Sky\ mean}$ of $H_{\rm obs}/H$ over $\ell$ (the model's normalisation)", HD),
 (r"$A_{\rm all}\approx1+\dfrac{\alpha GM}{2H^{2}}\left(\dfrac{1}{D^{3}}\right.$", EQ),
 (QQ + r"$\left.-\dfrac{3}{" + Sx + r"}\right)$", EQ),
]
pages[2] = [
 (r"$\bf{Page\ 2:\ }H_{\rm obs}/H\bf{\ to\ first\ order\ in\ }d/D$", HD), ("gap", 0.1),
 (r"$\dfrac{H_{\rm obs}}{H}\approx1-\dfrac{\alpha GM}{D^{3}H^{2}}$", EQ),
 (Q + r"$+\dfrac{3\alpha GM}{H^{2}}\left(\dfrac{1}{D^{3}}\right.$", EQ),
 (QQ + r"$\left.-" + iS + r"\right)\cos^{2}\ell$", EQ),
 (Q + r"$+\dfrac{3\alpha GM\,d}{2DH^{2}}\left[3\left(\dfrac{1}{D^{3}}\right.\right.$", EQ),
 (QQ + r"$\left.\left.-" + iS + r"\right)\cos\ell\right.$", EQ),
 (QQ + r"$\left.+\left(\dfrac{(1-f_u)\,D}{L\left[" + Sx + r"\right]}\right.\right.$", EQ),
 (QQ + r"$\left.\left.-\dfrac{5}{D^{3}}+\dfrac{5}{" + Sx + r"}\right)\cos^{3}\ell\right]$", EQ),
]
pages[3] = [
 (r"$\bf{Page\ 3:\ the\ maxima,\ shifts\ }\delta_1,\ \delta_2$", HD),
 (r"(radians; upper sign $\delta_1$, lower sign $\delta_2$)", NT), ("gap", 0.1),
 (r"$\delta_{1,2}\approx-\mathrm{Amp}_A\cos\phi_A\times\left\{\dfrac{H^{2}}{\alpha GM}+\dfrac{2}{D^{3}}\right.$", EQ),
 (QQ + r"$-\dfrac{3}{" + Sx + r"}$", EQ),
 (QQ + r"$\mp\dfrac{3d}{D}\left(\dfrac{1}{D^{3}}\right.$", EQ),
 (QQQ + r"$\left.\left.-\dfrac{2L+(1-f_u)D}{2L\left[" + Sx + r"\right]}\right)\right\}$", EQ),
 (Q + r"$\div\left\{\dfrac{6}{D^{3}}-\dfrac{6}{" + Sx + r"}\right.$", EQ),
 (QQ + r"$\mp\dfrac{9d}{2D}\left(\dfrac{4}{D^{3}}\right.$", EQ),
 (QQQ + r"$\left.\left.-\dfrac{4L+(1-f_u)D}{L\left[" + Sx + r"\right]}\right)\right\}$", EQ),
 ("gap", 0.2),
 (r"$\bf{Angles\ of\ the\ maxima}$", HD),
 (r"$\ell_{\max,1}\approx\dfrac{180^\circ}{\pi}\,\delta_1,\qquad\ell_{\max,2}\approx180^\circ-\dfrac{180^\circ}{\pi}\,\delta_2$", EQ),
]
pages[4] = [
 (r"$\bf{Page\ 4:\ heights\ of\ the\ maxima}$", HD),
 (r"(upper sign max,1 at $\approx 0^\circ$; lower sign max,2 at $\approx 180^\circ$)", NT), ("gap", 0.1),
 (r"$\left(\dfrac{A}{A_{\rm all}}\right)_{\max,1,2}\approx\left(1\mp\mathrm{Amp}_A\sin\phi_A-\dfrac{1}{2}\mathrm{Amp}_A\cos\phi_A\,\delta_{1,2}\right)$", EQ),
 (Q + r"$\times\dfrac{1}{A_{\rm all}}\left\{1+\dfrac{\alpha GM}{H^{2}}\left[\dfrac{2}{D^{3}}\right.\right.$", EQ),
 (QQ + r"$-\dfrac{3}{" + Sx + r"}$", EQ),
 (QQ + r"$\mp\dfrac{3d}{D}\left(\dfrac{1}{D^{3}}\right.$", EQ),
 (QQQ + r"$\left.\left.\left.-\dfrac{2L+(1-f_u)D}{2L\left[" + Sx + r"\right]}\right)\right]\right\}$", EQ),
 ("gap", 0.2),
 (r"$\delta_1$ and $\delta_2$ are on page 3; $A_{\rm all}$ is on page 1.", NT),
]
pages[5] = [
 (r"$\bf{Page\ 5:\ the\ minima}$", HD),
 (r"(radians; upper sign $\delta'_1$ / min,1 near $90^\circ$, lower sign $\delta'_2$ / min,2 near $270^\circ$)", NT), ("gap", 0.1),
 (r"$\delta'_{1,2}\approx-\mathrm{Amp}_A\left(\sin\phi_A\pm\dfrac{3d}{4D}\cos\phi_A\right)\left(\dfrac{H^{2}}{\alpha GM}-\dfrac{1}{D^{3}}\right)$", EQ),
 (Q + r"$\div\left(\dfrac{6}{D^{3}}-\dfrac{6}{" + Sx + r"}\right)$", EQ),
 ("gap", 0.2),
 (r"$\bf{Angles\ of\ the\ minima}$", HD),
 (r"$\ell_{\min,1}\approx90^\circ+\dfrac{135^\circ}{\pi}\,\dfrac{d}{D}+\dfrac{180^\circ}{\pi}\,\delta'_1$", EQ),
 (r"$\ell_{\min,2}\approx270^\circ-\dfrac{135^\circ}{\pi}\,\dfrac{d}{D}-\dfrac{180^\circ}{\pi}\,\delta'_2$", EQ),
 ("gap", 0.2),
 (r"$\bf{Heights\ of\ the\ minima}$", HD),
 (r"$\left(\dfrac{A}{A_{\rm all}}\right)_{\min,1,2}\approx\left(1\mp\mathrm{Amp}_A\cos\phi_A+\dfrac{3d}{4D}\mathrm{Amp}_A\sin\phi_A\right.$", EQ),
 (QQ + r"$\left.+\dfrac{1}{2}\mathrm{Amp}_A\sin\phi_A\,\delta'_{1,2}\right)\times\dfrac{1}{A_{\rm all}}\left(1-\dfrac{\alpha GM}{D^{3}H^{2}}\right)$", EQ),
]
pages[6] = [
 (r"$\bf{Page\ 6:\ notes}$", HD), ("gap", 0.1),
 (B + r"$\delta_1,\delta_2,\delta'_1,\delta'_2$ are in radians; $180^\circ/\pi$ turns them into degrees.", NT),
 (B + r"Amp$_A\cos\phi_A=-0.0564$,   Amp$_A\sin\phi_A=-0.0206$.", NT),
 (B + r"First order in $d/D$ and in Amp$_A$, keeping the Amp$_A\times d/D$ cross terms:", NT),
 (r"$\qquad$the $d$ parts of $\delta_{1,2}$ and $\delta'_{1,2}$, and $\frac{3d}{4D}$Amp$_A\sin\phi_A$ in the minimum heights.", NT),
 (B + r"Those cross terms bring $\ell_{\max,2}$ and both minima to within 0.1$^\circ$ at baseline.", NT),
 (B + r"$\ell_{\max,1}$ stays 1.1$^\circ$ high, from $(d/D)^2$ terms this order drops.", NT),
 (B + r"The $\frac{1}{2}$Amp$_A\,\delta$ terms in the heights: at the shifted peak the tilt factor", NT),
 (r"$\qquad$rises by Amp$_A\times$shift and $H_{\rm obs}/H$ falls by half of that.", NT),
 (B + r"They cut the height error at the maxima from 0.006 to 0.001 at baseline.", NT),
 (B + r"Sky mean: the tilt's $\sin\ell$ part averages to zero. Its product with the $\cos\ell$ part", NT),
 (r"$\qquad$of $H_{\rm obs}/H$ changes the mean of the tilted curve by", NT),
 (Q + r"$-\dfrac{9}{16}\,\mathrm{Amp}_A\sin\phi_A\,\dfrac{\alpha GM\,d}{DH^{2}}\left[\dfrac{L+(1-f_u)D}{L\left[" + Sx + r"\right]}\right.$", 16),
 (QQQ + r"$\left.-\dfrac{1}{D^{3}}\right]$", 16),
 (r"$\qquad$This is $%.1f\times10^{-5}$ at baseline (second order), and the model divides" % ((meanT - Aall_f)*1e5), NT),
 (r"$\qquad$by the mean of $H_{\rm obs}/H$, so $A_{\rm all}$ on page 1 is unchanged.", NT),
 (B + r"$\alpha$ only appears as $\alpha GM$, so the curve sees $\alpha M$, not $\alpha$ and $M$ separately.", NT),
]
TBL = ("table", None)
pages[7] = [
 (r"$\bf{Page\ 7:\ baseline\ check\ and\ validity}$", HD),
 (r"Baseline: $D/L$ = 5.9, $d/D$ = 0.1, $\alpha GM/D^{3}H^{2}$ = 0.2306", NT), ("gap", 0.1),
 TBL, ("gap", 0.1),
 (r"$A_{\rm all}$: formula %.4f, exact %.4f.   Exact = the no-boxcar curve with the tilt." % (Aall_f, hbar), NT),
 ("gap", 0.25),
 (r"$\bf{Validity}$", HD),
 (B + r"$D/L$: max,1 within 2$^\circ$ for $D/L\geq4.5$ and within 0.7$^\circ$ for $D/L\geq8$;", NT),
 (r"$\qquad$the other three within 0.7$^\circ$ for $D/L\geq4.3$; heights within 0.008.", NT),
 (B + r"$\alpha$ and $M$: max,1 within 2$^\circ$ for $\alpha GM/D^{3}H^{2}\geq0.16$, the others within 0.4$^\circ$,", NT),
 (r"$\qquad$heights within 0.012; fails for $\alpha GM>D^{3}H^{2}$ ($H_{\rm obs}<0$).", NT),
 (B + r"$d/D$ (at $\alpha$=1.3): max,1 within 2$^\circ$ up to $d/D\approx0.17$, fails near 0.3", NT),
 (r"$\qquad$(the $0^\circ$ peak flattens); max,2 within 1$^\circ$ up to 0.4, minima up to 0.5.", NT),
 (B + r"Below $D/L\approx2.36$ or $\alpha GM/D^{3}H^{2}\approx0.044$ the tilt wins: max,2 and min,1", NT),
 (r"$\qquad$merge and vanish, and the formulas no longer apply.", NT),
]

WMAX = 9.4
def render(n, items, dpi=250):
    fig = plt.figure(figsize=(10, 30)); r = fig.canvas.get_renderer()
    y = 30 - 0.35; x0 = 0.35; widest = 0
    for it in items:
        if it[0] == "gap": y -= it[1]; continue
        if it[0] == "table":
            rows = [("", r"$\ell$ formula", r"$\ell$ exact", r"$A/A_{\rm all}$ formula", r"$A/A_{\rm all}$ exact")]
            names = ["max,1", "max,2", "min,1", "min,2"]
            for i in range(4):
                rows.append((names[i], "%.2f$^\\circ$" % lf[i], "%.2f$^\\circ$" % le[i], "%.4f" % Af[i], "%.4f" % Ae[i]))
            cols = [x0, x0 + 1.2, x0 + 3.0, x0 + 4.8, x0 + 7.2]
            for j, row in enumerate(rows):
                for c, txt in zip(cols, row):
                    fig.text(c/10, y/30, txt, fontsize=16, va="top", fontweight="bold" if j == 0 else "normal")
                y -= 0.42 if j else 0.5
                if j == 0:
                    fig.add_artist(plt.Line2D([x0/10, 9.6/10], [(y + 0.06)/30]*2, color="k", lw=0.8))
            continue
        t, s = it
        ind = 0.55*(len(t) - len(t.lstrip("\x01"))); t = t.lstrip("\x01")
        tx = fig.text((x0 + ind)/10, y/30, t, fontsize=s, va="top")
        bb = tx.get_window_extent(r); w = bb.width/fig.dpi + ind; h = bb.height/fig.dpi
        widest = max(widest, w)
        if w > WMAX: print(f"  p{n} TOO WIDE {w:.2f}: {t[:70]}")
        y -= h + (0.28 if s >= EQ else 0.14)
    used = 30 - y + 0.25
    fig.set_size_inches(10, used, forward=False)
    for a in fig.texts: a.set_y(1 - (1 - a.get_position()[1])*30/used)
    for l in fig.artists if hasattr(fig, "artists") else []:
        yd = l.get_ydata(); l.set_ydata([1 - (1 - v)*30/used for v in yd])
    out = f"extrema-nobox-tilt-p{n}.png"
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    print(f"p{n}: 10 x {used:.2f} in, widest line {widest:.2f} in")
for n, items in pages.items():
    render(n, items)
