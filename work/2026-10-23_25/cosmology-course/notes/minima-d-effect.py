import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"mathtext.fontset":"cm","font.family":"serif"})
B  = r"\left[(1-f_u)(6L^3e^{D/L}-6L^3-6DL^2-3D^2L)+f_uD^3\right]"
Bm = r"(1-f_u)(6L^3e^{D/L}-6L^3-6DL^2-3D^2L-D^3)"
Q  = r"\dfrac{(D^3H^2-\alpha GM)\,%s}{6\,\alpha GM\,%s}" % (B, Bm)
T = [
 (r"$\bf{Effect\ of\ observation\ distance\ }d\bf{\ on\ the\ minima\ (no\ boxcar,\ with\ sine\ tilt)}$", 17),
 ("gap",0.3),
 (r"$\bf{Positions}$",17),
 (r"$\ell_{\min,1}\approx90^\circ+\dfrac{180^\circ}{\pi}\left[\dfrac{3d}{4D}-\mathrm{Amp}_A\left(\sin\phi_A+\dfrac{3d}{4D}\cos\phi_A\right)\right.$",18),
 (r"$\qquad\qquad\left.\times\,%s\right]$" % Q,18),
 ("gap",0.2),
 (r"$\ell_{\min,2}\approx270^\circ-\dfrac{180^\circ}{\pi}\left[\dfrac{3d}{4D}-\mathrm{Amp}_A\left(\sin\phi_A-\dfrac{3d}{4D}\cos\phi_A\right)\right.$",18),
 (r"$\qquad\qquad\left.\times\,%s\right]$" % Q,18),
 ("gap",0.35),
 (r"$\bf{Rate\ of\ change\ with\ }d$",17),
 (r"$\dfrac{\partial\ell_{\min,1}}{\partial d}\approx\dfrac{135^\circ}{\pi D}\left[1-\mathrm{Amp}_A\cos\phi_A\right.$",18),
 (r"$\qquad\qquad\left.\times\,%s\right]$" % Q,18),
 ("gap",0.2),
 (r"$\dfrac{\partial\ell_{\min,2}}{\partial d}\approx-\dfrac{135^\circ}{\pi D}\left[1+\mathrm{Amp}_A\cos\phi_A\right.$",18),
 (r"$\qquad\qquad\left.\times\,%s\right]$" % Q,18),
 ("gap",0.35),
 (r"$\bullet$  The $135^\circ/(\pi D)$ term is geometry: $\cos\ell_{\min}=-3d/(4D)$, so both minima swing toward the centre side by $\approx43^\circ$ per unit $d/D$.",15),
 (r"$\bullet$  The bracket is the sine tilt's correction; it does not depend on $d$.",15),
 (r"$\bullet$  $f_u=0.94$, Amp$_A=-0.06$, $\phi_A=0.35$ rad; $\alpha$ enters only as $\alpha GM$.",15),
 (r"$\bullet$  Baseline ($D/L=5.9$, $d/D=0.1$, $\alpha GM/D^3H^2=0.2306$):  $+46.8^\circ$ and $-39.1^\circ$ per unit $d/D$, within 0.5% of the exact curve.",15),
 (r"$\bullet$  First order in $d/D$: within $0.1^\circ$ of exact up to $d/D\approx0.2$, $0.9^\circ$ at $d/D=0.5$, $9.5^\circ$ at $d/D=0.9$.",15),
]
H={17:0.5,18:1.05,15:0.42}
tot=sum(v if t=="gap" else H[v] for t,v in T)
fig=plt.figure(figsize=(14,tot+0.4)); y=tot+0.2
for t,v in T:
    if t=="gap": y-=v; continue
    fig.text(0.02,y/(tot+0.4),t,fontsize=v,va="top"); y-=H[v]
fig.savefig("minima-d-effect.png",dpi=220,bbox_inches="tight",pad_inches=0.3)
