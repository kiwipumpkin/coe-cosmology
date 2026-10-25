"""Extremum angles and heights of A/A_all written in the physical inputs
D, L, d, M (with G, H, alpha). First order in d/D and Amp_A (eqs 6.4-6.8 of
extrema-equations.md, kappa=GM/(D^3H^2), x=D/L, eps=d/D substituted and
collected). Nothing fitted."""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"mathtext.fontset": "cm", "font.family": "serif"})
fu = 0.94; W = np.radians(90.75)
s1, s2, s3 = [np.sin(j*W/2)/(j*W/2) for j in (1, 2, 3)]
def ext(D, L, d, M, al=1.0, G=1.0, H=1.0, Amp=-0.06, ph=0.35, sgn=1):
    B = (1-fu)*(6*L**3*np.exp(D/L) - 6*L**3 - 6*D*L**2 - 3*D**2*L - D**3)
    gm = al*G*M; base = 2*D**3*(D**3*H**2 - gm)
    P = (2*D**3*H**2 + gm)*B + base
    def numax(dd):
        return (2*D**3*H**2 + gm*(1 + 3*s2 - 3*dd/(4*D)*(3*s1 + 5*s3)))*B + base \
               + 3*gm*dd*D**3/(4*L)*(1-fu)*(3*s1 + s3)
    def demax(dd):
        return gm*((12*s2 - 9*dd/(4*D)*(s1 + 15*s3))*B + 9*dd*D**3/(4*L)*(1-fu)*(s1 + 3*s3))
    numin = (2*D**3*H**2 + gm*(1 - 3*s2))*B + base
    dl = -3*d/(16*s2*D)*((s1 - 5*s3) - (1-fu)*(s1 - s3)*D**4/(L*B))
    l1 = -np.degrees(Amp*np.cos(ph)*numax(d)/demax(d))
    l3 = 180 + np.degrees(Amp*np.cos(ph)*numax(-d)/demax(-d))
    tm = np.degrees(Amp*np.sin(ph)*numin/(12*s2*gm*B))
    l2 = 90 + np.degrees(dl) - tm; l4 = 270 - np.degrees(dl) + tm
    A1 = (1-Amp*np.sin(ph))*numax(d)/P; A3 = (1+Amp*np.sin(ph))*numax(-d)/P
    A2 = (1-Amp*np.cos(ph))*numin/P; A4 = (1+Amp*np.cos(ph))*numin/P
    return (l1, l3, l2, l4), (A1, A3, A2, A4)
L_, A_ = ext(1.0, 1/5.9, 0.1, 0.230588)
print("s", s1, s2, s3, "\nl", L_, "\nA", A_)

Bx = r"(1-f_u)\left[6L^{3}e^{D/L}-6L^{3}-6DL^{2}-3D^{2}L-D^{3}\right]"
T = [
 (r"$\bf{Extremum\ angles\ and\ heights\ of\ }A/A_{\rm all}\bf{\ in\ the\ physical\ inputs\ }D,\,L,\,d,\,M$", 16),
 (r"$D$: CoE–observer distance;  $L$: core scale length of the density;  $d$: radius of the source shell;  $M$: mass inside $D$;  $G$: Newton's constant;  $H$: background expansion rate.", 12.5),
 ("gap", 0.25),
 (r"$\bf{Maximum\ near\ 10^\circ}$", 14),
 (r"$\ell_{\max,1}\approx-\dfrac{180^\circ}{\pi}\,\mathrm{Amp}_A\cos\phi_A\times$", 15),
 (r"$\dfrac{\left[2D^{3}H^{2}+\alpha GM\left(1+3s_2-\dfrac{3d}{4D}(3s_1+5s_3)\right)\right]" + Bx + r"+2D^{3}\left(D^{3}H^{2}-\alpha GM\right)+\dfrac{3\,\alpha GM\,d\,D^{3}}{4L}(1-f_u)(3s_1+s_3)}"
  r"{\alpha GM\left\{\left[12s_2-\dfrac{9d}{4D}(s_1+15s_3)\right]" + Bx + r"+\dfrac{9\,d\,D^{3}}{4L}(1-f_u)(s_1+3s_3)\right\}}$", 13.5),
 (r"$\left(\dfrac{A}{A_{\rm all}}\right)_{\max,1}\approx(1-\mathrm{Amp}_A\sin\phi_A)\times$", 15),
 (r"$\dfrac{\left[2D^{3}H^{2}+\alpha GM\left(1+3s_2-\dfrac{3d}{4D}(3s_1+5s_3)\right)\right]" + Bx + r"+2D^{3}\left(D^{3}H^{2}-\alpha GM\right)+\dfrac{3\,\alpha GM\,d\,D^{3}}{4L}(1-f_u)(3s_1+s_3)}"
  r"{\left(2D^{3}H^{2}+\alpha GM\right)" + Bx + r"+2D^{3}\left(D^{3}H^{2}-\alpha GM\right)}$", 13.5),
 (r"$\bf{Maximum\ near\ 170^\circ}$:  $\ell_{\max,2}\approx180^\circ+\dfrac{180^\circ}{\pi}\,\mathrm{Amp}_A\cos\phi_A\times\left[\mathrm{the\ }\ell_{\max,1}\mathrm{\ fraction\ with\ }d\to-d\right]$,", 14),
 (r"$\left(A/A_{\rm all}\right)_{\max,2}\approx(1+\mathrm{Amp}_A\sin\phi_A)\times\left[\mathrm{the\ }(A/A_{\rm all})_{\max,1}\mathrm{\ fraction\ with\ }d\to-d\right]$", 14),
 ("gap", 0.25),
 (r"$\bf{Minima\ near\ 95^\circ\ and\ 266^\circ}$  (upper sign $\ell_{\min,1}$, lower sign $\ell_{\min,2}$)", 14),
 (r"$\ell_{\min}\approx\left\{\begin{array}{c}90^\circ\\270^\circ\end{array}\right\}$", 15) if False else
 (r"$\ell_{\min,1,2}\approx\left(90^\circ\ \mathrm{or}\ 270^\circ\right)\mp\dfrac{180^\circ}{\pi}\,\dfrac{3d}{16\,s_2\,D}\left[s_1-5s_3-\dfrac{(1-f_u)(s_1-s_3)\,D^{4}}{L\," + Bx + r"}\right]$", 15),
 (r"$\qquad\qquad\mp\dfrac{180^\circ}{\pi}\,\mathrm{Amp}_A\sin\phi_A\,\dfrac{\left[2D^{3}H^{2}+\alpha GM(1-3s_2)\right]" + Bx + r"+2D^{3}\left(D^{3}H^{2}-\alpha GM\right)}{12\,s_2\,\alpha GM\;" + Bx + r"}$", 15),
 (r"$\left(\dfrac{A}{A_{\rm all}}\right)_{\min,1,2}\approx(1\mp\mathrm{Amp}_A\cos\phi_A)\,\dfrac{\left[2D^{3}H^{2}+\alpha GM(1-3s_2)\right]" + Bx + r"+2D^{3}\left(D^{3}H^{2}-\alpha GM\right)}{\left(2D^{3}H^{2}+\alpha GM\right)" + Bx + r"+2D^{3}\left(D^{3}H^{2}-\alpha GM\right)}$", 15),
 ("gap", 0.3),
 (r"Constants: $s_j=\sin(jW/2)\,/\,(jW/2)$ from the $W=90.75^\circ$ boxcar: $s_1=%.4f$, $s_2=%.4f$, $s_3=%.4f$;  $f_u=0.94$ (floor share of the density);  Amp$_A=-0.06$, $\phi_A=0.35$ rad." % (s1, s2, s3), 12.5),
 (r"$\alpha$ (peculiar-velocity scale) only ever multiplies $GM$, so the curve sees $\alpha M$, not $\alpha$ and $M$ separately.", 12.5),
 (r"The bracket $" + Bx + r"=6(1-f_u)L^{3}\left[e^{D/L}-1-\dfrac{D}{L}-\dfrac{D^{2}}{2L^{2}}-\dfrac{D^{3}}{6L^{3}}\right]$", 12.4),
 (r"is the exponential tail of the mass profile beyond its cubic part (it is $3M/4\pi\rho_D-D^{3}$, with $\rho_D$ the density at the observer).", 12.5),
 (r"Baseline (D/L = 5.9, d/D = 0.1, $\alpha GM/D^{3}H^{2}$ = 0.2306):  $\ell$ = %.1f$^\circ$, %.1f$^\circ$, %.1f$^\circ$, %.1f$^\circ$ (exact 10.4, 170.1, 95.3, 265.7);  $A/A_{\rm all}$ = %.3f, %.3f, %.3f, %.3f (exact 1.106, 1.076, 0.964, 0.862)." % (*L_, *A_), 12.5),
 (r"First order in $d/D$ and Amp$_A$ (it drops terms of order $(d/D)^2$); good near baseline, breaks down close to the merges (D/L $\approx$ 3.16, $\alpha GM/D^{3}H^{2}\approx$ 0.0705) and for $d/D$ near 1.", 12.5),
]
H = {16: 0.6, 15: 1.15, 14: 0.75, 13.5: 1.3, 12.5: 0.42, 12.4: 0.62}
tot = sum(v if t == "gap" else H[v] for t, v in T)
fig = plt.figure(figsize=(17, tot + 0.4)); y = tot + 0.2
for t, v in T:
    if t == "gap": y -= v; continue
    fig.text(0.01, y / (tot + 0.4), t, fontsize=v, va="top"); y -= H[v]
fig.savefig("extrema-physical.png", dpi=200, bbox_inches="tight", pad_inches=0.25)
