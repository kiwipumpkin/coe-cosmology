"""Extremum angles and heights of A/A_all(l) written directly in the inputs
(lowest order in d/D and Amp_A; eqs 6.7-6.8 of extrema-equations.md with a2
substituted and the fractions combined). Nothing fitted."""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"mathtext.fontset": "cm", "font.family": "serif"})
fu, Wd = 0.94, 90.75; W = np.radians(Wd); s = np.sin(W)/W
def ext(x, ak, Amp=-0.06, ph=0.35):
    N = 6*(1-fu)*(np.exp(x)-1-x-x**2/2) + fu*x**3
    den = 3*s*ak*(N-x**3); P = (2+ak)*N - 3*ak*x**3
    Fp = (2+ak+3*s*ak)*N - 3*ak*(1+s)*x**3
    Fm = (2+ak-3*s*ak)*N - 3*ak*(1-s)*x**3
    l1 = -(45/np.pi)*Amp*np.cos(ph)*Fp/den
    l2 = 90 - (45/np.pi)*Amp*np.sin(ph)*Fm/den
    A1 = (1-Amp*np.sin(ph))*Fp/P; A3 = (1+Amp*np.sin(ph))*Fp/P
    A2 = (1-Amp*np.cos(ph))*Fm/P; A4 = (1+Amp*np.cos(ph))*Fm/P
    return (l1, l2, 180-l1, 360-l2), (A1, A2, A3, A4)
L, A = ext(5.9, 0.230588)
print("s", s, "l", L, "A", A)
T = [
 (r"$\bf{Extremum\ angles\ and\ heights\ of\ }A/A_{\rm all}\bf{\ as\ functions\ of\ the\ inputs}$", 16),
 ("gap", 0.35),
 (r"$N(x)=6\,(1-f_u)\left(e^{x}-1-x-\dfrac{x^{2}}{2}\right)+f_u\,x^{3},\qquad s=\dfrac{\sin W}{W}=%.4f$" % s, 15),
 ("gap", 0.25),
 (r"$\bf{Maxima}$", 14),
 (r"$\ell_{\max,1}\approx-\dfrac{45}{\pi}\,\mathrm{Amp}_A\cos\phi_A\;\dfrac{(2+\alpha\kappa+3s\,\alpha\kappa)\,N(x)-3\,\alpha\kappa\,(1+s)\,x^{3}}{3s\,\alpha\kappa\,\left[N(x)-x^{3}\right]},\qquad\ell_{\max,2}\approx180^\circ-\ell_{\max,1}$", 15),
 ("gap", 0.15),
 (r"$(A/A_{\rm all})_{\max,1}\approx(1-\mathrm{Amp}_A\sin\phi_A)\;\dfrac{(2+\alpha\kappa+3s\,\alpha\kappa)\,N(x)-3\,\alpha\kappa\,(1+s)\,x^{3}}{(2+\alpha\kappa)\,N(x)-3\,\alpha\kappa\,x^{3}}$", 15),
 (r"$(A/A_{\rm all})_{\max,2}\approx(1+\mathrm{Amp}_A\sin\phi_A)\;\dfrac{(2+\alpha\kappa+3s\,\alpha\kappa)\,N(x)-3\,\alpha\kappa\,(1+s)\,x^{3}}{(2+\alpha\kappa)\,N(x)-3\,\alpha\kappa\,x^{3}}$", 15),
 ("gap", 0.25),
 (r"$\bf{Minima}$", 14),
 (r"$\ell_{\min,1}\approx90^\circ-\dfrac{45}{\pi}\,\mathrm{Amp}_A\sin\phi_A\;\dfrac{(2+\alpha\kappa-3s\,\alpha\kappa)\,N(x)-3\,\alpha\kappa\,(1-s)\,x^{3}}{3s\,\alpha\kappa\,\left[N(x)-x^{3}\right]},\qquad\ell_{\min,2}\approx360^\circ-\ell_{\min,1}$", 15),
 ("gap", 0.15),
 (r"$(A/A_{\rm all})_{\min,1}\approx(1-\mathrm{Amp}_A\cos\phi_A)\;\dfrac{(2+\alpha\kappa-3s\,\alpha\kappa)\,N(x)-3\,\alpha\kappa\,(1-s)\,x^{3}}{(2+\alpha\kappa)\,N(x)-3\,\alpha\kappa\,x^{3}}$", 15),
 (r"$(A/A_{\rm all})_{\min,2}\approx(1+\mathrm{Amp}_A\cos\phi_A)\;\dfrac{(2+\alpha\kappa-3s\,\alpha\kappa)\,N(x)-3\,\alpha\kappa\,(1-s)\,x^{3}}{(2+\alpha\kappa)\,N(x)-3\,\alpha\kappa\,x^{3}}$", 15),
 ("gap", 0.35),
 (r"$\bf{Notes}$", 14),
 (r"$\bullet$  $x=D/L$: distance from the centre of expansion (CoE) to the observer, divided by the core scale $L$.", 12.5),
 (r"$\bullet$  $\alpha$ is the velocity scale and $\kappa$ the mass; they enter only as the product $\alpha\kappa$.", 12.5),
 (r"$\bullet$  $f_u=0.94$ is the floor share of the density;  $W=%.2f^\circ$ is the boxcar width;  Amp$_A=-0.06$ and $\phi_A=0.35$ rad set the sine tilt." % Wd, 12.5),
 (r"$\bullet$  Which is which: max,1 $\approx10^\circ$, max,2 $\approx170^\circ$, min,1 $\approx95^\circ$, min,2 $\approx266^\circ$.", 12.5),
 (r"$\bullet$  Baseline ($x$=5.9, $\alpha\kappa$=0.2306) angles:  $\ell$ = %.1f$^\circ$, %.1f$^\circ$, %.1f$^\circ$, %.1f$^\circ$  (exact 10.4, 170.1, 95.3, 265.7)." % (L[0], L[2], L[1], L[3]), 12.5),
 (r"$\bullet$  Baseline heights:  $A/A_{\rm all}$ = %.3f, %.3f, %.3f, %.3f  (exact 1.106, 1.076, 0.964, 0.862)." % (A[0], A[2], A[1], A[3]), 12.5),
 (r"$\bullet$  Lowest order in $d/D$ and Amp$_A$: $d/D$ drops out here and enters at the next order (its sensitivities are on the equation card).", 12.5),
 (r"$\bullet$  Angles stay within 2$^\circ$ for D/L $\geq$ 4.9, $\alpha$ $\geq$ 0.73, $\kappa$ 0.17–0.70, d/D $\leq$ 0.32.", 12.5),
 (r"$\bullet$  Breaks down near the merges (D/L$\approx$3.16, $\alpha\kappa\approx$0.0705), where the second max and min vanish.", 12.5),
]
H = {16: 0.55, 15: 0.9, 14: 0.42, 12.5: 0.36}
tot = sum(t[1] if t[0] == "gap" else H[t[1]] for t in T)
fig = plt.figure(figsize=(13.5, tot + 0.4)); y = tot + 0.2
for t, v in T:
    if t == "gap": y -= v; continue
    fig.text(0.02, y / (tot + 0.4), t, fontsize=v, va="top"); y -= H[v]
fig.savefig("extrema-closedform.png", dpi=200, bbox_inches="tight", pad_inches=0.25)
