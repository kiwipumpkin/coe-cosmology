"""Pages: d(ell)/dd and the total differential of the four tilt-card angles, written explicitly
(no delta names). Forms checked against sympy in nobox_dl_dd_forms (to 3e-14)."""
import json, numpy as np
import pagekit as pk
from pagekit import render, note, EQ, NT, HD, TT, B
pk.PREFIX = "extrema-dl-dd"
J = json.load(open("nobox_tilt_deriv.json"))
from nobox_dl_dd_curves import formulas
Z = np.load("nobox_dl_dd_curves.npz"); ES, EX = Z["es"], Z["ex"]; DEX = np.gradient(EX, ES, axis=0)
Sx = r"(1-f_u)\left(6L^{3}e^{D/L}-6L^{3}-6DL^{2}-3D^{2}L\right)+f_uD^{3}"
SDx = r"(1-f_u)\left(6L^{2}e^{D/L}-6L^{2}-6DL\right)+3f_uD^{2}"
SLx = r"(1-f_u)\left(18L^{2}e^{D/L}-6DLe^{D/L}-18L^{2}-12DL-3D^{2}\right)"
BS, BSD, BSL = r"\left[" + Sx + r"\right]", r"\left[" + SDx + r"\right]", r"\left[" + SLx + r"\right]"
BT = r"\left[(1-f_u)\left(6L^{3}e^{D/L}-6L^{3}-6DL^{2}-3D^{2}L-D^{3}\right)\right]"
deg = r"\dfrac{180^\circ}{\pi}"
def QB(suffix=""):
    return [[r"@L[4LD" + BT, None], [r"\mp3d@L(4L" + BT, None], [r"-(1-f_u)D^{4}@R)@R]" + suffix, None]]
def div_Q(pre, suffix):   # "÷ { pre [Q] suffix }"
    q = QB(suffix + r"@R\}"); q[0][0] = r"\div@L\{" + pre + q[0][0]; q[0][1] = "!"; return q
def f4(v): return ("%#.4g" % v).rstrip(".").replace("-", "\u2212")
def ae(nm, k): return r"analytic $%s$, exact $%s$" % (f4(J["ana"][nm][k]).replace("\u2212", "-"), f4(J["efd"][k][nm]).replace("\u2212", "-"))
def DMIN(first, pre=r"\pm\dfrac{135^\circ}{\pi D}"):
    return [[pre + r"@L[" + first, None], [r"\mp\dfrac{\mathrm{Amp}_A\cos\phi_A\left(D^{3}H^{2}-\alpha GM\right)}{6\,\alpha GM}", None],
            [r"\times\dfrac{" + BS + "}{" + BT + r"}@R]", None]]
title = r"$\bf{How\ the\ angles\ move\ with\ }d\bf{,\ and\ the\ full\ differential\ (no\ boxcar,\ with\ tilt)}$"
P = {}
TILT = r"\dfrac{\mathrm{Amp}_A\left(4D\sin\phi_A\pm3d\cos\phi_A\right)\left(D^{3}H^{2}-\alpha GM\right)}{24D\,\alpha GM}"
P[1] = [
 (title, 19), (r"$\bf{Page\ 1:\ the\ four\ angles,\ written\ out}$", HD), ("gap", 0.05),
 *note(r"$\ell$ = angle of an extremum; $d$ = observation distance (source-shell radius); $D$ = distance to the centre of expansion (CoE); $L$ = core scale length; $M$ = mass inside $D$. $\ell=0^\circ$ points away from the CoE, $180^\circ$ towards it. Upper sign: index 1, lower sign: index 2."),
 (r"$\bf{Minima}$  ($\ell_{\min,1}-90^\circ$ and $270^\circ-\ell_{\min,2}$):", NT),
 ("eq", r"=", [[deg + r"@L[\dfrac{3d}{4D}", None], [r"-" + TILT, None], [r"\times\dfrac{" + BS + "}{" + BT + r"}@R]", None]]),
 (r"$\bf{Maxima}$  ($\ell_{\max,1}$ and $180^\circ-\ell_{\max,2}$):", NT),
 ("eq", r"=", [[r"-" + deg + r"\mathrm{Amp}_A\cos\phi_A", None], [r"@L\{2LD^{4}H^{2}" + BS, "!"],
               [r"+2\alpha GM\,LD@L(2" + BS + r"-3D^{3}@R)", "!"], [r"\mp3d\,\alpha GM@L(2L" + BT, "!"], [r"-(1-f_u)D^{4}@R)@R\}", None]]
              + div_Q(r"3\alpha GM", "")),
 *note(r"The second bracket, $(1-f_u)(6L^{3}e^{D/L}-6L^{3}-6DL^{2}-3D^{2}L-D^{3})$, is the first minus $D^{3}$; it is always positive."),
 *note(r"Baseline (code units): $D=1$, $L=D/5.9$, $d=0.1$, $M=0.2306$, $\alpha=G=H=1$, $f_u=0.94$, Amp$_A=-0.06$, $\phi_A=0.35$ rad. Derivatives are in degrees per unit of the input."),
]
corr = -0.06*np.cos(0.35)*(1 - 0.230588)*1.5372553/(6*0.230588*0.5372553)
P[2] = [
 (r"$\bf{Page\ 2:\ }\partial\ell_{\min}/\partial d\bf{\ (what\ Chris\ asked\ for)}$", HD), ("gap", 0.05),
 ("eq", r"\dfrac{\partial\ell_{\min,1,2}}{\partial d}=", DMIN("1")),
 *note(r"Leading term $\pm135^\circ/(\pi D)$ ($=\pm42.97^\circ$ per unit $d$ at $D=1$): with no tilt, $\ell_{\min,1}=90^\circ+(180^\circ/\pi)(3d/4D)$, i.e. to first order $\cos\ell_{\min}=-3d/4D$."),
 *note(r"Meaning: for a tiny shell the minima sit at $90^\circ$ and $270^\circ$, perpendicular to the observer–CoE line. As $d$ grows they swing towards the CoE side ($180^\circ$), so the direction's component along the line is $-3d/4D$: both minima close in on $180^\circ$ by the same amount."),
 *note(r"The tilt factor in brackets does not depend on $d$; at baseline it is $1\pm%.4f$: $\ell_{\min,1}$ moves %.0f%% faster, $\ell_{\min,2}$ %.0f%% slower." % (-corr, -100*corr, -100*corr)),
 *note(r"Baseline: $\partial\ell_{\min,1}/\partial d$ " + ae("lmin1", "d") + r"; $\partial\ell_{\min,2}/\partial d$ " + ae("lmin2", "d") + "."),
 ("gap", 0.05), (r"$\bf{The\ arccos\ form}$ (exact inverse of $\cos\ell_{\min}=-3d/4D$, same tilt term)", HD),
 ("eq", r"\ell_{\min,1}=", [[r"\arccos@L(-\dfrac{3d}{4D}@R)", None], [r"-" + deg + TILT, None], [r"\times\dfrac{" + BS + "}{" + BT + "}", None]]),
 *note(r"$\ell_{\min,2}=360^\circ-\arccos(-3d/4D)+(180^\circ/\pi)\times$(the same tilt term, lower sign).", bullet=False),
 ("eq", r"\dfrac{\partial\ell_{\min,1,2}}{\partial d}=", DMIN(r"\dfrac{1}{\sqrt{1-9d^{2}/16D^{2}}}")),
 *note(r"Against the exact curve (page 11, figure): no difference up to $d/D\approx0.2$ (both within $0.08^\circ$); beyond, the arccos form roughly halves the angle error ($0.33^\circ$ vs $0.87^\circ$ at $d/D=0.5$, $5.7^\circ$ vs $9.5^\circ$ at $0.9$) and tracks the slope better, but the exact slope still grows faster (104 vs 62 vs 47 deg at $d/D=0.9$)."),
]
PdCh = [[deg + r"\,2\mathrm{Amp}_A\cos\phi_A\,DL", None], [r"@L\{D^{3}H^{2}" + BS, "!"], [r"@L((1-f_u)D^{4}-4L" + BT + r"@R)", r"\times"],
        [r"-\alpha GM@L((1-f_u)D^{7}", "!"], [r"+4L" + BT, None], [r"@L(" + BS + r"-2D^{3}@R)@R)@R\}", r"\times"]]
P[3] = [
 (r"$\bf{Page\ 3:\ }\partial\ell_{\max}/\partial d$", HD), ("gap", 0.05),
 ("eq", r"\dfrac{\partial\ell_{\max,1,2}}{\partial d}=", PdCh + div_Q(r"\alpha GM", "^{2}")),
 *note(r"Same sign for both maxima: both move the same way in $\ell$ as $d$ grows ($\ell_{\max,1}$ away from $0^\circ$, $\ell_{\max,2}$ towards $180^\circ$)."),
 *note(r"$d$ enters the maxima only through the Amp$_A\times d/D$ cross terms (the $\mp3d$ terms on page 1, both multiplied by Amp$_A\cos\phi_A$). With Amp$_A=0$ the maxima sit at exactly $0^\circ$ and $180^\circ$ for every $d$, and this derivative is zero: the whole effect is $\propto$ Amp$_A\cos\phi_A$."),
 *note(r"The squared bracket vanishes at $d/D=0.399$ for $\ell_{\max,1}$ (upper sign): the first-order $\ell_{\max,1}$ has a pole there; it is already 14% high at $d/D=0.1$ and unreliable beyond $\approx0.15$ (figure). $\ell_{\max,2}$ has no pole."),
 *note(r"Baseline: $\partial\ell_{\max,1}/\partial d$ " + ae("lmax1", "d") + r"; $\partial\ell_{\max,2}/\partial d$ " + ae("lmax2", "d") + "."),
]
P[4] = [
 (r"$\bf{Page\ 4:\ the\ total\ differential;\ minima}$", HD), ("gap", 0.05),
 (r"$d\ell=\dfrac{\partial\ell}{\partial D}dD+\dfrac{\partial\ell}{\partial L}dL+\dfrac{\partial\ell}{\partial d}dd+\dfrac{\partial\ell}{\partial M}dM+\dfrac{\partial\ell}{\partial\alpha}d\alpha,\qquad\dfrac{\partial\ell}{\partial\alpha}=\dfrac{M}{\alpha}\dfrac{\partial\ell}{\partial M}$", 17),
 (r"$\partial\ell_{\min}/\partial d$: page 2. The other three for the minima:", NT),
 ("eq", r"\dfrac{\partial\ell_{\min,1,2}}{\partial M}=", [[r"\pm" + deg + r"\dfrac{\mathrm{Amp}_A\left(4D\sin\phi_A\pm3d\cos\phi_A\right)D^{2}H^{2}}{24\,\alpha GM^{2}}", None], [r"\times\dfrac{" + BS + "}{" + BT + "}", None]]),
 ("eq", r"\dfrac{\partial\ell_{\min,1,2}}{\partial L}=", [[r"\pm" + deg + r"\dfrac{\mathrm{Amp}_A D^{2}\left(4D\sin\phi_A\pm3d\cos\phi_A\right)\left(D^{3}H^{2}-\alpha GM\right)}{24\,\alpha GM}", None], [r"\times\dfrac{" + BSL + "}{" + BT + r"^{2}}", None]]),
 ("eq", r"\dfrac{\partial\ell_{\min,1,2}}{\partial D}=", [[r"\pm" + deg + r"@L\{-\dfrac{3d}{4D^{2}}", None], [r"-\dfrac{\mathrm{Amp}_A}{24\,\alpha GM}@L[@L(16D^{3}H^{2}\sin\phi_A-4\alpha GM\sin\phi_A\pm9dD^{2}H^{2}\cos\phi_A@R)", "!"],
     [r"\times\dfrac{" + BS + "}{D" + BT + "}", None],
     [r"+\dfrac{\left(4D\sin\phi_A\pm3d\cos\phi_A\right)\left(D^{3}H^{2}-\alpha GM\right)}{D^{2}" + BT + r"^{2}}", "!"],
     [r"\times@L(4D^{3}" + BS, "!"], [r"-" + BS + r"^{2}", None], [r"-D^{4}" + BSD + r"@R)@R]@R\}", "!"]]),
 *note(r"Baseline $\ell_{\min,1}$: $\partial_M$ " + ae("lmin1", "M") + r"; $\partial_L$ " + ae("lmin1", "L") + r"; $\partial_D$ " + ae("lmin1", "D") + r". All partials, all four angles: page 10."),
]
PLch = [[r"@L\{-4DL^{2}@L(D^{3}H^{2}-\alpha GM@R)", None], [BSL, "!"], [r"-9(1-f_u)\alpha GM\,d^{2}@L(D^{3}", "!"], [r"-L" + BSL, None], [r"-" + BS + r"@R)", None],
        [r"\pm3d@L[(1-f_u)D^{4}@L(H^{2}@L(L" + BSL, "!"], [r"+" + BS + r"@R)-\alpha GM@R)", None],
        [r"+4L^{2}@L(D^{3}H^{2}-\alpha GM@R)", "!"], [BSL + r"@R]@R\}", None]]
P[5] = [
 (r"$\bf{Page\ 5:\ maxima,\ }\partial/\partial M\bf{\ and\ }\partial/\partial L$", HD), ("gap", 0.05),
 ("eq", r"\dfrac{\partial\ell_{\max,1,2}}{\partial M}=", [[r"\pm" + deg + r"\,2\mathrm{Amp}_A\cos\phi_A\,LD^{4}H^{2}", None], [BS, None]] + div_Q(r"3\alpha GM^{2}", "")),
 ("gap", 0.1),
 ("eq", r"\dfrac{\partial\ell_{\max,1,2}}{\partial L}=", [[r"\mp" + deg + r"\dfrac{2}{3}\mathrm{Amp}_A\cos\phi_A\,D^{4}", None]] + PLch + div_Q(r"\alpha GM", "^{2}")),
 *note(r"Baseline $\ell_{\max,1}$: $\partial_M$ " + ae("lmax1", "M") + r"; $\partial_L$ " + ae("lmax1", "L") + "."),
 *note(r"Baseline $\ell_{\max,2}$: $\partial_M$ " + ae("lmax2", "M") + r"; $\partial_L$ " + ae("lmax2", "L") + "."),
]
PDch = [[r"@L\{4D^{4}L@L[D" + BSD, "!"], [r"@L(D^{3}H^{2}-\alpha GM@R)", None], [r"-3" + BS, "!"], [r"@L(H^{2}" + BS + r"-\alpha GM@R)@R]", r"\times"],
        [r"-9(1-f_u)\alpha GM\,d^{2}@L(D^{6}+D^{4}" + BSD, "!"], [r"-4D^{3}" + BS + r"@R)", None],
        [r"\mp3d@L[(1-f_u)D^{8}H^{2}" + BSD, "!"], [r"+4D^{4}L" + BSD, "!"], [r"@L(D^{3}H^{2}-\alpha GM@R)", None],
        [r"+4D^{3}H^{2}L" + BS, "!"], [r"@L(D^{3}-4" + BS + r"@R)", r"\times"],
        [r"-4\alpha GM\,L@L(2D^{6}-6D^{3}" + BS, "!"], [r"+" + BS + r"^{2}@R)@R]@R\}", None]]
P[6] = [
 (r"$\bf{Page\ 6:\ maxima,\ }\partial/\partial D$", HD), ("gap", 0.05),
 ("eq", r"\dfrac{\partial\ell_{\max,1,2}}{\partial D}=", [[r"\pm" + deg + r"\dfrac{2}{3}\mathrm{Amp}_A\cos\phi_A\,L", None]] + PDch + div_Q(r"\alpha GM", "^{2}")),
 *note(r"Baseline: $\ell_{\max,1}$ " + ae("lmax1", "D") + r"; $\ell_{\max,2}$ " + ae("lmax2", "D") + "."),
 *note(r"Check: $D\,\partial_D\ell+L\,\partial_L\ell+d\,\partial_d\ell+3M\,\partial_M\ell=0$ holds for all four angles (to $10^{-14}$), so $\partial_D\ell$ also follows from the other three (page 7)."),
]
def cu(nm): return (-J["base"]["L"]*J["ana"][nm]["L"], -J["base"]["L"]*J["efd"]["L"][nm])
P[7] = [
 (r"$\bf{Page\ 7:\ the\ compact\ complete\ form;\ minima}$", HD), ("gap", 0.05),
 *note(r"Each angle depends only on $D/L$, $d/D$ and $\alpha GM/D^{3}H^{2}$, so", bullet=False),
 (r"$d\ell=C_{D/L}\,d\ln\dfrac{D}{L}+C_{d/D}\,d\!\left(\dfrac{d}{D}\right)+C_{M}\,d\ln\dfrac{\alpha GM}{D^{3}H^{2}}$", EQ),
 (r"$C_{D/L}=-L\dfrac{\partial\ell}{\partial L},\qquad C_{d/D}=D\dfrac{\partial\ell}{\partial d},\qquad C_{M}=M\dfrac{\partial\ell}{\partial M}$", EQ),
 (r"$\bf{Minima}$", HD),
 ("eq", r"C_{D/L}=", [[r"\mp" + deg + r"\dfrac{\mathrm{Amp}_A D^{2}L\left(4D\sin\phi_A\pm3d\cos\phi_A\right)\left(D^{3}H^{2}-\alpha GM\right)}{24\,\alpha GM}", None], [r"\times\dfrac{" + BSL + "}{" + BT + r"^{2}}", None]]),
 ("eq", r"C_{d/D}=", DMIN("1", r"\pm\dfrac{135^\circ}{\pi}")),
 ("eq", r"C_{M}=", [[r"\pm" + deg + r"\dfrac{\mathrm{Amp}_A\left(4D\sin\phi_A\pm3d\cos\phi_A\right)D^{2}H^{2}}{24\,\alpha GM}", None], [r"\times\dfrac{" + BS + "}{" + BT + "}", None]]),
]
P[7] += [(r"$\bf{Maxima}$", HD),
 ("eq", r"C_{M}=", [[r"\pm" + deg + r"\,2\mathrm{Amp}_A\cos\phi_A\,LD^{4}H^{2}", None], [BS, None]] + div_Q(r"3\alpha GM", "")),
]
P[8] = [
 (r"$\bf{Page\ 8:\ compact\ form,\ maxima,\ }C_{D/L}$", HD), ("gap", 0.05),
 ("eq", r"C_{D/L}=", [[r"\pm" + deg + r"\dfrac{2}{3}\mathrm{Amp}_A\cos\phi_A\,D^{4}L", None]] + PLch + div_Q(r"\alpha GM", "^{2}")),
 *note(r"Equal to $-L$ times $\partial\ell_{\max}/\partial L$ (page 5)."),
]
P[9] = [
 (r"$\bf{Page\ 9:\ compact\ form,\ maxima,\ }C_{d/D}$", HD), ("gap", 0.05),
 ("eq", r"C_{d/D}=", [[PdCh[0][0].replace(r"\,DL", r"\,D^{2}L"), None]] + [list(c) for c in PdCh[1:]] + div_Q(r"\alpha GM", "^{2}")),
 *note(r"Equal to $D$ times $\partial\ell_{\max}/\partial d$ (page 3). Baseline values of all $C$: page 11."),
]
# tables
LAB = {"lmax1": r"$\ell_{\max,1}$", "lmax2": r"$\ell_{\max,2}$", "lmin1": r"$\ell_{\min,1}$", "lmin2": r"$\ell_{\min,2}$"}
ORD = ["lmin1", "lmin2", "lmax1", "lmax2"]
rows = [("", "", "analytic", "exact", "ratio")]
for nm in ORD:
    for i, (k, kl) in enumerate([("D", r"$\partial/\partial D$"), ("L", r"$\partial/\partial L$"), ("d", r"$\partial/\partial d$"), ("M", r"$\partial/\partial M$"), ("alpha", r"$\partial/\partial\alpha$")]):
        an, ex = J["ana"][nm][k], J["efd"][k][nm]
        rows.append((LAB[nm] if i == 0 else "", kl, f4(an), f4(ex), "%.2f" % (an/ex)))
P[10] = [(r"$\bf{Page\ 10:\ baseline\ partials,\ analytic\ vs\ exact}$", HD),
        *note(r"Degrees per unit input. exact = central difference (step $10^{-3}$) of the true extrema of the no-boxcar curve with tilt.", bullet=False),
        ("gap", 0.05), ("table", rows, [0, 1.4, 3.4, 5.2, 7.0], 15)]
rows2 = [("", r"$C_{D/L}$  an. / ex.", r"$C_{d/D}$  an. / ex.", r"$C_{M}$  an. / ex.")]
b = J["base"]
for nm in ORD:
    rows2.append((LAB[nm], f4(-b["L"]*J["ana"][nm]["L"]) + " / " + f4(-b["L"]*J["efd"]["L"][nm]),
                  f4(b["D"]*J["ana"][nm]["d"]) + " / " + f4(b["D"]*J["efd"]["d"][nm]),
                  f4(b["M"]*J["ana"][nm]["M"]) + " / " + f4(b["M"]*J["efd"]["M"][nm])))
rows3 = [(r"$d/D$", r"exact $\ell_{\min,1}$", "lin. err", "arccos err", "ex. slope", "lin. slope", "arccos sl.")]
for e in (0.1, 0.2, 0.3, 0.5, 0.7, 0.9):
    i = int(np.argmin(abs(ES - e))); lin, arc, dlin, darc = formulas(ES[i])
    rows3.append(("%.1f" % e, "%.2f" % EX[i, 2], ("%+.3f" % (lin[2] - EX[i, 2])).replace("-", "\u2212"), ("%+.3f" % (arc[0] - EX[i, 2])).replace("-", "\u2212"), "%.1f" % DEX[i, 2], "%.1f" % dlin[2], "%.1f" % darc[0]))
rows4 = [(r"$d/D$", r"exact $\ell_{\min,2}$", "lin. err", "arccos err", "ex. slope", "lin. slope", "arccos sl.")]
for e in (0.1, 0.3, 0.5, 0.7, 0.9):
    i = int(np.argmin(abs(ES - e))); lin, arc, dlin, darc = formulas(ES[i])
    rows4.append(("%.1f" % e, "%.2f" % EX[i, 3], ("%+.3f" % (lin[3] - EX[i, 3])).replace("-", "\u2212"), ("%+.3f" % (arc[1] - EX[i, 3])).replace("-", "\u2212"), ("%.1f" % DEX[i, 3]).replace("-", "\u2212"), ("%.1f" % dlin[3]).replace("-", "\u2212"), ("%.1f" % darc[1]).replace("-", "\u2212")))
C7 = [0, 0.75, 2.45, 3.85, 5.3, 6.7, 8.05]
P[11] = [(r"$\bf{Page\ 11:\ compact\ coefficients;\ arccos\ vs\ linear}$", HD),
        *note(r"Baseline $C$ coefficients (degrees per unit of $\ln(D/L)$, of $d/D$, of $\ln(\alpha GM/D^{3}H^{2})$):", bullet=False),
        ("table", rows2, [0, 1.3, 4.1, 6.7], 15), ("gap", 0.15),
        *note(r"Minima vs $d/D$ (other inputs at baseline). err = formula $-$ exact (deg); slope = $d\ell/d(d/D)$ (deg).", bullet=False),
        ("table", rows3, C7, 14.5), ("gap", 0.1), ("table", rows4, C7, 14.5)]
P[12] = [
 (r"$\bf{Page\ 12:\ notes}$", HD), ("gap", 0.05),
 *note(r"The figure extrema-dl-dd-plot.png shows all four angles and $d\ell/dd$ against $d/D$ from 0 to 0.9: exact curve, first-order formula, and the arccos form for the minima. Nothing is fitted."),
 *note(r"Minima: $d$ is the main driver. $\pm135^\circ/(\pi D)$ is 92% of $\partial\ell_{\min,1}/\partial d$ and 110% of $\partial\ell_{\min,2}/\partial d$ at baseline; the tilt factor makes up the rest, the same at every $d$."),
 *note(r"The exact minima accelerate as $d/D$ grows (slope 47 at $d/D=0.1$, 53 at 0.5, 104 at 0.9 deg per unit $d/D$ for $\ell_{\min,1}$). The arccos form captures part of that; the rest is beyond first order in $d/D$."),
 *note(r"Without the tilt the exact minimum obeys $\cos\ell_{\min}=-3d/4D$ to $0.03^\circ$ at $d/D=0.2$, $0.6^\circ$ at 0.5, but $6.9^\circ$ at 0.9: the geometric relation is a small-$d/D$ result."),
 *note(r"Maxima: $d$ acts only through Amp$_A\times d/D$. The first-order $\ell_{\max,1}$ is poor (its $d$-slope is 1.8$\times$ the exact one at baseline, with a pole at $d/D=0.399$); $\ell_{\max,2}$ is fair up to $d/D\approx0.2$."),
 *note(r"Every displayed partial matches sympy's derivative of the card formula to $3\times10^{-14}$ at random inputs; the analytic baseline values match a finite difference of the card formulas to $3\times10^{-9}$."),
]
for n, items in P.items():
    render(n, items)
