"""Pages: where does each tilt-card angle stop changing with an input?  d(ell)/dX = 0 solved
with sympy (nobox_tilt_stat_sym / _forms, checked to 1e-14); roots at baseline from
nobox_tilt_stat_roots.json; exact checks from nobox_tilt_stat_exact*.json."""
import json, numpy as np
import pagekit as pk
from pagekit import render, note, EQ, NT, HD, TT, B
pk.PREFIX = "extrema-nobox-tilt-stationary"
R = json.load(open("nobox_tilt_stat_roots.json")); E = json.load(open("nobox_tilt_stat_exact.json"))
E2 = json.load(open("nobox_tilt_stat_exact2.json"))
Sx = r"(1-f_u)\left(6L^{3}e^{D/L}-6L^{3}-6DL^{2}-3D^{2}L\right)+f_uD^{3}"
SDx = r"(1-f_u)\left(6L^{2}e^{D/L}-6L^{2}-6DL\right)+3f_uD^{2}"
SLx = r"(1-f_u)\left(18L^{2}e^{D/L}-6DLe^{D/L}-18L^{2}-12DL-3D^{2}\right)"
BS, BSD, BSL = r"\left[" + Sx + r"\right]", r"\left[" + SDx + r"\right]", r"\left[" + SLx + r"\right]"
BT = r"\left[(1-f_u)\left(6L^{3}e^{D/L}-6L^{3}-6DL^{2}-3D^{2}L-D^{3}\right)\right]"
deg = r"\dfrac{180^\circ}{\pi}"
Qch = [[r"\div@L\{\dfrac{6}{D^{3}}", None], [r"-\dfrac{6}{" + Sx + "}", None],
       [r"\mp\dfrac{9d}{2D}@L(\dfrac{4}{D^{3}}", None], [r"-\dfrac{4L+(1-f_u)D}{L\left[" + Sx + r"\right]}@R)@R\}", None]]
def f3(v): return ("%#.4g" % v).rstrip(".").replace("-", "\u2212")
title = r"$\bf{Where\ each\ tilt{-}card\ angle\ is\ stationary\ (no\ boxcar)}$"
P = {}
P[1] = [
 (title, TT), (r"$\bf{Page\ 1:\ set\ }\partial\ell/\partial X=0\bf{\ for\ }X=D,\,L,\,d,\,M\bf{\ and\ solve\ (sympy,\ exact)}$", HD), ("gap", 0.1),
 *note(" ".join([(r"The four angles on the card (already solved from $dA/d\ell=0$):")])),
 (r"$\ell_{\max,1}=" + deg + r"\delta_1,\quad\ell_{\max,2}=180^\circ-" + deg + r"\delta_2$", EQ),
 (r"$\ell_{\min,1}=90^\circ+" + deg + r"\left(\dfrac{3d}{4D}+\delta'_1\right),\quad\ell_{\min,2}=270^\circ-" + deg + r"\left(\dfrac{3d}{4D}+\delta'_2\right)$", EQ),
 ("eq", r"\delta'_{1,2}=", [[r"-\mathrm{Amp}_A@L(\sin\phi_A\pm\dfrac{3d}{4D}\cos\phi_A@R)@L(\dfrac{H^{2}}{\alpha GM}-\dfrac{1}{D^{3}}@R)", None],
                           [r"\div@L(\dfrac{6}{D^{3}}-\dfrac{6}{" + Sx + r"}@R)", None]]),
 *note(r"$\delta_{1,2}$ is the card's $-\mathrm{Amp}_A\cos\phi_A\{\ldots\}/\{\ldots\}$. Upper sign: index 1, lower sign: index 2.", bullet=False),
 *note(r"Physical range: $0<d<D$, $D/L>1$, $0<\alpha GM<D^{3}H^{2}$ (else $H_{\rm obs}/H\leq0$ at the minima)."),
 *note(r"Baseline (code units): $D=1$, $L=D/5.9$, $d=0.1$, $M=0.2306$, $\alpha=G=H=1$, $f_u=0.94$, Amp$_A=-0.06$, $\phi_A=0.35$ rad."),
 ("gap", 0.1), (r"$\bf{Three\ facts\ used\ throughout}$", HD),
 *note(r"(1) $\alpha$ follows $M$: $\partial\ell/\partial\alpha=(M/\alpha)\,\partial\ell/\partial M$, so $\partial\ell/\partial\alpha=0$ exactly where $\partial\ell/\partial M=0$; every $M$ result holds for $\alpha GM$."),
 *note(r"(2) Two signs that never change (series of $e^{D/L}$; the first is the bracket minus $D^{3}$):"),
 ("eq", BT + "=", [[r"6(1-f_u)L^{3}\left(e^{D/L}-1-\dfrac{D}{L}-\dfrac{D^{2}}{2L^{2}}-\dfrac{D^{3}}{6L^{3}}\right)>0", None]]),
 ("eq", r"\dfrac{\partial}{\partial L}" + BS + "=", [[r"-6(1-f_u)L^{2}\sum_{n\geq4}\dfrac{(n-3)}{n!}\left(\dfrac{D}{L}\right)^{n}<0", None]]),
 *note(r"(3) Scaling: $D\,\partial_D\ell+L\,\partial_L\ell+d\,\partial_d\ell+3M\,\partial_M\ell=0$ for all four angles (page 7)."),
]
P[2] = [
 (r"$\bf{Page\ 2:\ }\partial\ell/\partial M=0\ \ (\bf{and\ }\partial\ell/\partial\alpha=0)$", HD), ("gap", 0.1),
 (r"$\bf{Maxima}$", HD),
 ("eq", r"\dfrac{\partial\ell_{\max,1,2}}{\partial M}=", [[r"\pm" + deg + r"\dfrac{\mathrm{Amp}_A\cos\phi_A\,H^{2}}{\alpha GM^{2}}", None]] + Qch),
 *note(" ".join([(r"The numerator is never zero, so there is $\bf{no\ solution}$ for either maximum."), (r"($M$ sits only in the numerator brace of $\delta$, as $H^{2}/\alpha GM$.)")])),
 ("gap", 0.15), (r"$\bf{Minima}$", HD),
 ("eq", r"\dfrac{\partial\ell_{\min,1,2}}{\partial M}=", [[r"\pm" + deg + r"\dfrac{\mathrm{Amp}_A\left(4D\sin\phi_A\pm3d\cos\phi_A\right)D^{2}H^{2}}{24\,\alpha GM^{2}}", None],
                                                    [r"\times\dfrac{" + BS + "}{" + BT + "}", None]]),
 (r"Both brackets are positive (page 1), so this vanishes only where", NT),
 (r"$d=\mp\dfrac{4D}{3}\tan\phi_A$", EQ),
 *note(" ".join([(r"$\ell_{\min,1}$: $d=-\frac{4}{3}D\tan\phi_A<0$, so $\bf{no\ solution}$.")])),
 *note(" ".join([(r"$\ell_{\min,2}$: $d=\frac{4}{3}D\tan\phi_A$ ($d/D=0.487$), for every $D$, $L$, $M$. There $\delta'_2\equiv0$, so")])),
 (r"$\ell_{\min,2}=270^\circ-" + deg + r"\tan\phi_A=249.09^\circ$", EQ),
 *note(r"does not depend on $D$, $L$ or $M$ on that line ($\partial\ell_{\min,2}/\partial L=0$ there too, page 4).", bullet=False),
 *note(" ".join([(r"Exact check at $D$, $L$, $M$ baseline: the exact $\ell_{\min,2}$ turns in $M$ at $d=%s$" % f3(E["lmin2_M__d"]["exact_root"])), (r"(formula $0.4867$), and the exact $\ell_{\min,2}$ there is $248.2^\circ$.")])),
]
P[3] = [
 (r"$\bf{Page\ 3:\ }\partial\ell/\partial d=0$", HD), ("gap", 0.1),
 *note(" ".join([(r"Neither condition contains $d$, so it never fixes $d$; it fixes $M$ instead:"), (r"the numerator and denominator braces of $\delta$ are both linear in $d$, so $\partial\delta/\partial d$ has a"), (r"$d$-free numerator; and $\frac{3d}{4D}+\delta'$ is linear in $d$, so its $d$-derivative is $d$-free.")])),
 ("gap", 0.1), (r"$\bf{Maxima}$ (one relation for both)", HD),
 ("eq", r"M=", [[r"\dfrac{H^{2}}{\alpha G}" + BS, None], [r"\times@L\{(1-f_u)D^{7}+4D^{6}L", None], [r"-4D^{3}L" + BS + r"@R\}", None],
                [r"\div@L\{(1-f_u)D^{7}+8D^{6}L", "!"], [r"-12D^{3}L" + BS, None], [r"+4L" + BS + r"^{2}@R\}", None]]),
 *note(" ".join([(r"$\alpha GM/D^{3}H^{2}$ from this depends on $D/L$ (and $f_u$) only: $\geq1.10$ for $D/L<6.71$, infinite at"), (r"6.71, negative (towards $-1$) above; it is $%.2f$ at $D/L=5.9$. Never in $(0,1)$: $\bf{no\ physical\ solution}$." % R["info"]["kap_maxd_base"])])),
 ("gap", 0.1), (r"$\bf{Minima}$", HD),
 ("eq", r"M=", [[r"\pm\dfrac{D^{3}H^{2}}{\alpha G}\mathrm{Amp}_A\cos\phi_A" + BS, None], [r"\div@L\{\pm\mathrm{Amp}_A\cos\phi_A" + BS, "!"], [r"+6" + BT + r"@R\}", None]]),
 *note(" ".join([(r"That is $\alpha GM/D^{3}H^{2}=1/\left(1\pm6[\ldots-D^{3}]/(\mathrm{Amp}_A\cos\phi_A[\ldots])\right)$ with $[\ldots-D^{3}]>0$ and"), (r"Amp$_A\cos\phi_A<0$: for $\ell_{\min,1}$ it is never in $(0,1)$ ($\bf{none}$); for $\ell_{\min,2}$ it always is.")])),
 *note(r"$\ell_{\min,2}$ at baseline $D$, $L$: $M=%s$. (Read as a root in $D$ or $L$: $D=%s$, $L=%s$.)" % (
     f3(R["info"]["M_min2_d"]), f3(R["roots"]["lmin2_d"]["D"][0][0]), f3(R["roots"]["lmin2_d"]["L"][0][0]))),
 *note(r"Exact: the second minimum does not exist there. It merges with the second maximum below $M=%.4f$ (or above $D=%.3f$, $L=%.3f$), and just before that $\partial\ell_{\min,2}/\partial d\approx{-}%.1f$. So the exact $\ell_{\min,2}$ never turns in $d$; this formula root marks where the minimum disappears." % (
     E2["Mc_d0.1"], E2["Dc"], E2["Lc"], -E2["Mc_d0.1_dldd_near_0.01"])),
]
P[4] = [
 (r"$\bf{Page\ 4:\ }\partial\ell/\partial L=0$", HD), ("gap", 0.1),
 (r"$\bf{Minima}$", HD),
 ("eq", r"\dfrac{\partial\ell_{\min,1,2}}{\partial L}=", [[r"\pm" + deg + r"\dfrac{\mathrm{Amp}_A D^{2}\left(4D\sin\phi_A\pm3d\cos\phi_A\right)\left(D^{3}H^{2}-\alpha GM\right)}{24\,\alpha GM}", None],
                                                    [r"\times\dfrac{" + BSL + "}{" + BT + r"^{2}}", None]]),
 *note(" ".join([(r"The $L$-derivative of the bracket is never zero (page 1), so this vanishes only on $d=\mp\frac{4}{3}D\tan\phi_A$"), (r"(page 2: $\ell_{\min,2}$ only) or at $\alpha GM=D^{3}H^{2}$, the edge of the physical range. No root in $L$.")])),
 *note(" ".join([(r"Exact check ($D$, $L$, $M$ baseline): $\ell_{\min,2}$ turns in $L$ at $d=%s$ (formula 0.4867)." % f3(E["lmin2_L__d"]["exact_root"]))])),
 ("gap", 0.1), (r"$\bf{Maxima}$", HD),
 *note(r"Linear in $H^{2}/\alpha GM$, so explicit in $M$; transcendental in $D/L$ through $e^{D/L}$:", bullet=False),
 ("eq", r"M=", [[r"\dfrac{D^{3}H^{2}}{\alpha G}@L\{4DL^{2}" + BSL, None],
                [r"\mp3d@L[@L((1-f_u)D+4L@R)L" + BSL, "!"], [r"+(1-f_u)D" + BS + r"@R]@R\}", None],
                [r"\div@L\{4DL^{2}" + BSL, "!"], [r"\mp3d@L((1-f_u)D^{4}", None], [r"+4L^{2}" + BSL + r"@R)", None],
                [r"-9(1-f_u)d^{2}@L(D^{3}-L" + BSL, "!"], [r"-" + BS + r"@R)@R\}", None]]),
 *note(" ".join([(r"Root in $L$ (others at baseline): $\bf{none}$ for $0.02<L<1.5$, either maximum.")])),
 *note(r"$\ell_{\max,1}$, as a root in another input: $D=%s$, $d=%s$, $M=%s$. $\ell_{\max,2}$: $M=%s$ (just past the edge $\alpha GM=D^{3}H^{2}$), no physical $D$ or $d$." % (
     f3(R["roots"]["lmax1_L"]["D"][0][0]), f3(R["roots"]["lmax1_L"]["d"][0][0]), f3(R["roots"]["lmax1_L"]["M"][0][0]), "1.010")),
]
P[5] = [
 (r"$\bf{Page\ 5:\ }\partial\ell_{\max}/\partial D=0$", HD), ("gap", 0.1),
 *note(r"Linear in $H^{2}/\alpha GM$, so explicit in $M$; transcendental in $D/L$:", bullet=False),
 ("eq", r"M=", [[r"\dfrac{D^{3}H^{2}}{\alpha G}@L\{4DL@L(D^{4}" + BSD, None], [r"-3" + BS + r"^{2}@R)", "!"],
                [r"\mp3d@L((1-f_u)D^{5}" + BSD, "!"], [r"+4D^{4}L" + BSD, "!"], [r"+4D^{3}L" + BS, None], [r"-16L" + BS + r"^{2}@R)@R\}", None],
                [r"\div@L\{4D^{4}L@L(D" + BSD, "!"], [r"-3" + BS + r"@R)", None],
                [r"\mp12dL@L(2D^{6}+D^{4}" + BSD, "!"], [r"-6D^{3}" + BS, None], [r"+" + BS + r"^{2}@R)", None],
                [r"+9(1-f_u)d^{2}@L(D^{6}+D^{4}" + BSD, "!"], [r"-4D^{3}" + BS + r"@R)@R\}", None]]),
 *note(" ".join([(r"Root in $D$ (others at baseline): $\bf{none}$ for $0.2<D<3$, either maximum. $M$ from this is"), (r"negative at baseline ($%s$, $%s$), so no physical $M$ either." % (f3(-0.2257), f3(-0.3344)))])),
 *note(" ".join([(r"Root in $d$: $\ell_{\max,1}$ $d=%s$; $\ell_{\max,2}$ none ($d<0$)." % f3(R["roots"]["lmax1_D"]["d"][0][0]))])),
]
P[6] = [
 (r"$\bf{Page\ 6:\ }\partial\ell_{\min}/\partial D=0$", HD), ("gap", 0.1),
 (r"Linear in $d$ (and in $H^{2}/\alpha GM$), so explicit in $d$:", NT),
 ("eq", r"d=", [[r"-\dfrac{4}{3}\mathrm{Amp}_A D^{4}\sin\phi_A", None], [r"@L\{D@L(D^{3}H^{2}-\alpha GM@R)" + BSD, "!"],
                [r"-3" + BS, "!"], [r"\times@L(H^{2}" + BS + r"-\alpha GM@R)@R\}", None],
                [r"\div@L\{\pm\mathrm{Amp}_A\cos\phi_A", "!"], [r"\times@L[D^{4}@L(D^{3}H^{2}-\alpha GM@R)" + BSD, None],
                [r"-D^{3}@L(D^{3}H^{2}-4\alpha GM@R)" + BS, "!"], [r"-@L(2D^{3}H^{2}+\alpha GM@R)" + BS + r"^{2}@R]", "!"],
                [r"-6\alpha GM" + BT + r"^{2}@R\}", None]]),
 (r"At baseline, as a root in each input (the others at baseline):", NT),
]
def rr(k, Y, i=-1): return f3(R["roots"][k][Y][i][0])
P[6] += [(B + r"$\ell_{\min,1}$: $d=%s$, $D=%s$, $L=%s$, $M=%s$." % (rr("lmin1_D", "d"), rr("lmin1_D", "D"), rr("lmin1_D", "L"), rr("lmin1_D", "M")), NT),
         *note(r"(A second root in $D$, $D=%s$, has $\alpha GM>D^{3}H^{2}$: unphysical.)" % rr("lmin1_D", "D", 0), bullet=False),
         *note(" ".join([(r"$\ell_{\min,2}$: $d=%s$, $D=%s$, $L=%s$, $M=%s$." % (rr("lmin2_D", "d"), rr("lmin2_D", "D"), rr("lmin2_D", "L"), rr("lmin2_D", "M")))])),
         (B + r"Exact turning points: page 9.", NT)]
P[7] = [
 (r"$\bf{Page\ 7:\ which\ conditions\ are\ independent;\ common\ points}$", HD), ("gap", 0.1),
 (r"$D\,\dfrac{\partial\ell}{\partial D}+L\,\dfrac{\partial\ell}{\partial L}+d\,\dfrac{\partial\ell}{\partial d}+3M\,\dfrac{\partial\ell}{\partial M}=0$", EQ),
 *note(" ".join([(r"Each angle depends only on $D/L$, $d/D$ and $\alpha GM/D^{3}H^{2}$, so only three of the four"), (r"conditions are independent: any three imply the fourth.")])),
 *note(" ".join([(r"All four at once means $\ell$ is stationary in all three ratios.")])),
 ("gap", 0.1),
 *note(" ".join([(r"$\ell_{\max,1}$, $\ell_{\max,2}$: $\partial\ell/\partial M\neq0$ everywhere (page 2), so the $D$, $L$, $d$ conditions"), (r"can never hold together. $\bf{No\ common\ point.}$")])),
 *note(" ".join([(r"$\ell_{\min,1}$: $\partial\ell/\partial M\neq0$ for $d>0$. $\bf{No\ common\ point.}$")])),
 *note(" ".join([(r"$\ell_{\min,2}$: $d=\frac{4}{3}D\tan\phi_A$ gives $\partial_M\ell=\partial_L\ell=0$; the $M$ relation of page 3 gives $\partial_d\ell=0$;"), (r"then $\partial_D\ell=0$ follows from the identity. Common set (a one-parameter family in $D/L$):")])),
 (r"$d=\dfrac{4D}{3}\tan\phi_A$", EQ),
 ("eq", r"M=", [[r"-\dfrac{D^{3}H^{2}}{\alpha G}\mathrm{Amp}_A\cos\phi_A" + BS, None], [r"\div@L\{-\mathrm{Amp}_A\cos\phi_A" + BS, "!"], [r"+6" + BT + r"@R\}", None]]),
 *note(" ".join([(r"At $D/L=5.9$: $d/D=0.487$, $\alpha GM/D^{3}H^{2}=%.4f$, $\ell_{\min,2}=249.09^\circ$." % R["info"]["M_min2_d"])])),
 *note(" ".join([(r"Exact: on that line the second minimum vanishes below $\alpha GM/D^{3}H^{2}=%.4f$, 2%% above the" % E2["Mc_dline"]), (r"formula's point, so there is no exact common point; the formula's point lies just past the merger."), (r"Just above the merger the exact $\ell_{\min,2}=%.2f^\circ$ (formula $249.09^\circ$)." % E2["Mc_dline_lmin2_near_0.01"])])),
]
# ---------- summary tables ----------
def ex(key):
    e = E.get(key); return e
def row(lab, X, sol, Y, key=None, form=None, exact=None, slope=None):
    if key:
        e = E[key]; form = f3(R["checks"][key]["y0"])
        exact = f3(e["exact_root"]) if e["exact_root"] else ("none" if np.isfinite(e["at_root"]) else "no $\\ell$ there")
        slope = f3(e["at_root"]) if np.isfinite(e["at_root"]) else "\u2014"
    return (lab, X, sol, Y, form or "", exact or "", slope or "")
lmin2DM = E2["lmin2_D_M_exact_roots"][0]
T1 = [("angle", "set 0", "solution", "in", "formula", "exact", "exact slope")]
T1 += [row(r"$\ell_{\max,1}$", r"$\partial_D$", r"$M$ (p5)", r"$D$", form="none", exact="\u2014", slope=""),
       row("", "", "", r"$M$", form="<0", exact="\u2014"),
       row("", "", "", r"$d$", "lmax1_D__d"),
       row("", r"$\partial_L$", r"$M$ (p4)", r"$L$", form="none", exact="\u2014"),
       row("", "", "", r"$D$", "lmax1_L__D"), row("", "", "", r"$d$", "lmax1_L__d"), row("", "", "", r"$M$", "lmax1_L__M"),
       row("", r"$\partial_d$", r"$M$ (p3)", r"$M$", form="unphys.", exact="\u2014"),
       row("", r"$\partial_M$", "none", "", form="", exact=""),
       row(r"$\ell_{\max,2}$", r"$\partial_D$", r"$M$ (p5)", "all", form="none", exact="\u2014"),
       row("", r"$\partial_L$", r"$M$ (p4)", r"$M$", form="1.010 (unphys.)", exact="\u2014"),
       row("", r"$\partial_d$", r"$M$ (p3)", r"$M$", form="unphys.", exact="\u2014"),
       row("", r"$\partial_M$", "none", "", form="", exact="")]
T2 = [("angle", "set 0", "solution", "in", "formula", "exact", "exact slope")]
T2 += [row(r"$\ell_{\min,1}$", r"$\partial_D$", r"$d$ (p6)", r"$D$", "lmin1_D__D"), row("", "", "", r"$L$", "lmin1_D__L"),
       row("", "", "", r"$d$", "lmin1_D__d"), row("", "", "", r"$M$", "lmin1_D__M"),
       row("", r"$\partial_L$", "edge only", "", form="none", exact=""),
       row("", r"$\partial_d$", r"$M$ (p3)", r"$M$", form="unphys.", exact="\u2014"),
       row("", r"$\partial_M$", r"$d<0$", "", form="none", exact=""),
       row(r"$\ell_{\min,2}$", r"$\partial_D$", r"$d$ (p6)", r"$D$", "lmin2_D__D"), row("", "", "", r"$L$", "lmin2_D__L"),
       row("", "", "", r"$d$", "lmin2_D__d"),
       row("", "", "", r"$M$", form=rr("lmin2_D", "M"), exact=f3(lmin2DM), slope=f3(E["lmin2_D__M"]["at_root"])),
       row("", r"$\partial_L$", r"$d=\frac{4}{3}D\tan\phi_A$", r"$d$", "lmin2_L__d"),
       row("", r"$\partial_d$", r"$M$ (p3)", r"$M$", "lmin2_d__M"), row("", "", "", r"$D$", "lmin2_d__D"), row("", "", "", r"$L$", "lmin2_d__L"),
       row("", r"$\partial_M$", r"$d=\frac{4}{3}D\tan\phi_A$", r"$d$", "lmin2_M__d")]
COLS = [0, 0.95, 1.8, 4.0, 4.65, 6.55, 8.0]
P[8] = [(r"$\bf{Page\ 8:\ summary,\ maxima}$", HD),
        *note(r"Root of $\partial\ell/\partial X=0$ in the input named under 'in', the others at baseline. 'exact' = where the exact $\ell$ turns; 'exact slope' = exact $\partial\ell/\partial X$ (deg per unit) at the formula root.", bullet=False),
        ("gap", 0.05), ("table", T1, COLS, 15), ("gap", 0.2),
        *note(" ".join([(r"none = no root in the range searched; unphys. = root outside $0<\alpha GM<D^{3}H^{2}$ (pages 3, 4).")])),
        *note(" ".join([(r"$\ell_{\max,1}$ $d$-roots (0.286, 0.316) sit next to $d=0.399$, where the denominator brace"), (r"of $\delta_1$ is zero; the exact $\ell_{\max,1}$ does not turn within $\pm50\%$ of them.")])),
        (B + r"$\partial\ell/\partial\alpha=0$ wherever $\partial\ell/\partial M=0$: no separate rows.", NT)]
P[9] = [(r"$\bf{Page\ 9:\ summary,\ minima}$", HD), (r"Same columns as page 8.", NT), ("gap", 0.05), ("table", T2, COLS, 15), ("gap", 0.2),
        *note(" ".join([(r"no $\ell$ there = the exact second minimum does not exist at that input value (page 3).")])),
        *note(r"Exact zeros: Brent root of the exact $\partial\ell/\partial X$ (central difference, step $10^{-3}$), bracketed within $\pm50\%$ of the formula root ($\ell_{\min,2}$, $\partial_D$ in $M$: $0.04<M<0.95$).")]
P[10] = [
 (r"$\bf{Page\ 10:\ what\ stands\ out}$", HD), ("gap", 0.1),
 *note(" ".join([(r"$\partial\ell/\partial d=0$ never gives $d$, for any angle: the condition contains no $d$. It gives $M$, and"), (r"only $\ell_{\min,2}$ has a physical solution ($M=0.0262$ at baseline $D$, $L$).")])),
 *note(r"That $\ell_{\min,2}$ root sits where $\mathrm{Amp}_A\cos\phi_A\,(H^{2}/\alpha GM-1/D^{3})=-(6/D^{3}-6/[\ldots])$ (the tilt is no longer small): the exact second minimum merges with the second maximum before that (at $M=0.0439$), so the exact $\ell_{\min,2}$ never turns in $d$."),
 *note(" ".join([(r"No maximum is ever stationary in $M$ (or $\alpha$), and $\ell_{\max,2}$ has no physical stationary point"), (r"in any input at baseline.")])),
 *note(" ".join([(r"$\ell_{\min,2}$ is frozen at $249.09^\circ$ on $d=\frac{4}{3}D\tan\phi_A$, whatever $D$, $L$, $M$; the exact $\ell_{\min,2}$"), (r"is $248.2$–$249.1^\circ$ there, and turns in $M$ at $d=0.448$ and in $L$ at $d=0.434$.")])),
 *note(" ".join([(r"$\partial\ell_{\min}/\partial D=0$ is the best-tracked: the exact roots are within 5–6% ($D$), 10% ($L$), 11–13% ($d$);"), (r"the $M$ roots are poorer (0.133 vs 0.182 for $\ell_{\min,1}$; 0.125 vs 0.062 for $\ell_{\min,2}$).")])),
 *note(" ".join([(r"$\ell_{\max,1}$ is stationary in $L$ only near the edge $\alpha GM\approx D^{3}H^{2}$ ($M=0.966$ vs exact 0.984).")])),
 *note(" ".join([(r"Displayed solutions agree with sympy's own solutions to $10^{-13}$ at random inputs.")])),
]
for n, items in P.items():
    render(n, items)
