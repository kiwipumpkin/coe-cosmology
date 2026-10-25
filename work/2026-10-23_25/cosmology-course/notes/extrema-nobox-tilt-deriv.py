"""Pages: partial derivatives of the 8 no-boxcar tilt-card quantities wrt D, L, d, M (and alpha).
Derivatives come from nobox_tilt_deriv_sym (sympy); every displayed structure is checked
numerically in nobox_tilt_deriv_check. 10 in wide, dpi 250, equations 18 pt, notes 15 pt."""
import json, re, numpy as np, sympy as sp, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from nobox_tilt_deriv_sym import *
plt.rcParams.update({"mathtext.fontset": "cm", "font.family": "serif"})
J = json.load(open("nobox_tilt_deriv.json"))
EQ, NT, HD, TT = 18, 15, 17, 20
B = r"$\bullet$  "
Sx = r"(1-f_u)\left(6L^{3}e^{D/L}-6L^{3}-6DL^{2}-3D^{2}L\right)+f_uD^{3}"
SDx = r"(1-f_u)\left(6L^{2}e^{D/L}-6L^{2}-6DL\right)+3f_uD^{2}"
SLx = r"(1-f_u)\left(18L^{2}e^{D/L}-6DLe^{D/L}-18L^{2}-12DL-3D^{2}\right)"
BT = r"\left[(1-f_u)\left(6L^{3}e^{D/L}-6L^{3}-6DL^{2}-3D^{2}L-D^{3}\right)\right]"
BS, BSD, BSL = r"\left[" + Sx + r"\right]", r"\left[" + SDx + r"\right]", r"\left[" + SLx + r"\right]"
w = sp.Symbol("w", positive=True)
NAMES = {w: r"(1-f_u)", al: r"\alpha", a: r"\mathrm{Amp}_A", ph: r"\phi_A"}
def lx(e):
    t = sp.latex(e, symbol_names=NAMES, mul_symbol=None)
    return t.replace(r"\frac", r"\dfrac").replace(r"\left(1-f_u\right)", "(1-f_u)")

# ---------- term printer for collected (coef, monomial) lists ----------
def mono_parts(m):
    num, den = [], []
    for fac in sp.Mul.make_args(m):
        b, e = fac.as_base_exp()
        st = {S: BS, SD: BSD, SL: BSL}.get(b)
        if st is None: continue
        (num if e > 0 else den).append(st + ("^{%d}" % abs(e) if abs(e) != 1 else ""))
    return " ".join(num) or "1", " ".join(den)
def signed(c):
    """(sign_str, magnitude latex) for an s-free coefficient."""
    if c.could_extract_minus_sign(): return "-", -c
    return "+", c
def term_chunks(coef, m):
    coef = sp.expand(coef.subs(f, 1 - w))
    c0 = sp.factor(coef.subs(s, 0)); c1 = sp.factor(coef.coeff(s, 1))
    mn, md = mono_parts(m)
    out = []
    def merged(sg_sym, c):
        n_, d_ = sp.fraction(sp.factor(c))
        nn = lx(n_) if mn == "1" else (mn if n_ == 1 else lx(n_) + " " + mn)
        dd = " ".join(x for x in [lx(d_) if d_ != 1 else "", md] if x)
        return sg_sym + (r"\dfrac{%s}{%s}" % (nn, dd) if dd else nn)
    if m == 1 or c0 == 0 or c1 == 0:
        if c0 != 0:
            sg, mag = signed(c0); out.append([merged(sg, mag), None])
        if c1 != 0:
            sg, mag = signed(c1); out.append([merged(r"\mp" if sg == "-" else r"\pm", mag), None])
        return out
    sg0, m0 = signed(c0); sg1, m1 = signed(c1)
    if sg0 == "-": sg1 = "+" if sg1 == "-" else "-"
    pm = r"\mp" if sg1 == "-" else r"\pm"
    out.append([sg0 + r"\left(" + lx(m0) + pm + lx(m1) + r"\right)", None])
    out.append([r"\dfrac{%s}{%s}" % (mn, md) if md else mn, r"\times"])
    return out
def order_key(item):
    m = item[1]; nS = -sp.degree(m.subs(S, 1/sp.Symbol("z")), sp.Symbol("z")) if m != 1 else 0
    return (m.has(SD) or m.has(SL), m != 1, str(m))
def expr_chunks(expr):
    ch = []
    for c, m in sorted(collect_terms(expr), key=order_key):
        ch += term_chunks(c, m)
    if ch and ch[0][0].startswith("+"): ch[0][0] = ch[0][0][1:]
    return ch

# card's Q (denominator brace of delta), written as on the card
Qch = [[r"\dfrac{6}{D^{3}}", None], [r"-\dfrac{6}{" + Sx + "}", None],
       [r"\mp\dfrac{9d}{2D}@L(\dfrac{4}{D^{3}}", None], [r"-\dfrac{4L+(1-f_u)D}{L\left[" + Sx + r"\right]}@R)", None]]

# ---------- page layout with automatic wrapping ----------
WMAX = 9.25; IND = 0.55
fig0 = plt.figure(figsize=(30, 3)); R0 = fig0.canvas.get_renderer()
def width(t, sz):
    tx = fig0.text(0, .5, t, fontsize=sz); b = tx.get_window_extent(R0); tx.remove(); return b.width/fig0.dpi
def balance(body, carry):
    """Convert @L( / @R) markers; close/open across lines with \\left. / \\right."""
    toks = re.findall(r"@L(\\\{|\(|\[)|@R(\\\}|\)|\])", body)
    opens = carry
    for o, c in toks:
        opens += 1 if o else -1
    pre = r"\left." * carry
    body = re.sub(r"@L(\\\{|\(|\[)", r"\\left\1", body); body = re.sub(r"@R(\\\}|\)|\])", r"\\right\1", body)
    post = r"\right." * opens
    return pre + body + post, opens
def wrap(lhs, chunks, sz=EQ):
    lines, cur, carry, lvl = [], lhs, 0, 0
    def txt(body, c): return "$" + balance(body, c)[0] + "$"
    for ch, cont in chunks:
        trial = cur + " " + ch
        if cont != "!" and width(txt(trial, carry), sz) + lvl*IND <= WMAX:
            cur = trial
        else:
            t, carry_new = balance(cur, carry)
            lines.append(("\x01"*lvl + "$" + t + "$", sz)); carry = carry_new
            lvl = 1
            cur = (cont + " " if cont and cont != "!" else "") + ch
    t, _ = balance(cur, carry); lines.append(("\x01"*lvl + "$" + t + "$", sz))
    return lines
def render(n, items, dpi=250):
    fig = plt.figure(figsize=(10, 40)); r = fig.canvas.get_renderer()
    H0 = 40; y = H0 - 0.35; x0 = 0.35; widest = 0; extra = []
    flat = []
    for it in items:
        if it[0] == "eq": flat += wrap(it[1], it[2], it[3] if len(it) > 3 else EQ)
        else: flat.append(it)
    for it in flat:
        if it[0] == "gap": y -= it[1]; continue
        if it[0] == "table":
            rows, cols, fs = it[1], it[2], it[3]
            for j, row in enumerate(rows):
                for c, t in zip(cols, row):
                    fig.text((x0 + c)/10, y/H0, t, fontsize=fs, va="top", fontweight="bold" if j == 0 else "normal")
                y -= 0.36 if j else 0.42
                if j == 0: extra.append(y + 0.05)
            continue
        t, sz = it
        ind = IND*(len(t) - len(t.lstrip("\x01"))); t = t.lstrip("\x01")
        tx = fig.text((x0 + ind)/10, y/H0, t, fontsize=sz, va="top")
        bb = tx.get_window_extent(r); wd = bb.width/fig.dpi + ind; h = bb.height/fig.dpi
        widest = max(widest, wd)
        if wd > 9.4: print(f"  p{n} TOO WIDE {wd:.2f}: {t[:80]}")
        y -= h + (0.26 if sz >= EQ else 0.13)
    used = H0 - y + 0.25
    for yy in extra:
        fig.add_artist(plt.Line2D([x0/10, 9.65/10], [yy/H0]*2, color="k", lw=0.8))
    fig.set_size_inches(10, used, forward=False)
    for t in fig.texts: t.set_y(1 - (1 - t.get_position()[1])*H0/used)
    for l in fig.artists: l.set_ydata([1 - (1 - v)*H0/used for v in l.get_ydata()])
    fig.savefig(f"extrema-nobox-tilt-deriv-p{n}.png", dpi=dpi); plt.close(fig)
    print(f"p{n}: 10 x {used:.2f} in, widest {widest:.2f} in" + ("   <-- TALL" if used > 12.2 else ""))
    return used

def dpart(expr, X): return expr_chunks(pd(expr, X))
def ddelta_eq(X, lhs):
    ch = [[r"-@L\{\mathrm{Amp}_A\cos\phi_A@L[", None]]
    nch = dpart(N, X); ch[0][0] += " " + nch[0][0]; ch += nch[1:]
    ch[-1][0] += r"@R]"
    qch = dpart(Q, X)
    ch.append([r"+\,\delta_{1,2}@L[" + qch[0][0], "!"]); ch += qch[1:]
    ch[-1][0] += r"@R]@R\}"
    ch.append([r"\div@L\{" + Qch[0][0], "!"]); ch += [list(c) for c in Qch[1:]]
    ch[-1][0] += r"@R\}"
    return ("eq", lhs, ch)

title = r"$\bf{Derivatives\ of\ the\ tilt\ card\ (no\ boxcar)}$"
P = {}
P[1] = [
 (title, TT), (r"$\bf{Page\ 1:\ setup,\ shortcuts,\ building\ blocks}$", HD), ("gap", 0.1),
 (B + r"Eight quantities: $\ell_{\max,1},\ \ell_{\max,2},\ \ell_{\min,1},\ \ell_{\min,2}$ and their heights $A/A_{\rm all}$,", NT),
 (r"$\qquad$differentiated exactly as written on the card (extrema-nobox-tilt-p*.png).", NT),
 (B + r"Inputs $D$ (distance from the centre of expansion (CoE) to the observer),", NT),
 (r"$\qquad L$ (core scale length), $d$ (source-shell radius), $M$ (mass inside $D$), and $\alpha$.", NT),
 (B + r"$\delta_{1,2}$, $\delta'_{1,2}$, $A_{\rm all}$ and $(A/A_{\rm all})$ are the card's own quantities, used by name.", NT),
 (B + r"Units: the model code's units, $D=1$, $G=H=1$, so $L=D/5.9$, $d=0.1$, $M=0.2306\,D^{3}H^{2}/G$.", NT),
 (r"$\qquad$Angles in degrees; derivatives per unit of $D$ (length) and of $D^{3}H^{2}/G$ (mass).", NT),
 ("gap", 0.2),
 (r"$\bf{Shortcuts\ (exact\ for\ all\ eight)}$", HD),
 (r"$\dfrac{\partial}{\partial\alpha}=\dfrac{M}{\alpha}\dfrac{\partial}{\partial M}$   ($\alpha$ only appears as $\alpha GM$; likewise $\partial/\partial G$)", EQ),
 (r"$D\dfrac{\partial}{\partial D}+L\dfrac{\partial}{\partial L}+d\dfrac{\partial}{\partial d}+3M\dfrac{\partial}{\partial M}=0$", EQ),
 (r"(every card quantity depends only on $D/L$, $d/D$ and $\alpha GM/D^{3}H^{2}$)", NT),
 ("gap", 0.2),
 (r"$\bf{Building\ blocks}$", HD),
]
P[1] += [("eq", r"\dfrac{\partial}{\partial D}" + BS + "=", [[SDx, None]]),
         ("eq", r"\qquad=", [[r"3D^{2}\left[f_u+\dfrac{2(1-f_u)L^{2}}{D^{2}}\left(e^{D/L}-1-\dfrac{D}{L}\right)\right]", None]]),
         ("eq", r"\dfrac{\partial}{\partial L}" + BS + "=", [[SLx, None]]),
         (r"Check: $D\,\partial_D+L\,\partial_L$ of the bracket $=3\times$ the bracket.", NT)]
P[2] = [
 (r"$\bf{Page\ 2:\ }A_{\rm all}\bf{\ and\ the\ short\ cases}$", HD), ("gap", 0.1),
 ("eq", r"\dfrac{\partial A_{\rm all}}{\partial M}=", [[r"\dfrac{\alpha G}{2H^{2}}@L(\dfrac{1}{D^{3}}", None], [r"-\dfrac{3}{" + Sx + r"}@R)", None]]),
 (r"$\qquad=\dfrac{A_{\rm all}-1}{M}$", EQ),
 ("eq", r"\dfrac{\partial A_{\rm all}}{\partial D}=", [[r"\dfrac{3\alpha GM}{2H^{2}}@L(-\dfrac{1}{D^{4}}", None], [r"+\dfrac{" + BSD + "}{" + BS + r"^{2}}@R)", None]]),
 ("eq", r"\dfrac{\partial A_{\rm all}}{\partial L}=", [[r"\dfrac{3\alpha GM}{2H^{2}}", None], [r"\dfrac{" + BSL + "}{" + BS + r"^{2}}", r"\times"]]),
 (r"$\dfrac{\partial A_{\rm all}}{\partial d}=0$  exactly.", EQ),
 ("gap", 0.2),
 (r"$\bf{Short\ facts}$", HD),
 (B + r"Without the tilt the maxima sit at exactly $0^\circ$ and $180^\circ$, so $\partial\ell_{\max}/\partial d=0$;", NT),
 (r"$\qquad$here $\partial\ell_{\max}/\partial d$ comes only from the Amp$_A\times d/D$ cross terms in $\delta_{1,2}$.", NT),
 (B + r"$\partial\ell_{\min,1}/\partial d=135^\circ/(\pi D)+(180^\circ/\pi)\,\partial\delta'_1/\partial d$;", NT),
 (r"$\qquad\partial\ell_{\min,2}/\partial d=-135^\circ/(\pi D)-(180^\circ/\pi)\,\partial\delta'_2/\partial d$.", NT),
 (B + r"$M$ enters $\delta_{1,2}$ and $\delta'_{1,2}$ only through $H^{2}/\alpha GM$, and the heights only through $\alpha GM$.", NT),
 (B + r"$\partial\delta_{1,2}/\partial M$ has no denominator-derivative part; $\partial\delta'_{1,2}/\partial L$ comes only", NT),
 (r"$\qquad$from the bracket in its denominator.", NT),
]
P[3] = [
 (r"$\bf{Page\ 3:\ }\delta_{1,2}\bf{\ with\ respect\ to\ }M\bf{\ and\ }d$", HD),
 (r"(radians; upper sign $\delta_1$, lower sign $\delta_2$, as on the card)", NT), ("gap", 0.1),
 ("eq", r"\dfrac{\partial\delta_{1,2}}{\partial M}=", [[r"\dfrac{\mathrm{Amp}_A\cos\phi_A\,H^{2}}{\alpha GM^{2}}", None]] +
   [[r"\div@L\{" + Qch[0][0], None]] + Qch[1:-1] + [[Qch[-1][0] + r"@R\}", None]]),
 ("gap", 0.15),
 ddelta_eq(d, r"\dfrac{\partial\delta_{1,2}}{\partial d}="),
]
P[4] = [(r"$\bf{Page\ 4:\ }\delta_{1,2}\bf{\ with\ respect\ to\ }D$", HD),
        (r"(upper sign $\delta_1$, lower sign $\delta_2$)", NT), ("gap", 0.1),
        ddelta_eq(D, r"\dfrac{\partial\delta_{1,2}}{\partial D}=")]
P[5] = [(r"$\bf{Page\ 5:\ }\delta_{1,2}\bf{\ with\ respect\ to\ }L$", HD),
        (r"(upper sign $\delta_1$, lower sign $\delta_2$)", NT), ("gap", 0.1),
        ddelta_eq(L, r"\dfrac{\partial\delta_{1,2}}{\partial L}=")]
P[6] = [
 (r"$\bf{Page\ 6:\ }\delta'_{1,2}\bf{\ with\ respect\ to\ each\ input}$", HD),
 (r"(radians; upper sign $\delta'_1$, lower sign $\delta'_2$, as on the card)", NT), ("gap", 0.1),
 (r"$\dfrac{\partial\delta'_{1,2}}{\partial M}=-\dfrac{D^{3}H^{2}\,\delta'_{1,2}}{M\left(D^{3}H^{2}-\alpha GM\right)}$", EQ),
 (r"$\dfrac{\partial\delta'_{1,2}}{\partial d}=\pm\dfrac{3\cos\phi_A\,\delta'_{1,2}}{4D\sin\phi_A\pm3d\cos\phi_A}$", EQ),
 ("eq", r"\dfrac{\partial\delta'_{1,2}}{\partial D}=", [[r"\delta'_{1,2}@L[\mp\dfrac{3d\cos\phi_A}{D\left(4D\sin\phi_A\pm3d\cos\phi_A\right)}", None],
   [r"+\dfrac{3\alpha GM}{D\left(D^{3}H^{2}-\alpha GM\right)}", None],
   [r"+\dfrac{L" + BSL + r"}{D" + BT + r"}", None],
   [r"+\dfrac{" + BSD + "}{" + BS + r"}@R]", None]]),
 ("eq", r"\dfrac{\partial\delta'_{1,2}}{\partial L}=", [[r"\delta'_{1,2}@L[-\dfrac{" + BSL + "}{" + BT + r"}", None],
   [r"+\dfrac{" + BSL + "}{" + BS + r"}@R]", None]]),
 ("gap", 0.1),
 (B + r"The second bracket, $(1-f_u)\left(6L^{3}e^{D/L}-6L^{3}-6DL^{2}-3D^{2}L-D^{3}\right)$,", NT),
 (r"$\qquad$is the first bracket minus $D^{3}$.", NT),
]
P[7] = [
 (r"$\bf{Page\ 7:\ the\ four\ angles}$", HD), (r"($X$ is any of $D$, $L$, $d$, $M$; degrees per unit of $X$)", NT), ("gap", 0.1),
 (r"$\dfrac{\partial\ell_{\max,1}}{\partial X}=\dfrac{180^\circ}{\pi}\,\dfrac{\partial\delta_1}{\partial X},\qquad\dfrac{\partial\ell_{\max,2}}{\partial X}=-\dfrac{180^\circ}{\pi}\,\dfrac{\partial\delta_2}{\partial X}$", EQ),
 (r"$\dfrac{\partial\ell_{\min,1}}{\partial X}=\dfrac{180^\circ}{\pi}\left(\dfrac{\partial}{\partial X}\dfrac{3d}{4D}+\dfrac{\partial\delta'_1}{\partial X}\right)$", EQ),
 (r"$\dfrac{\partial\ell_{\min,2}}{\partial X}=-\dfrac{180^\circ}{\pi}\left(\dfrac{\partial}{\partial X}\dfrac{3d}{4D}+\dfrac{\partial\delta'_2}{\partial X}\right)$", EQ),
 (r"$\dfrac{\partial}{\partial d}\dfrac{3d}{4D}=\dfrac{3}{4D},\qquad\dfrac{\partial}{\partial D}\dfrac{3d}{4D}=-\dfrac{3d}{4D^{2}},\qquad\dfrac{\partial}{\partial L}\dfrac{3d}{4D}=\dfrac{\partial}{\partial M}\dfrac{3d}{4D}=0$", EQ),
 (r"$\partial\delta_{1,2}$ is on pages 3–5 and $\partial\delta'_{1,2}$ on page 6.", NT),
]
P[8] = [
 (r"$\bf{Page\ 8:\ heights\ of\ the\ maxima}$", HD), (r"(upper sign max,1, lower sign max,2; $X$ = $D$, $L$, $d$ or $M$)", NT), ("gap", 0.1),
 (r"$\dfrac{\partial}{\partial X}\left(\dfrac{A}{A_{\rm all}}\right)_{\max,1,2}=\left(\dfrac{A}{A_{\rm all}}\right)_{\max,1,2}\times\left[\dfrac{-\frac{1}{2}\mathrm{Amp}_A\cos\phi_A\,\partial\delta_{1,2}/\partial X}{1\mp\mathrm{Amp}_A\sin\phi_A-\frac{1}{2}\mathrm{Amp}_A\cos\phi_A\,\delta_{1,2}}\right.$", 17),
 (r"$\qquad\qquad\left.+\dfrac{\partial(H_{\rm obs}/H)_{\max,1,2}/\partial X}{(H_{\rm obs}/H)_{\max,1,2}}-\dfrac{\partial A_{\rm all}/\partial X}{A_{\rm all}}\right]$", EQ),
 ("gap", 0.15),
 (B + r"$(H_{\rm obs}/H)_{\max,1,2}$ is $H_{\rm obs}/H$ at $\ell=0^\circ$ / $180^\circ$: the curly-brace factor", NT),
 (r"$\qquad$$\left\{1+\frac{\alpha GM}{H^{2}}\left[\ldots\right]\right\}$ in the card's height formula.", NT),
 (B + r"For $X=D,\,L,\,d$:  $\partial(H_{\rm obs}/H)_{\max,1,2}/\partial X$ equals $\alpha GM/H^{2}$ times the square", NT),
 (r"$\qquad$bracket that multiplies Amp$_A\cos\phi_A$ in $\partial\delta_{1,2}/\partial X$ (pages 3–5).", NT),
 (r"$\qquad$(Exact, because $(H_{\rm obs}/H)_{\max}$ is $\alpha GM/H^{2}$ times the numerator brace of $\delta$.)", NT),
 (B + r"For $X=M$:", NT),
 (r"$\dfrac{\partial(H_{\rm obs}/H)_{\max,1,2}}{\partial M}=\dfrac{(H_{\rm obs}/H)_{\max,1,2}-1}{M}$", EQ),
 (B + r"$\partial A_{\rm all}/\partial X$ is on page 2.", NT),
]
P[9] = [
 (r"$\bf{Page\ 9:\ heights\ of\ the\ minima}$", HD), (r"(upper sign min,1, lower sign min,2)", NT), ("gap", 0.1),
 (r"$\dfrac{\partial}{\partial X}\left(\dfrac{A}{A_{\rm all}}\right)_{\min,1,2}=\left(\dfrac{A}{A_{\rm all}}\right)_{\min,1,2}$", EQ),
 ("\x01" + r"$\times\left[\mathrm{Amp}_A\sin\phi_A\left(\dfrac{\partial}{\partial X}\dfrac{3d}{4D}+\dfrac{1}{2}\dfrac{\partial\delta'_{1,2}}{\partial X}\right)\right.$", EQ),
 ("\x01\x01" + r"$\times\left(1\mp\mathrm{Amp}_A\cos\phi_A+\dfrac{3d}{4D}\mathrm{Amp}_A\sin\phi_A+\dfrac{1}{2}\mathrm{Amp}_A\sin\phi_A\,\delta'_{1,2}\right)^{-1}$", EQ),
 ("\x01" + r"$\left.+\dfrac{\partial(H_{\rm obs}/H)_{\min}/\partial X}{1-\alpha GM/D^{3}H^{2}}-\dfrac{\partial A_{\rm all}/\partial X}{A_{\rm all}}\right]$", EQ),
 ("gap", 0.15),
 (r"with $(H_{\rm obs}/H)_{\min}=1-\alpha GM/D^{3}H^{2}$, so", NT),
 (r"$\dfrac{\partial(H_{\rm obs}/H)_{\min}}{\partial D}=\dfrac{3\alpha GM}{D^{4}H^{2}},\qquad\dfrac{\partial(H_{\rm obs}/H)_{\min}}{\partial M}=-\dfrac{\alpha G}{D^{3}H^{2}},\qquad\dfrac{\partial(H_{\rm obs}/H)_{\min}}{\partial L}=\dfrac{\partial(H_{\rm obs}/H)_{\min}}{\partial d}=0$", 16.2),
 (B + r"$\partial(3d/4D)/\partial X$ is on page 7, $\partial\delta'_{1,2}/\partial X$ on page 6, $\partial A_{\rm all}/\partial X$ on page 2.", NT),
]
# ---------- baseline tables ----------
ORDER = ["lmax1", "lmax2", "lmin1", "lmin2", "Amax1", "Amax2", "Amin1", "Amin2"]
LAB = {"lmax1": r"$\ell_{\max,1}$", "lmax2": r"$\ell_{\max,2}$", "lmin1": r"$\ell_{\min,1}$", "lmin2": r"$\ell_{\min,2}$",
       "Amax1": r"$A_{\max,1}$", "Amax2": r"$A_{\max,2}$", "Amin1": r"$A_{\min,1}$", "Amin2": r"$A_{\min,2}$"}
INP = [("D", r"$D$"), ("L", r"$L$"), ("d", r"$d$"), ("M", r"$M$"), ("alpha", r"$\alpha$")]
base = {k: float(v) for k, v in J["base"].items()}
def fmt(v): return ("%+.3g" % v).replace("-", "\u2212")
def table(names):
    rows = [("", "", "analytic", "exact FD", "ratio", "elast. an.", "elast. ex.")]
    for nm in names:
        for i, (k, kl) in enumerate(INP):
            an = J["ana"][nm][k]; ex = J["efd"][k][nm]
            ea = base[k]*an/J["val"][nm]; ee = base[k]*ex/J["exact"][nm]
            rows.append((LAB[nm] if i == 0 else "", kl, fmt(an), fmt(ex), "%.2f" % (an/ex), fmt(ea), fmt(ee)))
    return ("table", rows, [0, 1.1, 1.75, 3.3, 4.85, 5.95, 7.6], 15)
P[10] = [(r"$\bf{Page\ 10:\ baseline\ derivatives\ of\ the\ angles}$", HD),
         (r"Degrees per unit of input. Baseline $D=1$, $L=D/5.9$, $d=0.1$, $M=0.2306$, $\alpha=1$ ($G=H=1$).", NT),
         ("gap", 0.05), table(ORDER[:4])]
P[11] = [(r"$\bf{Page\ 11:\ baseline\ derivatives\ of\ the\ heights}$", HD),
         (r"$A$ means $A/A_{\rm all}$; change per unit of input; same baseline.", NT), ("gap", 0.05), table(ORDER[4:])]
worst = sorted(((abs(J["ana"][nm][k]/J["efd"][k][nm] - 1), nm, k) for nm in ORDER for k, _ in INP), reverse=True)
P[12] = [
 (r"$\bf{Page\ 12:\ notes}$", HD), ("gap", 0.1),
 (B + r"analytic = derivative of the card formula. It matches a central finite difference", NT),
 (r"$\qquad$of the card formulas to $3\times10^{-9}$ (relative) for all 40 entries.", NT),
 (B + r"exact FD = central finite difference (step $10^{-3}$ of each input) of the true extrema of the", NT),
 (r"$\qquad$no-boxcar curve with the tilt; a $3\times10^{-4}$ step changes these by under 0.1%.", NT),
 (B + r"ratio = analytic / exact FD.  elast. = (input / quantity) $\times$ derivative:", NT),
 (r"$\qquad$the % change of the quantity per 1% change of the input. For angles it is relative", NT),
 (r"$\qquad$to the angle itself, so compare angle elasticities among angles only.", NT),
 (B + r"Each quantity's elasticities obey $e_D+e_L+e_d+3e_M=0$ and $e_\alpha=e_M$ (page 1 shortcuts).", NT),
 (B + r"Heights: the first-order derivatives track the exact ones to within 4% for $D$, $L$, $M$, $\alpha$.", NT),
 (r"$\qquad$The $d$-derivatives of the heights are small and off by 13–24%.", NT),
 (B + r"Poorly tracked: $\partial\ell_{\max,1}/\partial d$ is 1.8$\times$ too large and $\partial\ell_{\max,2}/\partial d$ 0.70$\times$;", NT),
 (r"$\qquad\partial\ell_{\min,1}/\partial D$ is 1.5$\times$ (a small number, $-1.1$ vs $-0.73$ deg per unit $D$).", NT),
 (r"$\qquad\ell_{\max,1}$'s other derivatives run 16–21% high, as $\ell_{\max,1}$ itself does (8.9$^\circ$ vs 7.8$^\circ$).", NT),
 (r"$\qquad\ell_{\min,2}$: $D$ 18% low; $L$, $M$, $\alpha$ 15% high.  $\ell_{\min,1}$: $L$, $M$, $\alpha$ 10% low.", NT),
 (B + r"$\ell_{\max,2}$ tracks within 7% for $D$, $L$, $M$, $\alpha$.", NT),
 (B + r"$\partial\ell_{\min}/\partial d$ is good (within 0.5%): the $135^\circ/(\pi D)$ term dominates it.", NT),
]
for n, items in P.items():
    render(n, items)
print("worst ratios:", [(round(r, 3), nm, k) for r, nm, k in worst[:8]])
