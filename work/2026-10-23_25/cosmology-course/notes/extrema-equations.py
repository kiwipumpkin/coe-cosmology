#!/usr/bin/env python3
"""
Equations for the extrema of the model A/A_all(l) curve, plus implicit-function
sensitivities, fold (merge) points and a small-eps closed form.

Model = af_an_alpha / scale_kappa / model_AB_fixed from
  /workspace/figs_letter_exact_python/makers/make_amp_phi_and_geometry_heatmaps.py
  /workspace/cosmology/density_local_match/run_density_local_match.py
re-written below as explicit formulas (h, m, m', boxcar as an integral) and
checked against the imported code at baseline.  NO fitting anywhere: every
number is an exact root, an implicit-function derivative, or a series term.

Outputs (in /workspace/cosmology-course/notes):
  extrema-equations-p*.png, extrema-equations.md, extrema-equations.json
"""
from __future__ import annotations
import importlib.util, json, re
from pathlib import Path
import numpy as np
from scipy.optimize import brentq, fsolve
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

NOTES = Path("/workspace/cosmology-course/notes")
FU, AMP, PHI = 0.94, -0.06, 0.35
K0 = 0.23058828928336284
W = np.radians(91 * 360 / 361)      # code: 91-sample boxcar on 361-pt grid
w = W / 2
KEYS = ["x", "e", "a", "k"]
BASE = dict(x=5.9, e=0.1, a=1.0, k=K0)
BASE13 = dict(x=5.9, e=0.1, a=1.3, k=K0)   # Fig A13 hold (code hold["alpha"]=1.3)
GLx, GLw = np.polynomial.legendre.leggauss(160)
TG = np.linspace(0, 2 * np.pi, 8192, endpoint=False)

# ----------------------------------------------------------------- model
def N(x): return 6 * (1 - FU) * (np.exp(x) - 1 - x - x**2 / 2) + FU * x**3
def m(r, x, k):
    return k * (6 * (1 - FU) * (np.exp(x) - np.exp(x * (1 - r)) * (1 + x * r + x**2 * r**2 / 2))
                + FU * x**3 * r**3) / N(x)
def h(t, x, e, a, k):
    c = np.cos(t); r = np.sqrt(1 + 2 * e * c + e * e)
    return 1 + (a / e) * (k * c - m(r, x, k) * (c + e) / r**3)
def hp(t, x, e, a, k):   # fully simplified dh/dtheta (floor term cancels in r m' - 3m)
    c, s = np.cos(t), np.sin(t); r = np.sqrt(1 + 2 * e * c + e * e); n = N(x)
    T1 = (6 * (1 - FU) * (np.exp(x) - np.exp(x * (1 - r)) * (1 + x * r + x**2 * r**2 / 2)) + FU * x**3 * r**3) / (r**3 * n)
    T2 = 18 * (1 - FU) * e * (c + e) * (np.exp(x) - np.exp(x * (1 - r)) * (1 + x * r + x**2 * r**2 / 2 + x**3 * r**3 / 6)) / (r**5 * n)
    return (a * k * s / e) * (T1 - 1 - T2)
def P(d): return tuple(d[q] for q in KEYS)
def S(v): return float(np.ravel(v)[0])
def hbar(x, e, a, k): return h(TG, x, e, a, k).mean()
def AF(t, x, e, a, k):
    t = np.atleast_1d(t)[:, None]; tt = t + w * GLx[None, :]
    return (h(tt, x, e, a, k) * GLw).sum(1) * w / (W * hbar(x, e, a, k))
def A(t, p): return S((1 - AMP * np.sin(t + PHI)) * AF(t, *p))
def At(t, p):
    hb = hbar(*p); afp = (h(t + w, *p) - h(t - w, *p)) / (W * hb)
    return S(-AMP * np.cos(t + PHI) * AF(t, *p) + (1 - AMP * np.sin(t + PHI)) * afp)
def Att(t, p):
    hb = hbar(*p); af = AF(t, *p)
    afp = (h(t + w, *p) - h(t - w, *p)) / (W * hb); afpp = (hp(t + w, *p) - hp(t - w, *p)) / (W * hb)
    return S(AMP * np.sin(t + PHI) * af - 2 * AMP * np.cos(t + PHI) * afp + (1 - AMP * np.sin(t + PHI)) * afpp)
def roots(p, n=1440):
    t = np.linspace(0, 2 * np.pi, n, endpoint=False); g = np.array([At(s, p) for s in t]); out = []
    for i in range(n):
        j = (i + 1) % n
        if g[i] * g[j] < 0:
            r = brentq(lambda s: At(s, p), t[i], t[i] + 2 * np.pi / n, xtol=1e-13) % (2 * np.pi)
            out.append(dict(t=r, l=np.degrees(r), kind="max" if Att(r, p) < 0 else "min", A=A(r, p)))
    return sorted(out, key=lambda d: d["l"])

# --------------------------------------------- implicit-function derivatives
def implicit(t0, base, q):
    dq = 1e-4 * max(1, abs(base[q])); dt = 1e-4
    def F(t, d): return At(t, P({**base, q: base[q] + d}))
    def Av(t, d): return A(t, P({**base, q: base[q] + d}))
    Ft = Att(t0, P(base))
    Fp = (F(t0, dq) - F(t0, -dq)) / (2 * dq)
    Fpp = (F(t0, dq) - 2 * F(t0, 0) + F(t0, -dq)) / dq**2
    Ftp = (Att(t0, P({**base, q: base[q] + dq})) - Att(t0, P({**base, q: base[q] - dq}))) / (2 * dq)
    Ftt = (Att(t0 + dt, P(base)) - Att(t0 - dt, P(base))) / (2 * dt)
    Ap = (Av(t0, dq) - Av(t0, -dq)) / (2 * dq)
    App = (Av(t0, dq) - 2 * Av(t0, 0) + Av(t0, -dq)) / dq**2
    t1 = -Fp / Ft
    t2 = -(Fpp + 2 * Ftp * t1 + Ftt * t1**2) / Ft
    return dict(dl=np.degrees(t1), d2l=np.degrees(t2), dA=Ap, d2A=App + Fp * t1)

# ------------------------------------------------- small-eps closed form
sk = lambda k: np.sin(k * w) / (k * w)
def harmonics(x, e, a, k):
    n = N(x)
    h0 = 1 - (a * k / 2) * (3 * x**3 / n - 1)
    h1 = -(9 * a * e * k / 8) * (1 - x**3 * (1 + x * (1 - FU)) / n)
    h2 = -(a * k / 2) * (3 * x**3 / n - 3)
    h3 = -(3 * a * e * k / 8) * (5 - x**3 * (5 + x * (1 - FU)) / n)
    return h0, (h1, h2, h3), tuple(sk(j + 1) * hh / h0 for j, hh in enumerate((h1, h2, h3)))
def closed_O_eps(x, e, a, k):
    _, _, (a1, a2, a3) = harmonics(x, e, a, k)
    AFf = lambda t: 1 + a1 * np.cos(t) + a2 * np.cos(2 * t) + a3 * np.cos(3 * t)
    AFpp = lambda t: -a1 * np.cos(t) - 4 * a2 * np.cos(2 * t) - 9 * a3 * np.cos(3 * t)
    disc = a2**2 - 3 * a3 * (a1 - 3 * a3)
    cs = [(-a2 + sg * np.sqrt(disc)) / (6 * a3) for sg in (1, -1)] if disc >= 0 else []
    cs = [c for c in cs if abs(c) <= 1]
    t0s = [0.0, np.pi] + [np.arccos(c) for c in cs] + [2 * np.pi - np.arccos(c) for c in cs]
    return [dict(l=np.degrees(t0 + AMP * np.cos(t0 + PHI) * AFf(t0) / AFpp(t0)) % 360,
                 kind="max" if AFpp(t0) < 0 else "min",
                 A=(1 - AMP * np.sin(t0 + PHI)) * AFf(t0)) for t0 in t0s]
def closed_eps0(x, e, a, k):
    _, _, (_, a2, _) = harmonics(x, e, a, k)
    d1 = np.degrees(AMP * np.cos(PHI) * (1 + a2) / (4 * a2))
    d2 = np.degrees(AMP * np.sin(PHI) * (1 - a2) / (4 * a2))
    return [dict(l=(-d1) % 360, kind="max", A=(1 - AMP * np.sin(PHI)) * (1 + a2)),
            dict(l=90 - d2, kind="min", A=(1 - AMP * np.cos(PHI)) * (1 - a2)),
            dict(l=180 + d1, kind="max", A=(1 + AMP * np.sin(PHI)) * (1 + a2)),
            dict(l=270 + d2, kind="min", A=(1 + AMP * np.cos(PHI)) * (1 - a2))]

def cdist(a, b): return abs((a - b + 180) % 360 - 180)

# ------------------------------------------------------------- compute all
def compute():
    out = {}
    # 0. faithfulness check vs imported code
    spec = importlib.util.spec_from_file_location("lv", "/workspace/figs_letter_exact_python/makers/make_one_at_a_time_levelsets.py")
    lv = importlib.util.module_from_spec(spec); spec.loader.exec_module(lv)
    lf, ma = lv.model_A_curve(5.9, 0.94, 1.0, 0.1, K0)
    mine = np.array([A(np.radians(l), P(BASE)) for l in lf])
    out["check_vs_code_maxabs"] = float(np.max(np.abs(mine - ma)))
    tchk = np.linspace(0, 2 * np.pi, 13); d = 1e-6
    out["check_hprime_maxabs"] = float(np.max(np.abs((h(tchk + d, *P(BASE)) - h(tchk - d, *P(BASE))) / (2 * d) - hp(tchk, *P(BASE)))))
    # 1. baseline extrema + sensitivities
    for name, base in (("base", BASE), ("base13", BASE13)):
        ex = roots(P(base))
        for k_, e in enumerate(ex):
            e["k"] = k_ + 1
            e["sens"] = {q: implicit(e["t"], base, q) for q in KEYS}
            hb = hbar(*P(base)); t0 = e["t"]
            e["A_from_eq3_2"] = S((1 - AMP * np.sin(t0 + PHI))**2 * (h(t0 + w, *P(base)) - h(t0 - w, *P(base))) / (AMP * np.cos(t0 + PHI) * W * hb))
        out[name] = ex
    out["hbar_base"] = hbar(*P(BASE)); out["hbar_base13"] = hbar(*P(BASE13))
    out["harmonics_base"] = harmonics(*P(BASE)); out["harmonics_base13"] = harmonics(*P(BASE13))
    # 2. folds
    folds = []
    for q, gl, gv, bb_ in (("x", 127, 3.16, BASE), ("a", 127, 0.31, BASE), ("k", 127, 0.07, BASE),
                           ("e", 135, 0.964, BASE13), ("e", 204, 0.967, BASE13)):
        def fn(z):
            b = P({**bb_, q: z[1]}); return [At(z[0], b), Att(z[0], b)]
        sol, info, ier, msg = fsolve(fn, [np.radians(gl), gv], full_output=True, xtol=1e-13)
        b = {**bb_, q: sol[1]}; tf = sol[0]
        dq = 1e-5; dt = 1e-5
        Fp = (At(tf, P({**b, q: sol[1] + dq})) - At(tf, P({**b, q: sol[1] - dq}))) / (2 * dq)
        Ftt = (Att(tf + dt, P(b)) - Att(tf - dt, P(b))) / (2 * dt)
        folds.append(dict(var=q, value=float(sol[1]), l=float(np.degrees(tf) % 360), A=A(tf, P(b)),
                          resid=[float(v) for v in fn(sol)], Fp=Fp, Ftt=Ftt,
                          puiseux_deg=float(np.degrees(np.sqrt(abs(2 * Fp / Ftt)))),
                          sign=float(np.sign(-2 * Fp / Ftt))))
    out["folds"] = folds
    # 3. accuracy against the exact (code) sweeps
    spec = importlib.util.spec_from_file_location("ev", NOTES / "extrema-vs-inputs.py")
    ev = importlib.util.module_from_spec(spec); spec.loader.exec_module(ev)
    SW = {"x": ("x", 3, 10, 5.9, {}, "base"), "e": ("eps", .06, .74, .1, {"alpha": 1.3}, "base13"),
          "a": ("alpha", .2, 2, 1.0, {}, "base"), "k": ("kappa", .15, .7, K0, {}, "base")}
    acc = {}
    for q, (key, lo, hi, b0, over, bn) in SW.items():
        T = ev.track(key, lo, hi, b0, over); v = T["vals"]
        for br in T["branches"]:
            e = min([e for e in out[bn] if e["kind"] == br["kind"]], key=lambda e: cdist(e["l"], br["l0"]))
            s = e["sens"][q]; dv = v - b0; ok = ~np.isnan(br["l"])
            L1 = e["l"] + s["dl"] * dv; L2 = L1 + .5 * s["d2l"] * dv**2
            A1 = e["A"] + s["dA"] * dv; A2 = A1 + .5 * s["d2A"] * dv**2
            el1 = np.abs(L1 - br["l"])[ok]; el2 = np.abs(L2 - br["l"])[ok]
            ea1 = np.abs(A1 - br["A"])[ok]; ea2 = np.abs(A2 - br["A"])[ok]
            ec, ec0 = [], []
            for i in np.where(ok)[0]:
                kw = {**T["hold"], key: float(v[i])}
                pp = (kw["x"], kw["eps"], kw["alpha"], kw["kappa"])
                for fun, lst in ((closed_O_eps, ec), (closed_eps0, ec0)):
                    C = [c for c in fun(*pp) if c["kind"] == br["kind"]]
                    lst.append(min(cdist(c["l"], br["l"][i] % 360) for c in C) if C else np.inf)
            vv = v[ok]
            def within(err, tol=1.0):
                g = vv[err <= tol]
                return [float(g.min()), float(g.max())] if len(g) else None
            # error at the sweep endpoints (as listed in extrema-vs-inputs.json)
            acc[f"{q}|k{e['k']}"] = dict(sweep=q, k=e["k"], kind=br["kind"], l0=e["l"],
                range=[float(vv[0]), float(vv[-1])],
                lin_l=float(el1.max()), quad_l=float(el2.max()), lin_A=float(ea1.max()), quad_A=float(ea2.max()),
                cOe_l=float(np.max(ec)), c0_l=float(np.max(ec0)),
                lin_l_ends=[float(el1[0]), float(el1[-1])], quad_l_ends=[float(el2[0]), float(el2[-1])],
                lin_1deg=within(el1), quad_1deg=within(el2), cOe_2deg=within(np.array(ec), 2.0),
                c0_2deg=within(np.array(ec0), 2.0))
    out["accuracy"] = acc
    out["closed_base"] = dict(Oe=closed_O_eps(*P(BASE)), e0=closed_eps0(*P(BASE)),
                              Oe13=closed_O_eps(*P(BASE13)), e013=closed_eps0(*P(BASE13)))
    return out

# ================================================================ rendering
PW, PH, DPI = 8.5, 11.0, 200
ML, MR, MT, MB = 0.5, 0.5, 0.55, 0.6
CW = PW - ML - MR
FS_P, FS_EQ, FS_H, FS_T = 13.5, 15.5, 18, 12.5

OVER = []
def tok(s): return re.findall(r"(?:\$[^$]+\$|[^\s$])+", s)

class Pager:
    def __init__(s): s.pages = []; s.new()
    def new(s):
        s.fig = plt.figure(figsize=(PW, PH), dpi=DPI); s.y = MT; s.pages.append(s.fig)
        s.r = s.fig.canvas.get_renderer()
    def fx(s, xin): return xin / PW
    def fy(s, yin): return 1 - yin / PH
    def txt(s, xin, yin, string, **kw):
        t = s.fig.text(s.fx(xin), s.fy(yin), string, va="top", **kw)
        bb = t.get_window_extent(s.r); return t, bb.width / DPI, bb.height / DPI
    def room(s, hin):
        if s.y + hin > PH - MB: s.new()
    def heading(s, string):
        s.room(1.6); s.y += 0.12
        t, wd, ht = s.txt(ML, s.y, string, fontsize=FS_H, fontweight="bold")
        s.y += ht + 0.06
        s.fig.add_artist(Line2D([s.fx(ML), s.fx(PW - MR)], [s.fy(s.y)] * 2, color="k", lw=0.8))
        s.y += 0.12
    def para(s, string, fs=FS_P, indent=0.0):
        words = tok(string); lines = []; cur = ""
        for wd_ in words:
            trial = (cur + " " + wd_).strip()
            t = s.fig.text(0, 0, trial, fontsize=fs); wpx = t.get_window_extent(s.r).width / DPI; t.remove()
            if wpx > CW - indent and cur:
                lines.append(cur); cur = wd_
            else:
                cur = trial
        if cur: lines.append(cur)
        for ln in lines:
            t = s.fig.text(0, 0, ln, fontsize=fs); ht = t.get_window_extent(s.r).height / DPI; t.remove()
            ht = max(ht, fs / 72 * 1.25)
            s.room(ht + 0.05)
            t, wd, _ = s.txt(ML + indent, s.y, ln, fontsize=fs)
            assert wd <= CW + 1e-3, ("para overflow", ln)
            s.y += ht + 0.07
        s.y += 0.08
    def eq(s, lines, label=None):
        hs = []
        for ln in lines:
            t = s.fig.text(0, 0, ln, fontsize=FS_EQ); bb = t.get_window_extent(s.r); t.remove()
            if bb.width / DPI > CW - 0.75 - (0.35 if hs else 0): OVER.append(f"({label}) {bb.width/DPI:.2f} in: {ln[:70]}")
            hs.append(bb.height / DPI)
        s.room(sum(hs) + 0.14 * len(lines) + 0.1)
        s.y += 0.06
        for i, ln in enumerate(lines):
            cont = i > 0 and re.match(r"^\$(=|\\times|\+|-|\\left\.)", ln) is not None
            t, wd, ht = s.txt(ML + 0.3 + (0.35 if cont else 0), s.y, ln, fontsize=FS_EQ)
            if label and i == len(lines) - 1:
                s.fig.text(s.fx(PW - MR), s.fy(s.y + ht / 2), f"({label})", va="center", ha="right", fontsize=FS_P)
            s.y += ht + 0.14
        s.y += 0.06
    def table(s, header, rows, widths, title=None, fs=FS_T):
        assert sum(widths) <= CW + 1e-6, ("table too wide", sum(widths))
        rh = fs / 72 * 1.75
        if title:
            s.room(0.5 + rh * min(len(rows) + 1, 5)); s.para(title)
        def line(yin, lw=0.8):
            s.fig.add_artist(Line2D([s.fx(ML), s.fx(ML + sum(widths))], [s.fy(yin)] * 2, color="k", lw=lw))
        allrows = [header] + rows
        def allrows_hdr(hdr):
            xx = ML
            for ci, (cell, wd) in enumerate(zip(hdr, widths)):
                s.fig.text(s.fx(xx + 0.05 if ci == 0 else xx + wd / 2), s.fy(s.y + rh / 2), cell, va="center",
                           ha="left" if ci == 0 else "center", fontsize=fs, fontweight="bold")
                xx += wd
            s.y += rh; line(s.y, 0.8)
        s.room(rh * min(len(allrows), 4) + 0.2)
        line(s.y, 1.1)
        for ri, row in enumerate(allrows):
            if s.y + rh > PH - MB:
                s.new(); line(s.y, 1.1); allrows_hdr(header)
            xx = ML
            for ci, (cell, wd) in enumerate(zip(row, widths)):
                ha = "left" if ci == 0 else "center"
                xpos = xx + 0.05 if ci == 0 else xx + wd / 2
                t = s.fig.text(s.fx(xpos), s.fy(s.y + rh / 2), cell, va="center", ha=ha, fontsize=fs,
                               fontweight="bold" if ri == 0 else "normal")
                cw_ = t.get_window_extent(s.r).width / DPI
                if cw_ > wd - 0.04: OVER.append(f"cell {cell!r} {cw_:.2f}>{wd}")
                xx += wd
            s.y += rh
            if ri == 0: line(s.y, 0.8)
        line(s.y, 1.1); s.y += 0.2
    def gap(s, g): s.y += g
    def finish(s, stem):
        n = len(s.pages); paths = []
        for i, f in enumerate(s.pages):
            f.text(0.5, 0.5 * MB / PH, f"Extremum equations for the A/A_all model — page {i+1} of {n}",
                   ha="center", va="center", fontsize=10, color="#555555")
            p = NOTES / f"{stem}-p{i+1}.png"; f.savefig(p, dpi=DPI); plt.close(f); paths.append(p)
        return paths

def f3(v, n=3): return f"{v:+.{n}f}".replace("-", "−")
def g4(v): return f"{v:.4g}".replace("-", "−")

def build(R):
    doc = []  # (kind, payload)
    H = lambda s: doc.append(("h", s)); Pp = lambda s: doc.append(("p", s))
    E = lambda lines, lab=None: doc.append(("eq", (lines, lab)))
    T = lambda hdr, rows, wid, title=None: doc.append(("t", (hdr, rows, wid, title)))
    b = R["base"]; b13 = R["base13"]; fx, fa, fk, fe, fe2 = R["folds"]
    Pp(r"$\mathbf{Extrema\ of\ the\ model\ } A/A_{\rm all}(\ell)\mathbf{:\ equations,\ sensitivities\ and\ merge\ points}$")
    Pp("This sheet writes the model behind Figs A10–A13 explicitly in its inputs, derives the condition that fixes the "
       "Galactic longitude $\\ell$ of each maximum and minimum, gives the value of $A/A_{\\rm all}$ there, the condition "
       "under which a maximum and a minimum merge, and the first- and second-order dependence of all four extrema on "
       "the four geometry inputs. Nothing here is fitted: every number is an exact root of the model equations, an "
       "implicit-function derivative, or a term of a series expansion.")
    # ---------------- section 1
    H("1. Model and AF written in the inputs")
    Pp("Symbols. $\\ell$: Galactic longitude [deg]; $\\theta=\\pi\\ell/180$ [rad]. $D$: distance from the center of "
       "expansion (CoE) to the observer. $L$: scale length of the exponential core of the density profile. "
       "$x=D/L$. $d$: radius of the sampled shell of sources around the observer; $\\varepsilon=d/D$. "
       "$\\alpha$: dimensionless scale of the Gauss peculiar-velocity ansatz. "
       "$\\kappa=GM(<D)/(D^3H^2)$: dimensionless mass inside the observer's radius. $H$: background (homologous) "
       "expansion rate; $G$: Newton's constant. $\\mathrm{Amp}_A=-0.06$ and $\\phi_A=0.35$ rad: amplitude and phase "
       "of the sine modulation. $f_u=0.94$: fixed share of the local density $\\rho(D)$ carried by the uniform floor. "
       "$W$: width of the $\\Delta\\ell\\approx90^\\circ$ boxcar smoothing. Units: $G=H=D=1$, so lengths are in $D$, "
       "speeds in $HD$, masses in $D^3H^2/G$.")
    E([r"$A(\ell)=\left[1-\mathrm{Amp}_A\,\sin\left(\theta+\phi_A\right)\right]\,AF(\theta),\qquad\theta=\dfrac{\pi\,\ell}{180}$"], "1.1")
    Pp("Geometry. The CoE sits at the origin and the observer at $\\vec D=(1,0,0)$. A line of sight has unit vector "
       "$\\hat n=(\\cos\\theta,\\sin\\theta,0)$, so $\\theta=0$ points directly away from the CoE and $\\theta=\\pi$ "
       "toward it (the code identifies this angle with $\\ell$ directly; the axis offset is carried by $\\phi_A$). "
       "The source on the shell in that direction is at")
    E([r"$\vec r=\vec D+\varepsilon\,\hat n=\left(1+\varepsilon\cos\theta,\ \varepsilon\sin\theta,\ 0\right),\qquad r(\theta)=\sqrt{1+2\,\varepsilon\cos\theta+\varepsilon^2}$"], "1.2")
    Pp("Mass. The density is an exponential core plus a uniform floor, $\\rho(r)=\\rho_{c0}\\,\\mathrm{e}^{-xr}+\\rho_u$. "
       "The code first sets $\\rho(1)=\\Omega_m\\rho_{\\rm crit}$ with the floor holding the share $f_u$, then rescales "
       "every mass so that $m(1)=\\kappa$. Integrating $4\\pi r^2\\rho$ and dividing by its value at $r=1$, the matter "
       "density $\\Omega_m$, $\\rho_{\\rm crit}$, $G$ and $H$ all cancel, leaving")
    E([r"$N(x)=6\,(1-f_u)\left(\mathrm{e}^{x}-1-x-\dfrac{x^2}{2}\right)+f_u\,x^3$"], "1.3")
    E([r"$m(r)=\dfrac{\kappa}{N(x)}\left\{6\,(1-f_u)\left[\mathrm{e}^{x}-\mathrm{e}^{x(1-r)}\left(1+x\,r+\dfrac{x^2r^2}{2}\right)\right]+f_u\,x^3r^3\right\}$"], "1.4")
    Pp("Velocity. Gauss peculiar-velocity ansatz (inward pull of the enclosed mass), in the units above:")
    E([r"$\vec u(\vec r)=-\alpha\,\dfrac{G\,m(r)}{r^3H}\,\vec r=-\alpha\,\dfrac{m(r)}{r^3}\,\vec r$"], "1.5")
    Pp("Observed expansion rate along $\\hat n$ (the differential radial velocity across the shell, divided by its "
       "radius). With $\\vec u(\\vec r)\\cdot\\hat n=-\\alpha\\,m(r)(\\cos\\theta+\\varepsilon)/r^3$ and "
       "$\\vec u(\\vec D)\\cdot\\hat n=-\\alpha\\kappa\\cos\\theta$:")
    E([r"$h(\theta)\equiv\dfrac{H_{\rm obs}(\theta)}{H}=1+\dfrac{\left[\vec u(\vec r)-\vec u(\vec D)\right]\cdot\hat n}{\varepsilon}$",
       r"$=1+\dfrac{\alpha}{\varepsilon}\left[\kappa\cos\theta-\dfrac{m(r)\,(\cos\theta+\varepsilon)}{r^3}\right]$"], "1.6")
    Pp("Substituting (1.4), fully in the inputs:")
    E([r"$h(\theta)=1+\dfrac{\alpha\kappa}{\varepsilon}\cos\theta-\dfrac{\alpha\kappa\,(\cos\theta+\varepsilon)}{\varepsilon\,r^3\,N(x)}$",
       r"$\times\left\{6\,(1-f_u)\left[\mathrm{e}^{x}-\mathrm{e}^{x(1-r)}\left(1+x\,r+\dfrac{x^2r^2}{2}\right)\right]+f_u\,x^3r^3\right\}$"], "1.7")
    Pp("Two consequences. (i) $\\alpha$ and $\\kappa$ enter only through the product $\\alpha\\kappa$, so the "
       "$\\alpha$ and $\\kappa$ sweeps trace the same family of curves. (ii) The code's AF mode \"hratio\" uses "
       "$A_{\\rm raw}=H_{\\rm obs}/H=h$; its floor $H_{\\rm obs}\\geq10^{-6}H$ never acts on the swept ranges "
       f"(smallest $h$ there is 0.30).")
    Pp("Smoothing and normalisation. The code divides $h$ by its sky mean and applies a wrap-around boxcar of 91 "
       "samples on its 361-point grid, i.e. width $W=91\\times360^\\circ/361=90.75^\\circ$. As an integral:")
    E([r"$AF(\theta)=\dfrac{1}{W\,\bar h}\int_{\theta-W/2}^{\theta+W/2}h(\theta^\prime)\,d\theta^\prime,\qquad"
       r"\bar h=\dfrac{1}{2\pi}\int_{0}^{2\pi}h(\theta)\,d\theta,\qquad W=\dfrac{182\,\pi}{361}$"], "1.8")
    Pp("The core term contains $\\mathrm{e}^{-x\\,r(\\theta)}$ with $r(\\theta)=\\sqrt{1+2\\varepsilon\\cos\\theta+\\varepsilon^2}$, "
       "which has no elementary antiderivative in $\\theta$ (only the floor term integrates in closed form). "
       "So $AF$ itself stays a one-dimensional integral. Its slope and curvature do not: by the fundamental theorem "
       "of calculus they are differences of $h$ and $h^\\prime$ at the two window edges.")
    E([r"$\dfrac{\partial AF}{\partial\theta}=\dfrac{h\left(\theta+\frac{W}{2}\right)-h\left(\theta-\frac{W}{2}\right)}{W\,\bar h},\qquad"
       r"\dfrac{\partial^2 AF}{\partial\theta^2}=\dfrac{h^\prime\left(\theta+\frac{W}{2}\right)-h^\prime\left(\theta-\frac{W}{2}\right)}{W\,\bar h}$"], "1.9")
    Pp("Differentiating (1.7), using $dr/d\\theta=-\\varepsilon\\sin\\theta/r$ and "
       "$dm/dr=3\\kappa x^3r^2\\left[(1-f_u)\\mathrm{e}^{x(1-r)}+f_u\\right]/N(x)$; in $r\\,dm/dr-3m$ the floor "
       "cancels exactly, leaving:")
    E([r"$h^\prime(\theta)=\dfrac{\alpha\kappa\sin\theta}{\varepsilon}\left\{\dfrac{6\,(1-f_u)\left[\mathrm{e}^{x}-\mathrm{e}^{x(1-r)}\left(1+x\,r+\frac{x^2r^2}{2}\right)\right]+f_u\,x^3r^3}{r^3\,N(x)}-1\right.$",
       r"$\left.-\dfrac{18\,(1-f_u)\,\varepsilon\,(\cos\theta+\varepsilon)}{r^5\,N(x)}\left[\mathrm{e}^{x}-\mathrm{e}^{x(1-r)}\left(1+x\,r+\dfrac{x^2r^2}{2}+\dfrac{x^3r^3}{6}\right)\right]\right\}$"], "1.10")
    Pp(f"Check: (1.1)–(1.8) reproduce the code's plotted curve at baseline to {R['check_vs_code_maxabs']:.1e} "
       f"(max abs difference over 361 points), and (1.10) matches a numerical derivative of (1.7) to "
       f"{R['check_hprime_maxabs']:.0e}.")
    # ---------------- section 2
    H("2. Extremum condition")
    Pp("Setting $dA/d\\ell=0$ (Chris's form), with $\\partial AF/\\partial\\ell=(\\pi/180)\\,\\partial AF/\\partial\\theta$ from (1.9):")
    E([r"$-\mathrm{Amp}_A\,\dfrac{\pi}{180}\cos\left(\dfrac{\pi\ell}{180}+\phi_A\right)AF+\left[1-\mathrm{Amp}_A\sin\left(\dfrac{\pi\ell}{180}+\phi_A\right)\right]\dfrac{\partial AF}{\partial\ell}=0$"], "2.1")
    Pp("Multiplying by $180\\,W\\bar h/\\pi$ and inserting (1.8)–(1.9), the normalisation $\\bar h$ cancels. The "
       "extremum angles $\\theta_k=\\pi\\ell_k/180$ ($k=1..4$) are the roots of")
    E([r"$\left[1-\mathrm{Amp}_A\sin\left(\theta_k+\phi_A\right)\right]\left[h\left(\theta_k+\dfrac{W}{2}\right)-h\left(\theta_k-\dfrac{W}{2}\right)\right]$",
       r"$=\ \mathrm{Amp}_A\cos\left(\theta_k+\phi_A\right)\int_{\theta_k-W/2}^{\theta_k+W/2}h(\theta^\prime)\,d\theta^\prime$"], "2.2")
    Pp("with $h$ from (1.7). Plain words: the slope of the smoothed curve is set by how different $h$ is at the two "
       "edges of the 90° window, and the sine factor tilts it in proportion to the area under $h$ inside the window. "
       "Because $h$ contains $\\mathrm{e}^{-x\\sqrt{1+2\\varepsilon\\cos\\theta+\\varepsilon^2}}$ and the integral has "
       "no elementary form, (2.2) is transcendental in $\\theta_k$: there is no exact closed form for $\\ell_k$. "
       "It is solved numerically (bracketed root on a 0.25° scan, then Brent to $10^{-13}$ rad). A closed form exists "
       "only as a series in small $\\varepsilon$ and small $\\mathrm{Amp}_A$ (Section 6).")
    Pp("Type of extremum, from the second derivative ($<0$ maximum, $>0$ minimum):")
    E([r"$\dfrac{\partial^2A}{\partial\theta^2}=\mathrm{Amp}_A\sin\left(\theta+\phi_A\right)AF-2\,\mathrm{Amp}_A\cos\left(\theta+\phi_A\right)\dfrac{\partial AF}{\partial\theta}$",
       r"$+\left[1-\mathrm{Amp}_A\sin\left(\theta+\phi_A\right)\right]\dfrac{\partial^2AF}{\partial\theta^2}$"], "2.3")
    # ---------------- section 3
    H("3. Extremum value")
    E([r"$A^*_k=A(\ell_k)=\left[1-\mathrm{Amp}_A\sin\left(\theta_k+\phi_A\right)\right]\dfrac{1}{W\bar h}\int_{\theta_k-W/2}^{\theta_k+W/2}h(\theta^\prime)\,d\theta^\prime$"], "3.1")
    Pp("Using (2.2) to replace the integral (valid wherever $\\cos(\\theta_k+\\phi_A)\\neq0$), the value needs only "
       "$h$ at the two window edges and $\\bar h$:")
    E([r"$A^*_k=\dfrac{\left[1-\mathrm{Amp}_A\sin\left(\theta_k+\phi_A\right)\right]^2\left[h\left(\theta_k+\frac{W}{2}\right)-h\left(\theta_k-\frac{W}{2}\right)\right]}{\mathrm{Amp}_A\cos\left(\theta_k+\phi_A\right)\,W\,\bar h}$"], "3.2")
    Pp(f"Check: (3.2) equals (3.1) at all four baseline extrema to "
       f"{max(abs(e['A_from_eq3_2']-e['A']) for e in b):.0e}. Baseline $\\bar h={R['hbar_base']:.5f}$ "
       f"($\\alpha=1$) and {R['hbar_base13']:.5f} ($\\alpha=1.3$). Because $\\partial A/\\partial\\theta=0$ at "
       "$\\theta_k$, the shift of the angle does not change the value at first order (envelope theorem):")
    E([r"$\dfrac{dA^*_k}{dp}=\left.\dfrac{\partial A}{\partial p}\right|_{\theta=\theta_k},\qquad p\in\{x,\ \varepsilon,\ \alpha,\ \kappa\}$"], "3.3")
    # ---------------- section 4
    H("4. Merge (fold) condition")
    Pp("A maximum and a minimum merge and vanish where the curve has a horizontal inflection: (2.2) and (2.3) vanish "
       "together. With one input $p$ free and the others at baseline, this is two equations in two unknowns "
       "$(\\theta_f,\\,p_f)$:")
    E([r"$\left.\dfrac{\partial A}{\partial\theta}\right|_{\theta_f,\,p_f}=0\quad\mathrm{[eq.\ 2.2]},\qquad\left.\dfrac{\partial^2A}{\partial\theta^2}\right|_{\theta_f,\,p_f}=0\quad\mathrm{[eq.\ 2.3]}$"], "4.1")
    T(["free input", "fold from (4.1)", r"$\ell_f$ [deg]", r"$A_f$", "sweep (bisection)"],
      [["D/L", f"{fx['value']:.5f}", f"{fx['l']:.2f}", f"{fx['A']:.5f}", "3.15674"],
       ["α", f"{fa['value']:.5f}", f"{fa['l']:.2f}", f"{fa['A']:.5f}", "0.30558"],
       ["κ", f"{fk['value']:.5f}", f"{fk['l']:.2f}", f"{fk['A']:.5f}", "outside 0.15–0.70"],
       ["d/D (α=1.3)", f"{fe['value']:.5f}", f"{fe['l']:.2f}", f"{fe['A']:.5f}", "outside 0.06–0.74"]],
      [1.4, 1.6, 1.2, 1.1, 2.2])
    Pp(f"The two systems close to residuals of order $10^{{-16}}$ and agree with the sweep values to $10^{{-4}}$ "
       f"(the sweep used the code's 1° sampled curve). The merging pair is always the maximum near 170° and the "
       f"minimum near 95°. Since only $\\alpha\\kappa$ enters, the $\\alpha$ fold sits at $\\alpha\\kappa="
       f"{fa['value']*K0:.4f}$, so at $\\alpha=1$ a $\\kappa$ fold would need $\\kappa={fa['value']*K0:.4f}$, below the "
       "swept $0.15$–$0.70$; that is why the $\\kappa$ sweep shows no merge (solving (4.1) directly in $\\kappa$ gives the same "
       f"{fk['value']:.4f}). In $d/D$ ($\\alpha=1.3$) no pair meets inside the swept 0.06–0.74; the maximum near "
       f"173° and the minimum near 94° drift together and merge at $d/D={fe['value']:.4f}$ ($\\ell_f={fe['l']:.1f}^\\circ$). "
       f"Just beyond, at $d/D={fe2['value']:.4f}$, a new max–min pair is born near $\\ell={fe2['l']:.0f}^\\circ$, "
       "because the sampled shell then nearly touches the CoE ($r\\to1-\\varepsilon$ at $\\theta=\\pi$) and the model "
       "becomes singular at $d/D=1$.")
    Pp("Near a fold the two merging extrema move as a square root, not linearly. Keeping the leading terms of "
       "$\\partial A/\\partial\\theta$ around $(\\theta_f,p_f)$, where $A_{\\theta\\theta}=0$: "
       "$A_{\\theta p}\\,\\delta p+\\frac{1}{2}A_{\\theta\\theta\\theta}\\,\\delta\\theta^2=0$, so")
    E([r"$\ell_\pm\approx\ell_f\pm\dfrac{180}{\pi}\sqrt{-\dfrac{2\,A_{\theta p}}{A_{\theta\theta\theta}}\,\left(p-p_f\right)}$"], "4.2")
    Pp(f"with $A_{{\\theta p}}=\\partial^2A/\\partial\\theta\\partial p$ and $A_{{\\theta\\theta\\theta}}=\\partial^3A/\\partial\\theta^3$ at the fold. "
       f"Numerically $\\ell_\\pm\\approx{fx['l']:.1f}^\\circ\\pm{fx['puiseux_deg']:.1f}^\\circ\\sqrt{{D/L-{fx['value']:.4f}}}$ and "
       f"$\\ell_\\pm\\approx{fa['l']:.1f}^\\circ\\pm{fa['puiseux_deg']:.1f}^\\circ\\sqrt{{\\alpha-{fa['value']:.4f}}}$ "
       "(upper sign: the maximum). This square-root behaviour is why straight Taylor series in $D/L$ and $\\alpha$ "
       "(Section 5) fail toward the low ends of those sweeps.")
    # ---------------- section 5
    H("5. Sensitivities around baseline (implicit-function derivatives)")
    Pp("Baseline: $x=5.9$, $\\varepsilon=0.1$, $\\alpha=1$, $\\kappa=0.230588$ (Figs A10–A12). The A13 hold "
       "($\\alpha=1.3$) is a second base point, used for the $d/D$ sweep. Diagonal Taylor forms (one input varied at a time):")
    E([r"$\ell_k\approx\ell_{k,0}+\sum_p\dfrac{\partial\ell_k}{\partial p}\left(p-p_0\right)+\dfrac{1}{2}\sum_p\dfrac{\partial^2\ell_k}{\partial p^2}\left(p-p_0\right)^2$",
       r"$A^*_k\approx A^*_{k,0}+\sum_p\dfrac{\partial A^*_k}{\partial p}\left(p-p_0\right)+\dfrac{1}{2}\sum_p\dfrac{\partial^2A^*_k}{\partial p^2}\left(p-p_0\right)^2$"], "5.1")
    Pp("Cross terms $\\partial^2/\\partial p\\,\\partial q$ are not included. With $F=\\partial A/\\partial\\theta$ "
       "from (2.1)–(2.2) and subscripts denoting partial derivatives at $(\\theta_k,p_0)$:")
    E([r"$\dfrac{\partial\theta_k}{\partial p}=-\dfrac{F_p}{F_\theta},\qquad\dfrac{\partial^2\theta_k}{\partial p^2}=-\dfrac{F_{pp}+2F_{\theta p}\,\theta_k^\prime+F_{\theta\theta}\,\theta_k^{\prime\,2}}{F_\theta}$",
       r"$\dfrac{\partial A^*_k}{\partial p}=A_p,\qquad\dfrac{\partial^2A^*_k}{\partial p^2}=A_{pp}+F_p\,\theta_k^\prime,\qquad\theta_k^\prime\equiv\dfrac{\partial\theta_k}{\partial p}$"], "5.2")
    Pp("$F_\\theta=\\partial^2A/\\partial\\theta^2$ is evaluated from (2.3) analytically; the $p$-partials by central "
       "differences of the exact model. Each first and second derivative was cross-checked against central "
       "differences of the exact roots (agreement to 4–5 significant figures). Because only $\\alpha\\kappa$ enters, "
       "$\\partial/\\partial\\alpha=\\kappa\\,\\partial/\\partial\\kappa$ at $\\alpha=1$, visible in the tables.")
    nm = lambda e: f"k={e['k']} {e['kind']}"
    T(["extremum k", r"$\ell_{k,0}$ (α=1)", r"$A^*_{k,0}$ (α=1)", r"$\ell_{k,0}$ (α=1.3)", r"$A^*_{k,0}$ (α=1.3)"],
      [[nm(e), f"{e['l']:.3f}°", f"{e['A']:.5f}", f"{e2['l']:.3f}°", f"{e2['A']:.5f}"] for e, e2 in zip(b, b13)],
      [1.5, 1.45, 1.5, 1.5, 1.55])
    cols = ["x = D/L", "ε = d/D", "α", "κ", "ε (α=1.3)"]
    wid = [1.35, 1.2, 1.2, 1.2, 1.2, 1.25]
    def tab(title, key, fmt):
        T(["extremum k"] + cols,
          [[nm(e)] + [fmt(e["sens"][q][key]) for q in KEYS] + [fmt(e2["sens"]["e"][key])] for e, e2 in zip(b, b13)], wid, title)
    tab("First order, $\\partial\\ell_k/\\partial p$ [degrees per unit of $p$]:", "dl", lambda v: f3(v, 3))
    tab("Second order, $\\partial^2\\ell_k/\\partial p^2$ [degrees per unit of $p^2$] (enters (5.1) with the factor 1/2):", "d2l", lambda v: f3(v, 2))
    tab("First order, $\\partial A^*_k/\\partial p$ [per unit of $p$]:", "dA", lambda v: f3(v, 5))
    tab("Second order, $\\partial^2A^*_k/\\partial p^2$ [per unit of $p^2$]:", "d2A", lambda v: f3(v, 4))
    acc = R["accuracy"]
    acc_title = ("Accuracy against the exact extrema of the code model on the 200-point sweeps of extrema-vs-inputs "
       "(largest $|\\Delta\\ell|$ over the whole range where the extremum exists; the $d/D$ rows use the "
       "$\\alpha=1.3$ base point; 'O(ε)' and 'ε→0' are the closed forms of Section 6):")
    SWN = {"x": "D/L 3–10", "e": "d/D .06–.74", "a": "α 0.2–2", "k": "κ .15–.70"}
    rows = []
    for q in KEYS:
        for kk in (1, 2, 3, 4):
            a_ = acc[f"{q}|k{kk}"]
            rows.append([f"{SWN[q]}, k={kk}", f"{a_['lin_l']:.2f}°", f"{a_['quad_l']:.2f}°", f"{a_['cOe_l']:.2f}°",
                         f"{a_['c0_l']:.2f}°", f"{a_['lin_A']:.4f}", f"{a_['quad_A']:.4f}"])
    T(["sweep, extremum", r"linear $\ell$", r"quadr. $\ell$", r"O(ε) $\ell$", r"ε→0 $\ell$", r"linear $A^*$", r"quadr. $A^*$"], rows,
      [1.7, 0.93, 0.93, 0.93, 0.93, 1.02, 1.02], acc_title)
    rng_title = ("Range of each input over which the error stays within 1° (linear / quadratic Taylor) and within 2° "
       "(O(ε) closed form):")
    rows = []
    fmtr = lambda r: "—" if r is None else f"{r[0]:.3g}–{r[1]:.3g}"
    for q in KEYS:
        for kk in (1, 2, 3, 4):
            a_ = acc[f"{q}|k{kk}"]
            rows.append([f"{SWN[q]}, k={kk}", fmtr(a_["lin_1deg"]), fmtr(a_["quad_1deg"]), fmtr(a_["cOe_2deg"])])
    T(["sweep, extremum", "linear ≤1°", "quadratic ≤1°", "O(ε) closed ≤2°"], rows, [1.9, 1.85, 1.85, 1.85], rng_title)
    def mx(q, key):
        v = [(acc[f"{q}|k{kk}"][key], kk) for kk in (1, 2, 3, 4)]; return max(v)
    lin = {q: mx(q, "lin_l") for q in KEYS}; qd = {q: mx(q, "quad_l") for q in KEYS}
    Pp("Reading: the Taylor forms are good only near baseline. Over the full swept ranges the linear form misses by "
       f"up to {lin['x'][0]:.1f}° (D/L, k={lin['x'][1]}), {lin['e'][0]:.1f}° (d/D, k={lin['e'][1]}), "
       f"{lin['a'][0]:.1f}° (α, k={lin['a'][1]}) and {lin['k'][0]:.1f}° (κ, k={lin['k'][1]}); the quadratic form "
       f"brings these to {qd['x'][0]:.1f}°, {qd['e'][0]:.1f}° and {qd['a'][0]:.1f}° for D/L, d/D and α, but "
       f"overshoots for κ ({qd['k'][0]:.1f}°), so it is not uniformly better. The D/L and α misses peak at the low "
       "end, where the square-root fold (4.2) takes over; the κ miss comes from curvature in $\\alpha\\kappa$ at the "
       "high end. The closed forms of Section 6 stay within 2.4° (O(ε)) and 3.3° (ε→0) over the whole κ range, "
       "and within 2° over the upper parts of the D/L and α ranges.")
    # ---------------- section 6
    H(r"6. Closed form to first order in $\varepsilon$ and $\mathrm{Amp}_A$")
    Pp("Expanding (1.6) in $\\varepsilon$ to first order (Taylor series of $\\vec u$ along $\\hat n$ about $\\vec D$) "
       "turns $h$ into four harmonics. With $\\mu_1\\equiv m^\\prime(1)/\\kappa=3x^3/N(x)$ and "
       "$\\mu_2\\equiv m^{\\prime\\prime}(1)/\\kappa=3x^3\\left[2-x(1-f_u)\\right]/N(x)$, expanded:")
    E([r"$h(\theta)\approx h_0+h_1\cos\theta+h_2\cos2\theta+h_3\cos3\theta+O(\varepsilon^2)$"], "6.1")
    E([r"$h_0=1-\dfrac{\alpha\kappa}{2}\left[\dfrac{3x^3}{N(x)}-1\right],\qquad h_2=-\dfrac{\alpha\kappa}{2}\left[\dfrac{3x^3}{N(x)}-3\right]$"], "6.2")
    E([r"$h_1=-\dfrac{9\,\alpha\varepsilon\kappa}{8}\left[1-\dfrac{x^3\left(1+x(1-f_u)\right)}{N(x)}\right]$",
       r"$h_3=-\dfrac{3\,\alpha\varepsilon\kappa}{8}\left[5-\dfrac{x^3\left(5+x(1-f_u)\right)}{N(x)}\right]$"], "6.3")
    Pp("The boxcar multiplies harmonic $j$ by $\\sin(jW/2)/(jW/2)$, and the mean is $h_0$ to this order, so")
    E([r"$AF(\theta)\approx1+a_1\cos\theta+a_2\cos2\theta+a_3\cos3\theta,\qquad a_j=\dfrac{h_j}{h_0}\,\dfrac{\sin\left(jW/2\right)}{jW/2}$"], "6.4")
    h0, hs, (a1, a2, a3) = R["harmonics_base"]
    Pp(f"Baseline ($\\alpha=1$): $h_0={h0:.5f}$, $a_1={a1:.5f}$, $a_2={a2:.5f}$, $a_3={a3:.5f}$; the quadrupole "
       "$a_2$ dominates. With $\\mathrm{Amp}_A=0$ the extrema of (6.4) are at $\\theta=0$, $\\theta=\\pi$ and at "
       "$\\pm\\theta_c$ with")
    E([r"$\cos\theta_c=\dfrac{-a_2+\sqrt{a_2^2-3\,a_3\left(a_1-3a_3\right)}}{6\,a_3}\qquad\left(\mathrm{the\ root\ with\ }\left|\cos\theta_c\right|\leq1\right)$"], "6.5")
    Pp("The sine factor then shifts each $\\theta_{k,0}\\in\\{0,\\ \\theta_c,\\ \\pi,\\ 2\\pi-\\theta_c\\}$ to first order in $\\mathrm{Amp}_A$:")
    E([r"$\ell_k\approx\dfrac{180}{\pi}\left[\theta_{k,0}+\dfrac{\mathrm{Amp}_A\cos\left(\theta_{k,0}+\phi_A\right)\,AF(\theta_{k,0})}{AF^{\prime\prime}(\theta_{k,0})}\right]$",
       r"$A^*_k\approx\left[1-\mathrm{Amp}_A\sin\left(\theta_{k,0}+\phi_A\right)\right]AF(\theta_{k,0})$",
       r"$AF^{\prime\prime}(\theta)=-a_1\cos\theta-4a_2\cos2\theta-9a_3\cos3\theta$"], "6.6")
    Pp("Simplest explicit form ($\\varepsilon\\to0$, so $a_1=a_3=0$ and $\\theta_c=\\pi/2$). The four extrema are then "
       "fully explicit in the inputs, with $a_2=\\dfrac{\\sin W}{W}\\,\\dfrac{\\alpha\\kappa\\left[3-3x^3/N(x)\\right]}{2-\\alpha\\kappa\\left[3x^3/N(x)-1\\right]}$:")
    E67 = ([r"$\ell_1\approx-\dfrac{180}{\pi}\,\dfrac{\mathrm{Amp}_A\cos\phi_A\,(1+a_2)}{4a_2}$",
       r"$\ell_3\approx180^\circ+\dfrac{180}{\pi}\,\dfrac{\mathrm{Amp}_A\cos\phi_A\,(1+a_2)}{4a_2}$",
       r"$\ell_2\approx90^\circ-\dfrac{180}{\pi}\,\dfrac{\mathrm{Amp}_A\sin\phi_A\,(1-a_2)}{4a_2}$",
       r"$\ell_4\approx270^\circ+\dfrac{180}{\pi}\,\dfrac{\mathrm{Amp}_A\sin\phi_A\,(1-a_2)}{4a_2}$"], "6.7")
    E([E67[0][0], E67[0][2], E67[0][1], E67[0][3]], "6.7")
    E([r"$A^*_1\approx(1-\mathrm{Amp}_A\sin\phi_A)(1+a_2),\quad A^*_2\approx(1-\mathrm{Amp}_A\cos\phi_A)(1-a_2)$",
       r"$A^*_3\approx(1+\mathrm{Amp}_A\sin\phi_A)(1+a_2),\quad A^*_4\approx(1+\mathrm{Amp}_A\cos\phi_A)(1-a_2)$"], "6.8")
    cb = R["closed_base"]
    rows = []
    for e, c1, c0 in zip(b, sorted(cb["Oe"], key=lambda d: d["l"]), sorted(cb["e0"], key=lambda d: d["l"])):
        rows.append([nm(e), f"{e['l']:.2f}° / {e['A']:.4f}", f"{c1['l']:.2f}° / {c1['A']:.4f}", f"{c0['l']:.2f}° / {c0['A']:.4f}"])
    T(["extremum k", "exact (2.2), (3.1)", "O(ε), (6.5)–(6.6)", "ε→0, (6.7)–(6.8)"], rows, [1.5, 2.0, 2.0, 2.0],
      "Baseline ($\\alpha=1$), exact vs closed forms ($\\ell_k$ / $A^*_k$):")
    Pp("Validity: the expansion assumes $\\varepsilon\\ll1$ (the $O(\\varepsilon^2)$ harmonics are dropped) and "
       "$|\\mathrm{Amp}_A\\,AF/AF^{\\prime\\prime}|\\ll1$ (the sine tilt is small next to the quadrupole curvature, "
       "i.e. $|\\mathrm{Amp}_A|\\ll4a_2$). It is within 2° of the exact extrema for D/L ≥ 4.9, α ≥ 0.73, the whole "
       "κ range 0.17–0.70, and d/D ≤ 0.32 (table above). It breaks down as $a_2$ shrinks toward $|\\mathrm{Amp}_A|/4=0.015$: at both folds "
       "$a_2\\approx0.025$ (it is 0.086 at baseline), and it also fails for large $d/D$.")
    return doc

def render(doc):
    pg = Pager()
    for kind, pl in doc:
        if kind == "h": pg.heading(pl)
        elif kind == "p": pg.para(pl)
        elif kind == "eq": pg.eq(*pl)
        elif kind == "t": pg.table(*pl)
    return pg.finish("extrema-equations")

def to_md(doc):
    out = ["<!-- generated by extrema-equations.py; same content as extrema-equations-p*.png -->", ""]
    for kind, pl in doc:
        if kind == "h": out += [f"## {pl}", ""]
        elif kind == "p": out += [pl, ""]
        elif kind == "eq":
            lines, lab = pl
            body = [ln.strip("$") for ln in lines]
            tag = f"\\tag{{{lab}}}" if lab else ""
            if len(body) == 1: out += [f"$$ {body[0]} {tag} $$", ""]
            else: out += ["$$ \\begin{aligned} " + " \\\\ ".join("&" + b for b in body) + f" \\end{{aligned}} {tag} $$", ""]
        elif kind == "t":
            hdr, rows, _, title = pl
            if title: out += [title, ""]
            out += ["| " + " | ".join(hdr) + " |", "|" + "---|" * len(hdr)] + ["| " + " | ".join(r) + " |" for r in rows] + [""]
    return "\n".join(out)

def main():
    plt.rcParams.update({"mathtext.fontset": "cm", "font.family": "serif", "font.serif": ["STIXGeneral", "DejaVu Serif"]})
    R = compute()
    (NOTES / "extrema-equations.json").write_text(json.dumps(R, indent=1, default=lambda o: float(o) if np.ndim(o) == 0 else list(o)))
    doc = build(R)
    paths = render(doc)
    if OVER:
        print("\n".join(OVER)); raise SystemExit("equations too wide")
    (NOTES / "extrema-equations.md").write_text(to_md(doc))
    for p in paths: print(p)
    print("folds", R["folds"])

if __name__ == "__main__" and "--onepage" not in __import__("sys").argv:
    main()


# ====================================================== one-page equation card
def render_onepage(R):
    OP = []
    fig = plt.figure(figsize=(8.5, 11), dpi=DPI); rr = fig.canvas.get_renderer()
    X0, GAP, TOP, BOT = 0.3, 0.22, 0.3, 0.25
    COLW = (8.5 - 2 * X0 - GAP) / 2
    FB, FE, FH = 9.3, 10, 10.5
    def meas(sx, fs, **kw):
        t = fig.text(0, 0, sx, fontsize=fs, **kw); bb = t.get_window_extent(rr); t.remove()
        return bb.width / DPI, bb.height / DPI
    class Col:
        def __init__(c, x0, w, y): c.x0, c.w, c.y = x0, w, y
        def put(c, sx, fs, dx=0.0, **kw):
            wd, ht = meas(sx, fs, **kw)
            if wd > c.w - dx + 1e-3: OP.append(f"wide {wd:.2f}>{c.w-dx:.2f}: {sx[:60]}")
            fig.text((c.x0 + dx) / 8.5, 1 - c.y / 11, sx, fontsize=fs, va="top", **kw)
            return wd, ht
        def head(c, sx):
            c.y += 0.07; _, ht = c.put(sx, FH, fontweight="bold"); c.y += ht + 0.03
            fig.add_artist(Line2D([c.x0 / 8.5, (c.x0 + c.w) / 8.5], [1 - c.y / 11] * 2, color="k", lw=0.6))
            c.y += 0.07
        def para(c, sx, fs=FB, indent=0.0):
            words = tok(sx); cur = ""; lines = []
            for wd_ in words:
                trial = (cur + " " + wd_).strip()
                if meas(trial, fs)[0] > c.w - indent and cur: lines.append(cur); cur = wd_
                else: cur = trial
            if cur: lines.append(cur)
            for ln in lines:
                _, ht = c.put(ln, fs, dx=indent); c.y += max(ht, fs / 72 * 1.2) + 0.035
            c.y += 0.04
        def eq(c, lines, lab=None, fs=FE):
            c.y += 0.03
            for i, ln in enumerate(lines):
                cont = i > 0 and re.match(r"^\$(=|\\times|\+|-|\\left\.)", ln) is not None
                dx = 0.08 + (0.3 if cont else 0)
                wd, ht = meas(ln, fs)
                labw = 0.3 if (lab and i == len(lines) - 1) else 0
                if wd > c.w - dx - labw: OP.append(f"eq wide ({lab}) {wd:.2f}>{c.w-dx-labw:.2f}: {ln[:60]}")
                fig.text((c.x0 + dx) / 8.5, 1 - c.y / 11, ln, fontsize=fs, va="top")
                if labw:
                    fig.text((c.x0 + c.w) / 8.5, 1 - (c.y + ht / 2) / 11, f"({lab})", fontsize=FB, va="center", ha="right")
                c.y += ht + 0.08
            c.y += 0.03
        def table(c, hdr, rows, widths, fs=FB):
            rh = fs / 72 * 1.65
            def line(lw): fig.add_artist(Line2D([c.x0 / 8.5, (c.x0 + sum(widths)) / 8.5], [1 - c.y / 11] * 2, color="k", lw=lw))
            if sum(widths) > c.w + 1e-6: OP.append(f"table wide {sum(widths)}>{c.w}")
            line(0.9)
            for ri, row in enumerate([hdr] + rows):
                xx = c.x0
                for ci, (cell, wd) in enumerate(zip(row, widths)):
                    t = fig.text((xx + (0.04 if ci == 0 else wd / 2)) / 8.5, 1 - (c.y + rh / 2) / 11, cell,
                                 fontsize=fs, va="center", ha="left" if ci == 0 else "center",
                                 fontweight="bold" if ri == 0 else "normal")
                    cw_ = t.get_window_extent(rr).width / DPI
                    if cw_ > wd - 0.03: OP.append(f"cell {cell!r} {cw_:.2f}>{wd}")
                    xx += wd
                c.y += rh
                if ri == 0: line(0.6)
            line(0.9); c.y += 0.08
    b = R["base"]; b13 = R["base13"]; fx, fa, fk, fe, fe2 = R["folds"]
    # title
    full = Col(X0, 8.5 - 2 * X0, TOP)
    full.put(r"Extrema of the model $A/A_{\rm all}(\ell)$: equation card", 15, fontweight="bold"); full.y += 0.32
    full.para("Model of preprint Figs A10–A13 (code: af_an_alpha, scale_kappa, model_AB_fixed). Nothing fitted: roots, "
              "implicit derivatives and series terms only. Full derivation: extrema-equations-p1…p10.", fs=9.5)
    fig.add_artist(Line2D([X0 / 8.5, (8.5 - X0) / 8.5], [1 - full.y / 11] * 2, color="k", lw=1.2))
    full.y += 0.08
    full.y += 0.02
    for sx in [r"$\ell$: Galactic longitude, $\theta=\pi\ell/180$ ($\theta=0$: away from the center of expansion, CoE)",
               r"$x=D/L$: CoE–observer distance $D$ / core scale length $L$; $\varepsilon=d/D$: radius $d$ of the source shell",
               r"$\alpha$: Gauss peculiar-velocity scale; $\kappa=GM(<D)/(D^3H^2)$: mass inside $D$",
               r"$\mathrm{Amp}_A=-0.06$, $\phi_A=0.35$ rad: sine tilt; $f_u=0.94$: floor share of $\rho(D)$; "
               r"$W=90.75^\circ$ boxcar; $\bar h$: sky mean of $h$; $h^\prime=dh/d\theta$; $G=H=D=1$"]:
        full.para(sx)
    ytop = full.y + 0.0
    L = Col(X0, COLW, ytop); Rc = Col(X0 + COLW + GAP, COLW, ytop)
    # ---- left column
    L.head("Model")
    L.eq([r"$A(\ell)=\left[1-\mathrm{Amp}_A\sin\left(\theta+\phi_A\right)\right]AF(\theta)$"], "1")
    L.eq([r"$AF(\theta)=\dfrac{1}{W\bar h}\int_{\theta-W/2}^{\theta+W/2}h(\theta^\prime)\,d\theta^\prime$"], "2")
    L.eq([r"$h(\theta)=1+\dfrac{\alpha\kappa}{\varepsilon}\cos\theta-\dfrac{\alpha\kappa\,(\cos\theta+\varepsilon)}{\varepsilon\,r^3N(x)}$",
          r"$\times\left\{6(1-f_u)\left[\mathrm{e}^{x}-\mathrm{e}^{x(1-r)}\left(1+xr+\dfrac{x^2r^2}{2}\right)\right]\right.$",
          r"$\left.+f_u\,x^3r^3\right\}$"], "3")
    L.eq([r"$r=\sqrt{1+2\varepsilon\cos\theta+\varepsilon^2}$",
          r"$N(x)=6(1-f_u)\left(\mathrm{e}^{x}-1-x-\dfrac{x^2}{2}\right)+f_ux^3$"], "4")
    L.para(r"$h=H_{\rm obs}/H$ along $\theta$ (Gauss pull $\vec u=-\alpha m(r)\vec r/r^3$ differenced across the "
           r"shell); $\alpha,\kappa$ enter only as $\alpha\kappa$.")
    L.head(r"Extremum condition $dA/d\ell=0$")
    L.eq([r"$\dfrac{\partial AF}{\partial\ell}=\dfrac{\pi}{180}\,\dfrac{h\left(\theta+\frac{W}{2}\right)-h\left(\theta-\frac{W}{2}\right)}{W\bar h}$"], "5")
    L.eq([r"$\dfrac{dA}{d\ell}=\dfrac{\pi}{180}\left\{-\mathrm{Amp}_A\cos\left(\theta+\phi_A\right)AF(\theta)\right.$",
          r"$\left.+\left[1-\mathrm{Amp}_A\sin\left(\theta+\phi_A\right)\right]\dfrac{h\left(\theta+\frac{W}{2}\right)-h\left(\theta-\frac{W}{2}\right)}{W\bar h}\right\}=0$"], "6")
    L.para(r"Transcendental in $\theta$ ($\mathrm{e}^{-x\,r(\theta)}$ inside the integral (2)): no closed form; "
           r"roots $\ell_k$, $k=1..4$, found numerically.")
    # ---- right column
    Rc.head("Extremum value")
    Rc.eq([r"$A^*_k=\dfrac{\left[1-\mathrm{Amp}_A\sin\left(\theta_k+\phi_A\right)\right]^2}{\mathrm{Amp}_A\cos\left(\theta_k+\phi_A\right)}\,\dfrac{h\left(\theta_k+\frac{W}{2}\right)-h\left(\theta_k-\frac{W}{2}\right)}{W\bar h}$"], "7")
    Rc.para(r"(= (1) at $\theta_k$, with the integral eliminated by (6)).")
    Rc.head("Merge condition")
    Rc.eq([r"$\dfrac{dA}{d\ell}=0\quad\mathrm{and}\quad\dfrac{d^2A}{d\ell^2}=0$"], "8")
    Rc.eq([r"$\dfrac{d^2A}{d\theta^2}=\mathrm{Amp}_A\sin\left(\theta+\phi_A\right)AF-2\,\mathrm{Amp}_A\cos\left(\theta+\phi_A\right)AF^\prime$",
           r"$+\left[1-\mathrm{Amp}_A\sin\left(\theta+\phi_A\right)\right]\dfrac{h^\prime\left(\theta+\frac{W}{2}\right)-h^\prime\left(\theta-\frac{W}{2}\right)}{W\bar h}$"], "9")
    Rc.table(["free input", "merge at", r"$\ell_f$", r"$A_f$"],
             [["D/L", f"{fx['value']:.3f}", f"{fx['l']:.1f}°", f"{fx['A']:.3f}"],
              ["α", f"{fa['value']:.3f}", f"{fa['l']:.1f}°", f"{fa['A']:.3f}"],
              ["κ", f"{fk['value']:.4f}", f"{fk['l']:.1f}°", f"{fk['A']:.3f}"],
              ["d/D (α=1.3)", f"{fe['value']:.3f}", f"{fe['l']:.1f}°", f"{fe['A']:.3f}"]],
             [1.2, 0.85, 0.8, 0.75])
    Rc.para(r"(9) is $(180/\pi)^2\,d^2A/d\ell^2$. The max near 170° and min near 95° merge (173°/94° for d/D). "
            "Only the D/L and α merges lie inside the swept ranges; the κ merge is the α merge at the same ακ = 0.0705.")
    Rc.head("Small-input limit")
    Rc.eq([r"$\alpha\kappa\to0:\ AF\to1,\ A\to1-\mathrm{Amp}_A\sin\left(\theta+\phi_A\right)$"], "10")
    Rc.eq([r"$\ell_{\rm peak}=90^\circ-\dfrac{180}{\pi}\phi_A=69.9^\circ,\quad A=1.06$",
           r"$\ell_{\rm trough}=270^\circ-\dfrac{180}{\pi}\phi_A=249.9^\circ,\quad A=0.94$"], "11")
    Rc.para(r"(for $\mathrm{Amp}_A<0$). Below the merges one max and one min remain and head here "
            r"(D/L = 3: 28.9°, 260.3°; α = 0.2: 40.2°, 259.2°).")
    # ---- bottom full width table
    full.y = max(L.y, Rc.y) + 0.05
    full.head(r"Sensitivities at baseline ($x=5.9,\ \varepsilon=0.1,\ \alpha=1,\ \kappa=0.230588$)")
    nm = lambda e: f"{e['k']} {e['kind']}"
    hdr = ["k", r"$\ell_{k,0}$", r"$A^*_{k,0}$", r"$\partial\ell/\partial x$", r"$\partial\ell/\partial\varepsilon$",
           r"$\partial\ell/\partial\alpha$", r"$\partial\ell/\partial\kappa$", r"$\partial A^*/\partial x$",
           r"$\partial A^*/\partial\varepsilon$", r"$\partial A^*/\partial\alpha$", r"$\partial A^*/\partial\kappa$"]
    rows = [[nm(e), f"{e['l']:.2f}°", f"{e['A']:.4f}"] + [f3(e["sens"][q]["dl"], 2) for q in KEYS]
            + [f3(e["sens"][q]["dA"], 4) for q in KEYS] for e in b]
    full.table(hdr, rows, [0.68, 0.7, 0.66] + [0.69] * 4 + [0.705] * 4)
    acc = R["accuracy"]
    def inter(q):
        rs = [acc[f"{q}|k{kk}"]["lin_1deg"] for kk in (1, 2, 3, 4)]
        return max(r[0] for r in rs), min(r[1] for r in rs)
    def worst(q): return max(acc[f"{q}|k{kk}"]["lin_l"] for kk in (1, 2, 3, 4))
    ix, ie, ia, ik = inter("x"), inter("e"), inter("a"), inter("k")
    full.para(r"$\partial\ell/\partial p$ in degrees per unit $p$. Implicit differentiation of (6): "
              r"$\partial\ell_k/\partial p=-\dfrac{\partial^2A/\partial\ell\,\partial p}{\partial^2A/\partial\ell^2}$, "
              r"$\partial A^*_k/\partial p=\partial A/\partial p$ at $\ell_k$; $\partial/\partial\alpha=\kappa\,\partial/\partial\kappa$.", fs=9.5)
    full.para("A13 base (α = 1.3): ℓ = " + ", ".join(f"{e['l']:.2f}°" for e in b13) + "; ∂ℓ/∂ε = "
              + ", ".join(f3(e["sens"]["e"]["dl"], 2) for e in b13) + ". Linear form within 1° (all four extrema) for "
              f"D/L {ix[0]:.2f}–{ix[1]:.2f}, d/D {ie[0]:.2f}–{ie[1]:.2f} (α=1.3), α {ia[0]:.2f}–{ia[1]:.2f}, κ {ik[0]:.3f}–{ik[1]:.3f}; "
              f"over the full swept ranges it misses by up to {worst('x'):.1f}°, {worst('e'):.1f}°, {worst('a'):.1f}°, {worst('k'):.1f}°.", fs=9.5)
    if full.y > 11 - BOT: OP.append(f"page overflow y={full.y:.2f}")
    p = NOTES / "extrema-equations-onepage.png"
    fig.savefig(p, dpi=DPI); plt.close(fig)
    return p, OP, (L.y, Rc.y, full.y)

if __name__ == "__main__" and "--onepage" in __import__("sys").argv:
    plt.rcParams.update({"mathtext.fontset": "cm", "font.family": "serif", "font.serif": ["STIXGeneral", "DejaVu Serif"]})
    R = json.loads((NOTES / "extrema-equations.json").read_text())
    p, OP, ys = render_onepage(R)
    print(p, ys); print("\n".join(OP))
