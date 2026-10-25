"""Solve d(ell)/d(input) = 0 for the four tilt-card angles; numeric roots at baseline.
Writes nobox_tilt_stat_roots.json (formula roots + list of exact checks)."""
import json, numpy as np, sympy as sp
from scipy.optimize import brentq
from nobox_tilt_stat_sym import *
from nobox_check import K0
NUM = {al: 1, G: 1, H: 1, f: 0.94, a: -0.06, ph: 0.35}
BV = {D: 1.0, L: 1/5.9, d: 0.1, M: K0}
NAME = {("max", 1): "lmax1", ("max", -1): "lmax2", ("min", 1): "lmin1", ("min", -1): "lmin2"}
def cond(kind, X, sv):
    e = (Fmax if kind == "max" else Fmin)(X).subs(s, sv)
    return sp.lambdify([D, L, d, M], full(e).subs(NUM), "numpy")
def phys(v):   # physical range: 0<d<D, 0<alpha G M<D^3 H^2, D/L>1
    return 0 < v[d] < v[D] and 0 < v[M] < v[D]**3 and v[D]/v[L] > 1
RANGE = {D: (0.3, 3.0), L: (0.02, 1.5), d: (1e-3, 0.999), M: (1e-4, 0.9999)}
QF = {sv: sp.lambdify([D, L, d, M], full(Q.subs(s, sv)).subs(NUM), "numpy") for sv in (1, -1)}
DL = {}
def dl_fn(nm, X):
    if (nm, X) not in DL:
        DL[(nm, X)] = sp.lambdify([D, L, d, M], full(pd(quantities()[nm], X)).subs(NUM), "numpy")
    return DL[(nm, X)]
def classify(nm, X, v):
    """'ok' | 'unphys' | 'edge' (alpha G M = D^3 H^2) | 'pole' (card denominator of delta = 0)"""
    if abs(v[M] - v[D]**3) < 1e-7: return "edge"
    if nm.startswith("lmax") and abs(QF[1 if nm == "lmax1" else -1](v[D], v[L], v[d], v[M])) < 1e-6: return "pole"
    return "ok" if phys(v) else "unphys"
def roots_in(fn, Y, extra_hi=None):
    lo, hi = RANGE[Y]
    if Y == D: lo, hi = 0.2, 3.0
    g = np.linspace(lo, hi, 4001); out = []
    def F(y):
        v = dict(BV); v[Y] = y; return fn(v[D], v[L], v[d], v[M])
    ys = np.array([F(y) for y in g])
    for i in np.where(np.sign(ys[:-1]) != np.sign(ys[1:]))[0]:
        r = brentq(F, g[i], g[i+1], xtol=1e-13)
        v = dict(BV); v[Y] = r; out.append((r, v))
    return out
res = {}; checks = {}
for kind in ("max", "min"):
    for sv in (1, -1):
        nm = NAME[(kind, sv)]
        for X in INPUTS:
            fn = cond(kind, X, sv); key = f"{nm}_{X}"
            res[key] = {}
            for Y in INPUTS:
                lst = []
                for r, v in roots_in(fn, Y):
                    c = classify(nm, X, v); dv = float(dl_fn(nm, X)(v[D], v[L], v[d], v[M]))
                    lst.append((r, c, dv))
                res[key][str(Y)] = lst
            res[key]["base_value"] = float(fn(*[BV[k] for k in (D, L, d, M)]))
            print(key, res[key])
            for Y in INPUTS:
                for r, c, dv in res[key][str(Y)]:
                    if c == "ok": checks[f"{key}__{Y}"] = dict(name=nm, X=str(X), Y=str(Y), y0=r)
# explicit closed-form values at baseline (for the pages)
S0 = float(Sfun.subs(NUM).subs(BV)); w0 = 0.06; D0, L0 = 1.0, 1/5.9
kap_maxd = lambda x: (lambda Sx: Sx*(w0 + 4/x - 4*Sx/x)/(w0 + 8/x - 12*Sx/x + 4*Sx**2/x))(
    (1-.94)*(6/x**3*np.exp(x) - 6/x**3 - 6/x**2 - 3/x) + .94)
xs = np.linspace(1.01, 20, 4000); kv = kap_maxd(xs)
info = dict(S0=S0, kap_maxd_base=float(kap_maxd(5.9)), kap_maxd_min=float(kv[kv > 0].min()) if (kv > 0).any() else None,
            kap_maxd_minx=float(xs[kv > 0][np.argmin(kv[kv > 0])]) if (kv > 0).any() else None,
            d_line=float(4/3*np.tan(0.35)), lmin2_on_line=float(270 - np.degrees(np.tan(0.35))))
ac = -0.06*np.cos(0.35)
info["kap_min2_d"] = float(-ac*S0/(-ac*S0 + 6*(S0 - 1)))
info["M_min2_d"] = info["kap_min2_d"]
print(info)
checks["common_min2"] = dict(name="lmin2", X="d", Y="M", y0=info["M_min2_d"], over={"d": info["d_line"]})
json.dump(dict(roots=res, checks=checks, info=info), open("nobox_tilt_stat_roots.json", "w"), indent=1)
